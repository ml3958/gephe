#!/usr/bin/env python


# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: June 25, 2021

import argparse
from Bio import SeqIO
import re, plotnine, sklearn,json, itertools, os,  pickle, sys,glob
import pandas as pd

from datetime import datetime
from statsmodels import stats

from scipy.stats import zscore,pearsonr,spearmanr

import numpy as np
#from numba import jit,njit,cuda
import multiprocessing as mp
from itertools import islice


from plotnine import *
from matplotlib import gridspec
import seaborn
import matplotlib.pyplot as plt
from IPython import display

from dfply import *
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import AgglomerativeClustering


import scipy.cluster.hierarchy as shc
import dfply
import tqdm


def mi_null(metadata, times=1000):

    n = len(metadata)
    mat = pd.DataFrame([[1]*i + [0]*(n-i) for i in range(1,n+1)])

    mi_result = np.zeros((times,n)) # mat to store mutual information

    for i in range(times):
        mat_shuffled = mat.copy()
        mat_shuffled = mat_shuffled.apply(lambda x: pd.Series(np.random.choice(x,len(x),replace=False)),1)

        r_m = mat_shuffled.apply(lambda x: sklearn.metrics.mutual_info_score(x,metadata), axis=1).to_list()
        r_s = mat_shuffled.apply(lambda x: spearmanr(x,metadata)[0], axis=1)

        mi_result[i] =  np.sign(r_s) * r_m

    mu = mi_result.mean(axis=0)
    sd = mi_result.std(axis=0)

    return({i+1:[mu[i],sd[i]] for i in range(n)})


def mi_zscore(files,background):
    proteins = []
    mi_combine=pd.DataFrame()
    for f in files:
        mi = pd.read_pickle(open(f,'rb'))
        mi_combine = mi_combine.append(mi)
    mi_combine['mi_z'] = (mi_combine.mi - [background[i][0] for i in mi_combine.freq])/ [background[i][1] for i in mi_combine.freq]
    return(mi_combine)

if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('dir_association',help='path to association folder')
    parser.add_argument('file_metadata',help='metadata')
    parser.add_argument('metadata_colname',help='colname for phenotype')
    parser.add_argument('diamond_k',help='${ALIGNMENT_MAX}')
    parser.add_argument('-b','--extra_breaks', dest = 'extra_breaks' ,nargs='+', help='<Required> Set flag')

    args = parser.parse_args()
    dir_association = args.dir_association
    metadata_colname = args.metadata_colname
    diamond_k = args.diamond_k
    extra_breaks = args.extra_breaks


    metadata = pd.read_csv(args.file_metadata,sep="\t")

    # compute null distribution
    if os.path.exists(dir_association+'/summary/null.pickle'):
        print('loading percomputed null distribution:' + dir_association+'/summary/null.pickle')
        background = pickle.load(open(dir_association+'/summary/null.pickle','rb'))
    else:
        print('calculating null distribution:' + dir_association+'/summary/null.pickle')
        background = mi_null(metadata[metadata_colname])
        pickle.dump(background,open(dir_association+'/summary/null.pickle','wb'))

    # compute mi Z score
    result_combine = mi_zscore(glob.glob(dir_association+'/*mi.pickle'),background)
    pickle.dump(result_combine,open(dir_association+'/summary/mi_zscore.pickle','wb'))

    # print(extra_breaks)
    result_combine = result_combine[~result_combine.mi_z.isna()]

    # plot mi
    g=(ggplot(result_combine,aes('mi')) +
        geom_histogram(data=result_combine, bins=100, fill="red", alpha=0.5)+
        theme(figure_size=(5,5))+
        theme_bw())
    g.save(dir_association+'/summary/distribution_mi.png',dpi=300)

    breaks = [0.99,0.995,0.996,0.997,0.998,0.999]
    if extra_breaks:
        breaks = breaks + extra_breaks
    breaks = [float(i) for i in breaks]
    breaks = pd.DataFrame({"breaks": np.quantile(result_combine.mi_z,breaks),
                          "y": np.Inf-10,
                          "label": [(str(round((1-i)*100,6)) +'%') for i in breaks]})

    g=(ggplot(result_combine,aes('mi_z')) +
        geom_vline(xintercept= breaks.breaks,color="red", linetype='dotted',size=0.5) +
        geom_text(aes(x='breaks',y='y',label='label'),data = breaks,va="top",size=10) +
        geom_histogram(data=result_combine, bins=100, fill="red", alpha=0.5)+
          labs(x="Mutual information zscore",y='count') +
          theme(figure_size=(8,5))+
          theme_bw())
    g.save(dir_association+'/summary/mi_z.png',dpi=300)

    # plot mi_zscore
    g=(ggplot(result_combine[result_combine.mi_z > np.quantile(result_combine.mi_z,0.9)],
              aes('mi_z')) +
        geom_vline(xintercept= breaks.breaks,color="red", linetype='dotted',size=0.5) +
        geom_text(aes(x='breaks',y='y',label='label'),data = breaks,va="top",size=10) +
        geom_histogram(bins=100, fill="red", alpha=0.5)+
       labs(x="Mutual information zscore",y='count') +
      theme(figure_size=(8,5))+
  theme_bw())
    g.save(dir_association+'/summary/mi_z2.png',dpi=300)

    # plot HSP count relative to diamond_k
    HSP_count = pd.DataFrame()
    for f in tqdm.tqdm(glob.glob(dir_association + '/*HSP_count.pickle')):
        HSP_count = pd.concat((HSP_count, pd.read_pickle(f)))
    fig, ax = plt.subplots()
    HSP_count.hist(ax=ax)
    fig.savefig(dir_association + '/summary/hsp_count.png')
    plt.axvline(diamond_k, color='r') # vertical
    fig.savefig(dir_association + '/summary/hsp_count_with_k.png')
    plt.xlim([diamond_k-10,diamond_k+10]) # vertical
    plt.ylim([0,100]) # vertical
    fig.savefig(dir_association + '/summary/hsp_count_with_k_zoom.png')
