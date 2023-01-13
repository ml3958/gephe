#!/usr/bin/env python


# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: July 20, 2021

import argparse
from Bio import SeqIO
import re, plotnine, sklearn,json, itertools, os, pickle, sys,glob
# import skbio
import pandas as pd

from datetime import datetime
from statsmodels import stats

from scipy.stats import zscore,pearsonr,spearmanr

import numpy as np
from numba import jit,njit,cuda
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


def read_diamond_output(diamond_file):
    diamond = pd.read_table(diamond_file,
                            header=None,sep='\t',
                            names=['qseqid','sseqid','pident','length','mismatch','gapopen','qlen','qstart','qend','slen','start','send','evalue','bitscore'] # column already exists for customized diamond output format
                            )
    return(diamond)

def diamond_to_pp(diamond,genomes):
    pp = pd.DataFrame({'query':diamond.qseqid,
                       'genome':diamond.sseqid.map(lambda x: x.split('|')[0]),
                       'value':1})
    pp = pp.drop_duplicates()
    pp = pp.pivot_table(index='query',columns='genome',values='value',fill_value=0)

    # add space for genomes that are absent from the alignment
    genomes_absent = np.setdiff1d(genomes, pp.columns)
    pp_absent = pd.DataFrame(np.zeros((pp.shape[0],len(genomes_absent))), index=pp.index, columns = genomes_absent)
    pp2 = pd.concat((pp,pp_absent),axis=1)

    pp2 = pp2.loc[:,genomes]
    print(pp2.head())
    print(pp2.shape)
    return(pp2)

def pp_metadata_association(pp,metadata):
    r_s =pp.apply(lambda x: spearmanr(x,metadata) ,axis = 1)
    r_m = pp.apply(lambda x: sklearn.metrics.mutual_info_score(x,metadata) , axis = 1)
    result = pd.concat([pp.sum(axis=1),
                        r_s.map(lambda x: x[0]),
                        r_s.map(lambda x: x[1]),
                        r_m,
                        np.sign(r_s.map(lambda x: x[0])) * r_m],
                       axis=1)
    result.columns = ['freq','spearman_r','spearman_p','mi_orig','mi']
    result.sort_values('mi',ascending=False,inplace=True)
    return(result)

if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('file_diamond', help='diamond output of pickle')
    parser.add_argument('file_metadata', help='metadata')
    parser.add_argument('output_dir', help='Output directory')
    args = parser.parse_args()

    
    file_diamond = args.file_diamond
    file_metadata = args.file_metadata
    prefix = os.path.basename(file_diamond).replace('.diamond.out.pickle','')
    fileout_pp = args.output_dir + '/' + prefix+ '_pp.pickle'
    fileout_hsp_count = args.output_dir + '/' + prefix+ '_HSP_count.pickle'
    
    metadata = pd.read_table(file_metadata)
    genomes = metadata.taxon_oid
    
    if not os.path.exists(fileout_pp):

        #diamond = pd.read_pickle(open(file_diamond,'rb'))
        diamond=read_diamond_output(file_diamond)
        diamond[['qseqid']].value_counts().to_pickle(fileout_hsp_count)

        pp = diamond_to_pp(diamond,genomes)
        pickle.dump(pp, open(fileout_pp,'wb'))
