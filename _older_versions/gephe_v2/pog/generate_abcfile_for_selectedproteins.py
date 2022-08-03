
#!/usr/bin/env python


# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: June 25, 2021

import argparse
from Bio import SeqIO
import re, plotnine, sklearn,json, itertools, os, skbio, pickle, sys,glob
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

def get_top_proteins(files, cutoff):
    proteins = []
    mi_combine=pd.DataFrame()
    background_combine = []
    for f in files:
        mi, background = pd.read_pickle(open(f,'rb'))
        mi_combine = mi_combine.append(mi)
        background_combine = np.append(background_combine,background)

    mean = mi_combine.mi.mean()
    std = mi_combine.mi.std()
    mi_combine['mi_z'] = mi_combine.mi.map(lambda x: (x-mean)/std)
    mi_combine['mi_z_log10'] = np.log10(mi_combine['mi_z'])

    background_combine = pd.Dataframe({'mi':background_combine,'mi_z':(background_combine - background_combine.mean())/background_combine.std()})
    background_combine['mi_z_log10'] = np.log10(background_combine['mi_z'])

    return(mi_combine,background_combine)

def get_alignment(protein,files):
    result = pd.DataFrame()
    for f in files:
        # print(f)
        diamond = pd.read_pickle(open(f,'rb'))
        result = result.append(diamond[(diamond.qseqid.isin(proteins) & diamond.sseqid.isin(proteins))])
    return(result)

def get_pp(proteins,files):
    pp_combine = pd.DataFrame()
    for f in files:
        pp = pd.read_pickle(open(f,'rb'))
        pp_combine = pp_combine.append(pp[pp.index.isin(proteins)])
    return(pp_combine)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='produce pairwise alignment (.abc) for a provided list of proteins')
    parser.add_argument('dir_alignment',help='alignment folder')
    parser.add_argument('dir_association',help='association folder')
    parser.add_argument('dir_pog',help='pog folder')
    parser.add_argument('-p', dest = 'prefix', help='prefix')
    args = parser.parse_args()


    DIR_ALIGNMENT = args.dir_alignment
    DIR_ASSOCIATION  =args.dir_association
    DIR_POG =args.dir_pog
    PREFIX = args.prefix

    input,output = DIR_POG + '/' + PREFIX + '.csv', DIR_POG + '/' + PREFIX + '.abc'

    if os.path.exists(output):
        print("<" + output + "> exists, skipping...")

    else:
        print('Processing :', input)
        print('Saving to :', output)

        proteins = pd.read_csv(input)['query']
        alignment = get_alignment(proteins,glob.glob(DIR_ALIGNMENT+'/*diamond.out.pickle'))
        alignment[['qseqid','sseqid','evalue']].to_csv(output,sep="\t",header=False,index=False)
