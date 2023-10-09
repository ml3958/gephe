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


def mi_null_shuffle_gen(metadata, times=1000):

    n = len(metadata)
    mat = pd.DataFrame([[1]*i + [0]*(n-i) for i in range(1,n)])
#     print(n, mat.shape)

    mi_result = pd.DataFrame(np.zeros((n-1,times))) # mat to store mutual information
    sp_result = pd.DataFrame(np.zeros((n-1,times))) # mat to store spearman correlation
#     print(mi_result.shape)
    mi_result.index = sp_result.index = ['freq'+str(i) for i in range(1,n)]


    for i in tqdm.tqdm(range(times)):
        mat_shuffled = mat.copy()
        mat_shuffled = mat_shuffled.apply(lambda x: pd.Series(np.random.choice(x,len(x),replace=False)),1)

        r_m = mat_shuffled.apply(lambda x: sklearn.metrics.mutual_info_score(x,metadata), axis=1).to_list()
        r_s = mat_shuffled.apply(lambda x: spearmanr(x,metadata)[0], axis=1).to_list()
#         print(len(r_m), len(r_s))

        mi_result[i] =  r_m
        sp_result[i] =  r_s

    return({'mi':pd.DataFrame(mi_result),'spearman':pd.DataFrame(sp_result)})


def mi_null_shuffle_phe(metadata, times=1000):

    n = len(metadata)
    mat = pd.DataFrame([[1]*i + [0]*(n-i) for i in range(1,n)])
#     print(n, mat.shape)

    mi_result = pd.DataFrame(np.zeros((n-1,times))) # mat to store mutual information
    sp_result = pd.DataFrame(np.zeros((n-1,times))) # mat to store spearman correlation
    phe_shuffled = pd.DataFrame(np.zeros((len(metadata),times)))
    mi_result.index = sp_result.index = ['freq'+str(i) for i in range(1,n)]

    for i in tqdm.tqdm(range(times)):
        metadata_shuffled = np.random.choice(metadata,len(metadata),replace=False)

        r_m = mat.apply(lambda x: sklearn.metrics.mutual_info_score(x,metadata_shuffled), axis=1).to_list()
        r_s = mat.apply(lambda x: spearmanr(x,metadata_shuffled)[0], axis=1).to_list()
#         print(len(r_m), len(r_s))

        mi_result[i] =  r_m
        sp_result[i] =  r_s
        phe_shuffled[i] = metadata_shuffled

    return({'mi':pd.DataFrame(mi_result),'spearman':pd.DataFrame(sp_result),'shuffled_phenotype':phe_shuffled})


if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('file_metadata',help='metadata')
    parser.add_argument('metadata_colname',help='colname for phenotype')
    parser.add_argument('prefix',help='output_dir')
    parser.add_argument('dir_out',help='output_dir')
    parser.add_argument('-s','--shuffle_type', dest = 'shuffle_type' , default = 'phenotype')
    parser.add_argument('-t','--times_of_simulation', dest = 'times' , default = 1000, type=int)


    args = parser.parse_args()
    file_metadata = args.file_metadata
    metadata_colname = args.metadata_colname
    dir_out = args.dir_out
    prefix = args.prefix
    times = args.times
    shuffle_type = args.shuffle_type


    metadata = pd.read_csv(args.file_metadata,sep="\t")
    # print(metadata)

    if shuffle_type=='phenotype':
        background = mi_null_shuffle_phe(metadata[metadata_colname],times)
    else:
        background = mi_null_shuffle_gen(metadata[metadata_colname],times)

    pickle.dump(background,open(dir_out+prefix+'.pickle','wb'))
