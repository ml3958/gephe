#!/usr/bin/env python


# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: June 6, 2023

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

def diamond_select_max_bitscore_per_sseqid(diamond):

    # Find the maximum value in for each sseqid
    max_values = diamond.loc[diamond.groupby('sseqid')['bitscore'].idxmax(),['sseqid','bitscore']]
    print(max_values.head())

    # # Select rows with maximum value for each group
    diamond_sseq_max_value = max_values.merge(diamond,  how='left')

    return(diamond_sseq_max_value[diamond.columns])


if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('file_diamond', help='Process a diamond output to keep only HSP with highest bitscore for each sseqid')
    args = parser.parse_args()
    parser.add_argument('file_diamond', help='diamond output of pickle')
    parser.add_argument('output_dir', help='Output directory')
    args = parser.parse_args()


    file_diamond = args.file_diamond
    output_dir = args.output_dir
    file_out = output_dir + '/' + os.path.basename + '.maxbitscore.txt'

    diamond=read_diamond_output(file_diamond)

    diamond_select = diamond_select_max_bitscore_per_sseqid(diamond)

    diamond_select.to_csv(file_out,sep='\t',index=False)
