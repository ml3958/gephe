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


if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('dir_alignment',help='alignment folder')
    parser.add_argument('dir_association',help='association folder')
    parser.add_argument('dir_pog',help='pog folder')
    parser.add_argument('-p', dest = 'prefix', help='prefix')
    parser.add_argument('-percentile',  dest= 'percentile', help='perentile for selecting top proteins')


    args = parser.parse_args()
    DIR_ASSOCIATION = args.dir_association
    DIR_POG = args.dir_pog
    PERCENTILE = float(args.percentile)
    PREFIX = args.prefix

    output= DIR_POG + '/' + PREFIX + '.csv'

    if os.path.exists(output):
        print("<" + output + "> exists, skipping...")

    else:
        # import results with mi zscore
        result_combine = pd.read_pickle(DIR_ASSOCIATION+'/summary/mi_zscore.pickle')
        result_combine = result_combine[~result_combine.mi_z.isna()]
        # select top proteins
        print('saving result to ' + output)
        # print(np.quantile(result_combine.mi_z, float(PERCENTILE)))
        result_combine[result_combine.mi_z >= np.quantile(result_combine.mi_z,1-1*PERCENTILE*0.01)].to_csv(output)
