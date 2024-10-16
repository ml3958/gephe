#!/usr/bin/env python

# Author: Menghan Liu
# Date: August 13, 2021
# produce two phylogenetic matrix from alignments


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

# def merge_pp(path):
#
#     all_files = glob.glob(path+'/*.csv')
#     df_from_each_file = (pd.read_csv(f) for f in all_files)
#     df_merged  = pd.concat(df_from_each_file,ignore_index=False,axis=1)
#     return(df_merged)

def pog_pp(file_pp, file_pog, cutoff):

    pp = pd.read_csv(file_pp,index_col=['query'])
    pog = pickle.load(open(file_pog,'rb'))
    protein_pog_dic = pog[1]
    protein_pog_dic['query'] = protein_pog_dic.protein
    protein_pog_dic.set_index('query')

    pog = pog[0]
    pp_pog = pd.DataFrame.merge(pp,protein_pog_dic, on='query')
    pp_pog = pp_pog.groupby(['pog']).sum()
    pp_pog_prop = pd.DataFrame(pp_pog.values/[[len(i)] for i in pog.values()],columns=pp_pog.columns,index = pp_pog.index)
    pp_pog_binary = (pp_pog_prop>float(cutoff)).astype(int)

    return(pp_pog_prop,pp_pog_binary)



if __name__ == '__main__':

    # merge protein-level pp to pog-level

    parser = argparse.ArgumentParser(description='collapse all single-sample pog pp')
    parser.add_argument('dir_pog', help='pog folder')
    parser.add_argument('-prefix_protein', help='${PREFIX_PROTEIN}')
    parser.add_argument('-prefix_pog', help='${PREIFX_POG}')
    parser.add_argument('-c', type=float, help='cutoff for binary phylogenetic profile')

    args = parser.parse_args()
    DIR_POG = args.dir_pog
    PREFIX_PROTEIN = args.prefix_protein
    PREFIX_POG = args.prefix_pog
    CUTOFF = args.c

    file_pp = DIR_POG + '/' + PREFIX_PROTEIN + '.pp'
    file_pog = DIR_POG + '/' + PREFIX_POG + '.mcl.pickle'
    fileout_pp_pog = DIR_POG + '/' + PREFIX_POG + '.pp'
    fileout_pp_pog_binary = DIR_POG + '/' + PREFIX_POG + '_binary' + str(CUTOFF) + '.pp'

    if os.path.exists(fileout_pp_pog):
        print("<" + fileout_pp_pog + "> exists, skipping...")

    else:
        print("processing: " + file_pp)
        pp_pog,pp_pog_binary = pog_pp(file_pp, file_pog, CUTOFF)

        print("Saving to: " + fileout_pp_pog)
        pp_pog.to_csv(fileout_pp_pog)
        print("Saving to: " + fileout_pp_pog_binary)
        pp_pog_binary.to_csv(fileout_pp_pog_binary)
