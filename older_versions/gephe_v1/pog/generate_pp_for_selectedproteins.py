#!/usr/bin/env python

# Author: Menghan Liu
# Date: August 13, 2021
# produce two phylogenetic matrix from alignments


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


def pp(file_protein, dir_association):

    # proteins list
    proteins = pd.read_csv(file_protein).loc[:,'query'].to_list()

    # protein-level pp
    pp = pd.DataFrame([])
    for f in glob.glob(dir_association + '/*_pp.pickle'):
        p = pd.read_pickle(f)
        pp = pd.concat([pp,p[p.index.isin(proteins)]])

    return(pp)



if __name__ == '__main__':
    # gather the phylogenetic profile for selected proteins
    parser = argparse.ArgumentParser(description='collapse all single-sample pog pp')
    parser.add_argument('dir_association', help='association folder (contain sample-level protein)')
    parser.add_argument('dir_pog', help='pog folder')
    parser.add_argument('-prefix_protein', help='${PREFIX_PROTEIN}')

    args = parser.parse_args()


    DIR_ASSOCIATION = args.dir_association
    DIR_POG = args.dir_pog
    PREFIX_PROTEIN = args.prefix_protein

    FILE_PROTEIN = DIR_POG + '/' + PREFIX_PROTEIN + '.csv'
    FILEOUT_PP = DIR_POG + '/' + PREFIX_PROTEIN + '.pp'

    print("processing: " + FILE_PROTEIN)
    pp = pp(FILE_PROTEIN,DIR_ASSOCIATION)

    print("Saving protein pp to: " + FILEOUT_PP)
    pp.to_csv(FILEOUT_PP)
