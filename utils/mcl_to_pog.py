#!/usr/bin/env python
import argparse
from Bio import SeqIO
import re, plotnine, sklearn,json, itertools, os,  pickle, sys,glob
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
from dplython import *
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import AgglomerativeClustering


def mcl_to_pog(mcl_file, dic_file):
    dic = pd.read_csv(dic_file,sep="\t",header=None,)
    dic = dic[1].to_dict()
    with open(mcl_file) as f:
        text = f.read().split('begin\n')[1].split('$')
        text = list(map(lambda x: x.replace('\n','').split(), text))[:-1]
        text = {int(i[0]):[dic[int(x)] for x in i[1:]] for i in text}
        f.close()
    protein_pog_dic = pd.DataFrame({'protein':list(itertools.chain.from_iterable(list(text.values()))),
                                    'pog' : np.repeat(list(text.keys()),[len(i) for i in text.values()])})
    protein_pog_dic.set_index('protein')
    return(text,protein_pog_dic)

if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('file_mcl', help='mcl file')
    parser.add_argument('file_dic', help='dic file')
    args = parser.parse_args()

    pickle.dump(mcl_to_pog(args.file_mcl, args.file_dic), open(args.file_mcl+'.pickle','wb'))


