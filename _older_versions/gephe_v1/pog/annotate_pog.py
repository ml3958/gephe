#!/usr/bin/env python

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
from dplython import *
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import AgglomerativeClustering


import scipy.cluster.hierarchy as shc
import dfply

def read_diamond_output(diamond_file):
    diamond = pd.read_table(diamond_file,
                            header=None,sep='\t',
                            names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','start','send','evalue','bitscore'])
    return(diamond)

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


def pog_annotate_summary(protein_annot, annot_N = 2, column_for_sorting ='N'):
    protein_annot.keywords.fillna('',inplace=True)
    # protein_annot =
    pog_N = DplyFrame(protein_annot)>> group_by(X.pog) >> summarize(N = X.pog.count())
    pog_keywords = DplyFrame(protein_annot) >>  group_by(X.pog,X.keywords) >> summarize(keyword_N = X.pog.count()) >> \
        sift(X.keyword_N>annot_N) >> mutate(keywords = X.keywords + ' (n='+X.keyword_N.astype(str) +')') >> \
        ungroup() >> arrange(-X.keyword_N)>> group_by(X.pog) >>summarize(keywords = X.keywords.unique()) >> \
        mutate(keywords = X.keywords.map(lambda x: dplyr_paste(x,'; ')))

    result = full_join(pog_N,pog_keywords).sort_values(column_for_sorting,ascending=False)
    return(result)


def dplyr_paste(x,sep=','):
    if type(x)==str: return x
    else: return sep.join(x)




if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='pog functional annotation')
    parser.add_argument('dir_alignment',help='alignment folder')
    parser.add_argument('dir_pog',help='pog folder')
    parser.add_argument('prefix_protein',help='${PREFIX_PROTEIN}')
    parser.add_argument('prefix_pog',help='${PREFIX_POG}')
    parser.add_argument('-f',dest='pog_protein_fraction',type=float,help='${PREFIX_POG}')


    args = parser.parse_args()
    DIR_ALIGNMENT = args.dir_alignment
    DIR_POG = args.dir_pog
    PREFIX_PROTEIN = args.prefix_protein
    PREFIX_POG = args.prefix_pog
    REP_PROTEIN_FRAC = args.pog_protein_fraction

    # file_mcl = DIR_POG + '/' + PREFIX_POG + '.mcl'
    # file_protein_dic = DIR_POG + '/' + PREFIX + '.dic'

    # input
    FILE_POG = DIR_POG + PREFIX_POG + '.mcl'
    FILE_PROTEIN_DIC = DIR_POG + '/' + PREFIX_PROTEIN + '.dic'
    FAA = DIR_ALIGNMENT + '/input.fasta'

    # output
    FILE_REPSEQ = DIR_POG + PREFIX_POG  + '.repseq'
    FILE_REPSEQ_FAA = DIR_POG + PREFIX_POG + '.repseq.faa'

    FILE_REPSEQ_ALIGNMENT_META = DIR_POG + PREFIX_POG + '.repseq.faa_meta.diamond.out'
    # FILE_REPSEQ_ALIGNMENT_META_ANNOTATION = DIR_POG + PREFIX_POG + '.annot_meta'

    # FILE_REPSEQ_ALIGNMENT_ECOLI = DIR_POG + PREFIX_POG + '.repseq.faa_ecoli.diamond.out'
    # FILE_REPSEQ_ALIGNMENT_ECOLI_ANNOTATION = DIR_POG + PREFIX_POG + '.annot_ecoli'

    out_full = DIR_POG + '/' + PREFIX_PROTEIN + '.annot'
    out_summary = DIR_POG + '/' + PREFIX_POG + '.annot'

    print('Processing :' + FILE_POG)


    # # step 1
    pog = mcl_to_pog(FILE_POG,FILE_PROTEIN_DIC)
    pickle.dump(pog,open(FILE_POG + '.pickle','wb'))

    # step 2. select representative proteins
    # 1st sample - 1% per group
    pog = pog[1]
    temp = DplyFrame(pog) >> group_by(X.pog) >> sample_frac(REP_PROTEIN_FRAC)

    # 2nd sample - 1 per group for those 1%<1
    if temp.shape[0] < pog.shape[0]:
        temp2 = DplyFrame(pog) >> \
        sift(X.pog.isin(np.setdiff1d(pog.pog.unique(),temp.pog.unique()))) >> \
        group_by(X.pog) >> sample_n(1)
        repseq = pd.concat((temp,temp2),axis=0)
    else:
        repseq = temp

    # write and select sequences
    repseq.protein.to_csv(FILE_REPSEQ,sep='\t',index=False)
    command = ['seqtk subseq',
               FAA,
               FILE_REPSEQ,
               '>',
               FILE_REPSEQ_FAA]
    command = ' '.join(command)
    os.system(command)

    # step 3. diamond alignment
    # meta
    command = ['time diamond blastp',
              '--db /mnt/data1/menghanliu/data/biocyc_annotation/tier1/meta/25.1/data/protseq.fsa',
              '--query', FILE_REPSEQ_FAA, '--out', FILE_REPSEQ_ALIGNMENT_META,
               '--max-target-seqs','1',
              # '--evalue 1e-10',
              '--very-sensitive']
    command = ' '.join(command)
    os.system(command)

    # step 4. annotate diamond alignment

    out = read_diamond_output(FILE_REPSEQ_ALIGNMENT_META)
    out.loc[:,'protein'] = out.qseqid

    # annot  = left_join(repseq.loc[:,['pog','protein']],out.loc[:,['protein','sseqid']])
    annot = out.loc[:,['protein','sseqid']
    print(head(annot))
    annot = left_join(,
                      pd.read_table('/mnt/data1/menghanliu/data/biocyc_annotation/tier1/meta/25.1/protein_product.txt',
                                  names = ['sseqid','keywords']))
    print(head(print(head(annot))))
    annot.to_csv(out_full, index=False,sep="\t")
    print('saved to :' + out_full)
    annot = pd.read_table(out_full)

    annot_sum = pog_annotate_summary(annot)
    annot_sum.to_csv(out_summary, index=False,sep="\t")
    print('saved to :' + out_summary)
