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


import scipy.cluster.hierarchy as shc
import dfply


def read_diamond_output(diamond_file):
    diamond = pd.read_table(diamond_file,
                            header=None,sep='\t',
                            names=['qseqid','sseqid','pident','length','mismatch','gapopen','qlen','qstart','qend','slen','start','send','evalue','bitscore'])
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
    protein_annot.protein_name.fillna('',inplace=True)
    # protein_annot =
    pog_N = DplyFrame(protein_annot)>> group_by(X.pog) >> summarize(N = X.pog.count())
    pog_protein_name = DplyFrame(protein_annot) >>  group_by(X.pog,X.protein_name) >> summarize(keyword_N = X.pog.count()) >> \
        sift(X.keyword_N>annot_N) >> mutate(protein_name = X.protein_name + ' (n='+X.keyword_N.astype(str) +')') >> \
        ungroup() >> arrange(-X.keyword_N)>> group_by(X.pog) >>summarize(protein_name = X.protein_name.unique()) >> \
        mutate(protein_name = X.protein_name.map(lambda x: dplyr_paste(x,'; ')))

    result = full_join(pog_N,pog_protein_name).sort_values(column_for_sorting,ascending=False)
    return(result)


def dplyr_paste(x,sep=','):
    if type(x)==str: return x
    else: return sep.join(x)




if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='pog functional annotation')
    parser.add_argument('dir_alignment_master',help='alignment folder')
    parser.add_argument('dir_pog',help='pog folder')
    parser.add_argument('prefix_protein',help='${PREFIX_PROTEIN}')
    parser.add_argument('prefix_pog',help='${PREFIX_POG}')
    parser.add_argument('-f',dest='pog_protein_fraction',type=float,help='${PREFIX_POG}')


    args = parser.parse_args()

    DIR_ALIGNMENT_MASTER = args.dir_alignment_master
    DIR_POG = args.dir_pog
    PREFIX_PROTEIN = args.prefix_protein
    PREFIX_POG = args.prefix_pog
    REP_PROTEIN_FRAC = args.pog_protein_fraction

    # input
    FILE_POG = DIR_POG + PREFIX_POG + '.mcl'
    FILE_PROTEIN_DIC = DIR_POG + '/' + PREFIX_PROTEIN + '.dic'
    FAA = DIR_ALIGNMENT_MASTER + '/input.fasta'

    # output
    FILE_PROTSEQ = DIR_POG + PREFIX_PROTEIN  + '.protseq'
    FILE_PROTSEQ_FAA = DIR_POG + PREFIX_PROTEIN + '.protseq.faa'

    FILE_PROTSEQ_ALIGNMENT_META = DIR_POG + PREFIX_PROTEIN + '.protseq.faa_biocycmeta.diamond.out'
    FILE_PROTSEQ_ALIGNMENT_UNIREF = DIR_POG + PREFIX_PROTEIN + '.protseq.faa_uniref50.diamond.out'
    # FILE_PROTSEQ_ALIGNMENT_META_ANNOTATION = DIR_POG + PREFIX_POG + '.annot_meta'
    # FILE_PROTSEQ_ALIGNMENT_ECOLI = DIR_POG + PREFIX_POG + '.protseq.faa_ecoli.diamond.out'
    # FILE_PROTSEQ_ALIGNMENT_ECOLI_ANNOTATION = DIR_POG + PREFIX_POG + '.annot_ecoli'

    out_full = DIR_POG + '/' + PREFIX_PROTEIN + '.annot'
    out_summary = DIR_POG + '/' + PREFIX_POG + '.annot'

    print('Processing :' + FILE_POG)

    # # step 1
    pog = mcl_to_pog(FILE_POG,FILE_PROTEIN_DIC)
    if not os.path.exists(FILE_POG + '.pickle'):
        pickle.dump(pog,open(FILE_POG + '.pickle','wb'))

    # step 2. map selected proteins to biocyc tier 1 database (skip if done previously)
    if os.path.exists(out_full):
        print(out_full+' exists --> skipping the subseq, alignment, protein annotation')
    else:
        pog[1].protein.to_csv(FILE_PROTSEQ,sep='\t',index=False)
        command = ['seqtk subseq',
           FAA,
           FILE_PROTSEQ,
           '>',
           FILE_PROTSEQ_FAA]
        command = ' '.join(command)
        os.system(command)

        # diamond alignment to
        if not os.path.exists(FILE_PROTSEQ_ALIGNMENT_META):
            command = ['time diamond blastp','--db /mnt/data1/menghanliu/data/biocyc_annotation/tier1/combined/combined_protseq_diamond.dmnd',
            '--query', FILE_PROTSEQ_FAA, '--out', FILE_PROTSEQ_ALIGNMENT_META,
               '--max-target-seqs','1',
               '--query-cover', '66',
              '--evalue 1e-10',
               '--outfmt','6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore',
              '--very-sensitive']
            command = ' '.join(command)
            os.system(command)

        if not os.path.exists(FILE_PROTSEQ_ALIGNMENT_UNIREF):
            command = ['time diamond blastp',
              '--db /mnt/data1/menghanliu/data/protein_database/uniref/uniref50.fasta',
              '--query', FILE_PROTSEQ_FAA, '--out', FILE_PROTSEQ_ALIGNMENT_UNIREF,
               '--max-target-seqs','1',
               '--query-cover', '66',
              '--evalue 1e-10',
               '--outfmt','6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore',
              '--very-sensitive']
            command = ' '.join(command)
            os.system(command)

        # annotate hits
        out = read_diamond_output(FILE_PROTSEQ_ALIGNMENT_META)

        annot = pd.DataFrame()
        annot['protein'] = out.qseqid
        temp = out.sseqid.str.split('|',expand=True)
        annot[['species','proteinID']] = temp[[1,2]]


        annot = left_join(DplyFrame(annot),
                  pd.read_table('/mnt/data1/menghanliu/data/biocyc_annotation/tier1/combined/combined_protein-links.dat'))
        annot = annot.drop_duplicates()
        print('saved to :' + out_full)
        annot.to_csv(out_full, index=False,sep="\t")


    # step 3. pog annotation
    if os.path.exists(out_summary):
        print(out_summary+' exists --> skipping the subseq, alignment, protein annotation')
    else:
        annot = DplyFrame(pd.read_table(out_full))
        annot_sum = pog_annotate_summary(left_join(annot,pog[1]),  annot_N=0)
        # annot_sum
        annot_sum.to_csv(out_summary, index=False,sep="\t")
        print('saved to :' + out_summary)
