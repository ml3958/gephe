#!/usr/bin/env python

# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: Feb 21, 2022
# Change the header of .faa files for all .faa files in a folder
# produce


import argparse
from Bio import SeqIO
import re, plotnine, sklearn,json, itertools, os, pickle, sys,glob
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

def change_genome_protein(input_dir,output_dir,prefix):

    # input_dir = './faa/'
    # output_dir = './faa_corrected/'

    if not os.path.exists(output_dir): os.mkdir(output_dir)
    print(' '.join(['Converting files in:',input_dir]))
    print(' '.join(['Saving converted files to:',output_dir]))
    print(' '.join(['prefix','is:',prefix]))

    # index genome names
    files = glob.glob(input_dir+'/*')
    genomes = [f.replace(input_dir,'') for f in files]
    genome_dic = { g : prefix+str(i) for i,g in enumerate(genomes)}
    pd.DataFrame([genome_dic]).transpose().reset_index().iloc[:,[1,0]].to_csv(output_dir + '/genomes.dic',sep='\t', index=False)

    # correct genome names and protein names
    # save to new directory
    for i,g in enumerate(genomes):
        g = output_dir + '/' + genome_dic[g] + '.faa'
        command = ['seqtk','rename',files[i],
                   ''.join(['"',genome_dic[genomes[i]]+'|p','"']),
                   ">",output_dir + '/' + genome_dic[genomes[i]] + '.faa']
        command = " ".join(command)
        os.system(command)

    # protein dictionary for each file
    # write to new directory
    for i,f in enumerate(files):
        original = f
        corrected = output_dir + '/' + genome_dic[genomes[i]] + '.faa'
        original = SeqIO.parse(original,'fasta')
        corrected = SeqIO.parse(corrected,'fasta')
        output = output_dir + '/' + genome_dic[genomes[i]] + '.proteindic'
        pd.DataFrame([[i.name for i in corrected],
                      [i.name for i in original]]).transpose().to_csv(output,sep='\t',index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('dir_input',help='folder containing faa files')
    parser.add_argument('dir_output',help='folder to save processed faa files')
    parser.add_argument('prefix',help='prefix')

    args = parser.parse_args()
    change_genome_protein(args.dir_input,args.dir_output,args.prefix)
