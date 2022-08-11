#!/usr/bin/env python
import argparse
import os, pickle, sys, glob
import pandas as pd
import glob
import os
import numpy as np
import itertools
import pandas as pd
from tqdm import tqdm

# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: July 8, 2022
# USAGE: python $gephe_dir/alignment/diamond_to_pickle.py $METADATA_POS $DIR_INPUT $DIR_FAA_MERGE



def divide_diamondout(FILE, DIR_OUT, split_char='_', chunksize = 10 ** 6, test_only = False,
                     diamond_colnames = ['qseqid','sseqid','pident','length','mismatch','gapopen','qlen','qstart','qend','slen','start','send','evalue','bitscore']):

    # this script divide a LARGE diamond output (usually generated from 200 genomes) into

    i=0
    for chunk in pd.read_table(FILE, chunksize=chunksize, header=None, names=diamond_colnames):
        dat = chunk
        index = dat.qseqid.map(lambda x: x.split(split_char)[0]) # index is the genome names
        for gem in index.unique():
            f = DIR_OUT + gem + '.diamond.out'
            df = dat[index==gem]
            if not os.path.isfile(f):
                df.to_csv(f, sep='\t', index=False)
            else: # else it exists so append without writing the header
                df.to_csv(f, sep='\t', mode='a',header=False, index=False)  # when appending, header is discarded
        # for i in dat.qseqid.map(lambda x: x.split('|')[0]).unique()
        i+=1

        if test_only and i==2:
            break


if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Merge faa files')
    parser.add_argument('diamond_file_merged', help='file name of merged diamond output')
    parser.add_argument('dir_output', help='Path to save merged .faa files')
    parser.add_argument('-sep', default='|', help='Seperator [genome|protein]')
    args = parser.parse_args()

    diamond_file_merged = args.diamond_file_merged
    dir_output = args.dir_output
    sep = args.sep

    divide_diamondout(diamond_file_merged,dir_output, split_char=sep)
