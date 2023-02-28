#!/usr/bin/env python
import argparse
import os, pickle, sys, glob
import pandas as pd

# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: June 6, 2021

# def read_diamond_output(diamond_file):
#     diamond = pd.read_table(diamond_file,
#                             header=None,sep='\t'
#                             # names=['qseqid','sseqid','pident','length','mismatch','gapopen','qlen','qstart','qend','slen','start','send','evalue','bitscore'] # column already exists for customized diamond output format
#                             )
#     return(diamond)


if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('file_diamond', help='.diamond.out')
    args = parser.parse_args()

    file_diamond  = args.file_diamond

    if os.path.exists(file_diamond + '.pickle'):
        print(file_diamond +'.pickle exists, skipping')
        # diamond = read_diamond_output(f)
    else:
        diamond = pd.read_table(file_diamond)
        print('saving ' + file_diamond+ '.pickle...')
        pickle.dump(diamond,open(file_diamond+'.pickle','wb'))
