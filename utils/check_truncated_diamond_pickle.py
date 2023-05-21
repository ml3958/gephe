#!/usr/bin/env python
import argparse
import os, pickle, sys, glob
import pandas as pd
import subprocess

# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: May 21, 2023

# def read_diamond_output(diamond_file):
#     diamond = pd.read_table(diamond_file,
#                             header=None,sep='\t'
#                             # names=['qseqid','sseqid','pident','length','mismatch','gapopen','qlen','qstart','qend','slen','start','send','evalue','bitscore'] # column already exists for customized diamond output format
#                             )
#     return(diamond)



if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('file_pickle', help='.diamond.out')
    args = parser.parse_args()

    file  = args.file_pickle

    try:
        # Command that may raise an error
        result = pd.read_pickle(open('test.pickle','rb'))
        print(file + ' ok')
    except Exception as e:
        # Handle the error
        # print("An error occurred:", str(e))
        print(file + ' truncated, removing....')
        command = "rm -rf " + file
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
