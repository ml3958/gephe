#!/usr/bin/env python
import argparse
import os, pickle, sys, glob
import pandas as pd

# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: June 6, 2021

def read_diamond_output(diamond_file):
    diamond = pd.read_table(diamond_file,
                            header=None,sep='\t',
                            names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','start','send','evalue','bitscore'])
    return(diamond)


if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Change fasta name for jgi [*.genes.faa] files ')
    parser.add_argument('dir', help='In/output directory')
    args = parser.parse_args()

    dir = args.dir

    files = glob.glob(dir + '/*diamond.out')

    for f in files:
        print(f)
        diamond = read_diamond_output(f)
        pickle.dump(diamond,open(f+'.pickle','wb'))
