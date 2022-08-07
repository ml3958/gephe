#!/usr/bin/env python
import argparse
import os, pickle, sys, glob
import pandas as pd

# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: July 8, 2022
# USAGE: python $gephe_dir/alignment/diamond_to_pickle.py $METADATA_POS $DIR_INPUT $DIR_FAA_MERGE

if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Merge faa files')
    parser.add_argument('metadata_pos', help='metadata_pos')
    parser.add_argument('dir_input', help='Path to .faa files')
    parser.add_argument('dir_output', help='Path to save merged .faa files')
    parser.add_argument('-n_faa_to_merge', default=200, help='${N_FAA_TO_MERGE}')
    args = parser.parse_args()

    metadata_pos = args.metadata_pos
    dir_input = args.dir_input
    dir_output = args.dir_output
    n_faa_to_merge = args.n_faa_to_merge

    files = glob.glob(dir + '/*diamond.out')

    for f in files:
        if os.path.exists(f+'.pickle'):
            print(f+'.pickle exists, skipping')
            next
        diamond = read_diamond_output(f)
        print('saving ' + f+ '.pickle...')
        pickle.dump(diamond,open(f+'.pickle','wb'))
