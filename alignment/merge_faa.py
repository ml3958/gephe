#!/usr/bin/env python
import argparse
import os, pickle, sys, glob
import pandas as pd

# Author: Menghan Liu @ Tavazoie lab at Columbia University
# Date: July 8, 2022
# USAGE: python $gephe_dir/alignment/diamond_to_pickle.py $METADATA_POS $DIR_INPUT $DIR_FAA_MERGE

def merge_faa(dir_faa, dir_faa_merge, n_file_to_merge = 200):

    os.chdir(dir_faa)

    # all faa files
    # --------------
    faa_files = glob.glob(dir_faa + '/*faa')
    # print(faa_files)

    # group faa files into chuncks, with each chunck roughly 200 files
    # --------------
    n_chucks = round(len(faa_files)/n_file_to_merge)
    print(n_chucks)
    faa_files_split = np.array_split(faa_files, n_chucks)
    # print([len(i) for i in faa_files_split])

    # make new folder to store merged files
    # --------------
    if not os.path.exists(dir_faa_merge):
        os.makedirs(dir_faa_merge)

    # save the map from faa to merged faa
    # --------------
    faa_dic = pd.DataFrame(dict(V1=faa_files,
                                V2=list(itertools.chain.from_iterable([[dir_faa_merge + '/merged_'+ str(i) + '.faa'] * len(files) for i,files in enumerate(faa_files_split)]))))
    faa_dic.to_csv(dir_faa_merge + '/faa.map')

    # merge faa
    # --------------
    merge_command = ['cat ' + ' '.join(list(files)) + '>' + dir_faa_merge + '/merged_'+ str(i) + '.faa' for i,files in enumerate(faa_files_split)]
    # print(merge_command[0])

    for i in tqdm(merge_command):
        os.system(i)


if __name__ == '__main__':
    # parse input argument
    parser = argparse.ArgumentParser(description='Merge faa files')
    parser.add_argument('dir_input', help='Path to .faa files')
    parser.add_argument('dir_output', help='Path to save merged .faa files')
    parser.add_argument('-n_faa_to_merge', default=200, help='${N_FAA_TO_MERGE}')
    args = parser.parse_args()

    dir_input = args.dir_input
    dir_output = args.dir_output
    n_faa_to_merge = args.n_faa_to_merge

    merge_faa(dir_input,dir_output,n_faa_to_merge)
