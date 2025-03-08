#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import argparse
import glob
import os

    
def main():
    parser = argparse.ArgumentParser(description='A function to combine sgRNA and clonal barcode information')
    parser.add_argument("--o", required=True, help="This is the prefix of output file")
    args = parser.parse_args()
    temp_output_prefix = args.o
    MYDIR = temp_output_prefix
    folder_paths = []
    for entry_name in os.listdir(MYDIR):
        entry_path = os.path.join(MYDIR, entry_name)
        if os.path.isdir(entry_path):
            folder_paths.append(entry_path)
    temp_final_df_list = []
    for temp_folder in folder_paths:
        temp_sampleID = temp_folder.split('/')[-1]
        temp_o2 = temp_output_prefix + temp_sampleID + '/Combined_deduplexed_df_full.csv'
        temp_df2 = pd.read_csv(temp_o2)
        temp_final_df_list.append(temp_df2)
    temp_final_df = pd.concat(temp_final_df_list)
    temp_final_df.to_csv(temp_output_prefix + 'gRNA_clonalbarcode_combined_full.csv', index=False)
if __name__ == "__main__":
    main()  


