#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import argparse
import glob
import os


def Combine_sgRNA_barcode_from_the_Same_mouse(input_folder_address):
    # find all file for each sgRNA
    temp_pattern = '/**/Clonal_barcode/*_cluster.csv' # When recursive is set, ** followed by a path separator matches 0 or more subdirectories.
    barcode_address_list = glob.glob(input_folder_address+temp_pattern, recursive=True)
    # address for sgRNA
    temp_ref_address = input_folder_address+'/Intermediate_df.csv'
    temp_bc_df_list = [] # store all df for each sgRNA into one list 
    for temp_a in barcode_address_list:
        t1 = temp_a.replace('cluster', 'barcode') # address for cluster.csv
        t2 = temp_a # address for barcode.csv
        t3 = temp_a.replace('_cluster.csv', '.bartender') # # address for bartender input
        temp_df = merge_barcode_and_sgRNA_output(t1,t2,t3)
        temp_bc_df_list.append(temp_df)
    temp_total_bc_df = pd.concat(temp_bc_df_list).reset_index(drop = True)
    temp_sgRNA_ref_df = pd.read_csv(temp_ref_address)
    temp_final = temp_total_bc_df.merge(temp_sgRNA_ref_df, on=['Read_ID','Clonal_barcode'], how = 'left')
    # deduplex
    temp_name_list = ['gRNA_combination', 'Clonal_barcode_center',
                  'gRNA1','gRNA2','gRNA3', 'Clonal_barcode', 'Sample_ID']
    temp_final_raw = temp_final.groupby(temp_name_list, as_index=False)['Read_ID'].count()
    temp_final_raw.rename(columns={'Read_ID':'Count'}, inplace=True)
    temp_final_true_sgRNA = temp_final # I only keep those perfect match, which I already did for cas12a
    temp_final_completely_deduplexed = temp_final_true_sgRNA.groupby(['gRNA_combination', 'Clonal_barcode_center','gRNA1','gRNA2','gRNA3', 'Sample_ID'], as_index=False)['Read_ID'].count()
    temp_final_completely_deduplexed.rename(columns={'Read_ID':'Count',
                                                    'Clonal_barcode_center':'Clonal_barcode'}, inplace=True)
    return (temp_final_raw, temp_final_completely_deduplexed)

def merge_barcode_and_sgRNA_output(input_barcode_address,input_cluster_address,bartender_input):
    temp_df1 = pd.read_csv(input_barcode_address).drop(columns = ['Frequency'])
    temp_df2 = pd.read_csv(input_cluster_address).drop(columns = ['Cluster.Score','time_point_1'])
    temp_df3 = pd.read_csv(bartender_input, sep =',', names=['Clonal_barcode','Read_ID'], header=None)
    temp_merge=pd.merge(temp_df1, temp_df2, how='inner', on=['Cluster.ID'],
         left_index=False, right_index=False, sort=True, copy=True, indicator=False,
         validate=None)
    temp_merge.rename(columns={'Unique.reads':'Clonal_barcode',
                              'Center':'Clonal_barcode_center'}, inplace=True)
    temp_merge = temp_merge.drop(columns = ['Cluster.ID']).merge(temp_df3, on='Clonal_barcode', how='right')
    return(temp_merge)
    
#!/usr/bin/env python
# coding: utf-8

import argparse
import os

def main():
    """
    Combines sgRNA and clonal barcode information from the input sample folder 
    and saves the results to specified output files.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Combine sgRNA and clonal barcode information.')
    parser.add_argument("--a", required=True, help="Input sample folder containing data to process.")
    parser.add_argument("--o", required=True, help="Output folder to save combined results.")
    args = parser.parse_args()

    input_folder = args.a
    output_prefix = args.o

    # Validate input folder
    if not os.path.isdir(input_folder):
        raise FileNotFoundError(f"Error: The input folder '{input_folder}' does not exist.")

    # Ensure the output directory exists
    os.makedirs(output_prefix, exist_ok=True)
    
    # Extract sample ID from the input folder
    sample_id = os.path.basename(os.path.normpath(input_folder))
    print(f"Processing sample: {sample_id}")

    # Define output file paths
    combined_nd_file = os.path.join(output_prefix, 'Combined_ND_df.csv')
    combined_deduplexed_file = os.path.join(output_prefix, 'Combined_deduplexed_df.csv')

    # Call the data processing function
    try:
        df1, df2 = Combine_sgRNA_barcode_from_the_Same_mouse(input_folder)
    except Exception as e:
        raise RuntimeError(f"Error during data processing for sample {sample_id}: {e}")

    # Save results to CSV files
    try:
        df1.to_csv(combined_nd_file, index=False)
        print(f"Saved Combined ND DataFrame to: {combined_nd_file}")
        df2.to_csv(combined_deduplexed_file, index=False)
        print(f"Saved Combined Deduplexed DataFrame to: {combined_deduplexed_file}")
    except Exception as e:
        raise RuntimeError(f"Error while saving output files: {e}")

if __name__ == "__main__":
    main()



