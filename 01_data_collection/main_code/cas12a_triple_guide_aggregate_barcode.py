#!/usr/bin/env python
# coding: utf-8
"""
README / Documentation
----------------------
This script combines sgRNA and clonal barcode information from a single sample folder.
It processes data files produced in earlier steps (such as sgRNA extraction, clustering, and bartender input)
to generate consolidated outputs that merge barcode and sgRNA details.

Usage:
------
Run the script from the command line with the following required arguments:
    --a : Path to the input sample folder containing data files.
    --o : Path to the output folder where the combined results will be saved.

Example:
--------
    python Combine_sgRNA_and_barcode.py --a path/to/sample_folder --o path/to/output_folder

Outputs:
--------
- Combined_ND_df.csv: A CSV file with the merged data from sgRNA and clonal barcode outputs, containing detailed
  information for each record. This file includes all entries from the merging process.
  
- Combined_deduplexed_df.csv: A deduplicated CSV file that groups records by unique sgRNA combinations,
  clonal barcode, and other identifiers, providing a count for each group.

Processing Overview:
--------------------
1. Search the input folder recursively for barcode-related CSV files.
2. For each identified file, determine the corresponding cluster and bartender input files based on naming conventions.
3. Merge the barcode, cluster, and bartender data using helper functions.
4. Concatenate all merged data into a consolidated DataFrame.
5. Merge the consolidated DataFrame with the sgRNA reference data from "Intermediate_df.csv".
6. Group and deduplicate the data to generate both detailed and summarized outputs.
7. Save the combined outputs as CSV files in the specified output folder.

File Dependencies:
------------------
- The input sample folder must contain subdirectories with the following files:
    * Barcode files: CSV files with barcode data (derived by replacing 'cluster' with 'barcode').
    * Cluster files: CSV files with a filename ending in "_cluster.csv".
    * Bartender input files: .bartender files providing mapping from Read_ID (identifier from the sequencer) to clonal barcode 
    * Intermediate_df.csv: A reference file generated from previous analysis that contains all the information for each read.


"""

import pandas as pd
import argparse
import glob
import os

def Combine_sgRNA_barcode_from_the_Same_mouse(input_folder_address):
    """
    Combines sgRNA and clonal barcode data from files within the specified sample folder.

    Steps:
    1. Recursively search for barcode-related CSV files matching the specified pattern.
    2. For each file, determine the corresponding cluster and bartender input files.
    3. Merge data from these files using the merge_barcode_and_sgRNA_output helper function.
    4. Concatenate all merged DataFrames into one.
    5. Merge the combined data with the sgRNA reference DataFrame from "Intermediate_df.csv".
    6. Group and deduplicate data to generate both detailed and summarized outputs.

    Parameters:
    -----------
    input_folder_address : str
        Path to the input sample folder containing the necessary data files.

    Returns:
    --------
    tuple of pandas.DataFrame:
        - The first DataFrame contains the raw combined data.
        - The second DataFrame contains the deduplicated and summarized data.
    """
    # Define the pattern to locate cluster CSV files in any subdirectory under Clonal_barcode.
    temp_pattern = '/**/Clonal_barcode/*_cluster.csv'  # '**' matches 0 or more subdirectories.
    barcode_address_list = glob.glob(input_folder_address + temp_pattern, recursive=True)

    # Path to the sgRNA reference file
    temp_ref_address = os.path.join(input_folder_address, 'Intermediate_df.csv')
    temp_bc_df_list = []  # List to store merged DataFrames for each sgRNA.

    for temp_a in barcode_address_list:
        # Determine file paths based on naming conventions:
        # t1: Path for the barcode data file (replace 'cluster' with 'barcode').
        t1 = temp_a.replace('cluster', 'barcode')
        # t2: Path for the cluster file (original file).
        t2 = temp_a
        # t3: Path for the bartender input file (replace '_cluster.csv' with '.bartender').
        t3 = temp_a.replace('_cluster.csv', '.bartender')
        
        # Merge the three files into one DataFrame.
        temp_df = merge_barcode_and_sgRNA_output(t1, t2, t3)
        temp_bc_df_list.append(temp_df)

    # Concatenate all merged DataFrames.
    temp_total_bc_df = pd.concat(temp_bc_df_list).reset_index(drop=True)
    
    # Read the sgRNA reference data.
    temp_sgRNA_ref_df = pd.read_csv(temp_ref_address)
    
    # Merge the combined data with the sgRNA reference DataFrame based on 'Read_ID' and 'Clonal_barcode'.
    temp_final = temp_total_bc_df.merge(temp_sgRNA_ref_df, on=['Read_ID', 'Clonal_barcode'], how='left')
    
    # Group by key columns to obtain a detailed merged output.
    temp_name_list = ['gRNA_combination', 'Clonal_barcode_center',
                      'gRNA1', 'gRNA2', 'gRNA3', 'Clonal_barcode', 'Sample_ID']
    temp_final_raw = temp_final.groupby(temp_name_list, as_index=False)['Read_ID'].count()
    temp_final_raw.rename(columns={'Read_ID': 'Count'}, inplace=True)
    
    # Group to generate the deduplicated summary (perfect matches as pre-filtered for Cas12a).
    temp_final_true_sgRNA = temp_final
    temp_final_completely_deduplexed = temp_final_true_sgRNA.groupby(
        ['gRNA_combination', 'Clonal_barcode_center', 'gRNA1', 'gRNA2', 'gRNA3', 'Sample_ID'],
        as_index=False)['Read_ID'].count()
    temp_final_completely_deduplexed.rename(columns={'Read_ID': 'Count',
                                                       'Clonal_barcode_center': 'Clonal_barcode'}, inplace=True)
    return (temp_final_raw, temp_final_completely_deduplexed)

def merge_barcode_and_sgRNA_output(input_barcode_address, input_cluster_address, bartender_input):
    """
    Merges barcode, cluster, and bartender input files into a consolidated DataFrame.

    Steps:
    1. Read the barcode data and drop the 'Frequency' column.
    2. Read the cluster data and remove unnecessary columns.
    3. Read the bartender input file and assign column names.
    4. Merge the barcode and cluster DataFrames on 'Cluster.ID'.
    5. Rename columns for clarity and merge with the bartender data.

    Parameters:
    -----------
    input_barcode_address : str
        Path to the CSV file containing barcode data.
    input_cluster_address : str
        Path to the CSV file containing cluster data.
    bartender_input : str
        Path to the bartender input file containing clonal barcode and read ID pairs.

    Returns:
    --------
    pandas.DataFrame:
        Merged DataFrame with combined information from all input sources.
    """
    # Read barcode data and drop the 'Frequency' column.
    temp_df1 = pd.read_csv(input_barcode_address).drop(columns=['Frequency'])
    
    # Read cluster data and drop unneeded columns.
    temp_df2 = pd.read_csv(input_cluster_address).drop(columns=['Cluster.Score', 'time_point_1'])
    
    # Read bartender input data; assign column names.
    temp_df3 = pd.read_csv(bartender_input, sep=',', names=['Clonal_barcode', 'Read_ID'], header=None)
    
    # Merge barcode and cluster data on 'Cluster.ID'.
    temp_merge = pd.merge(temp_df1, temp_df2, how='inner', on=['Cluster.ID'], sort=True)
    
    # Rename columns for clarity.
    temp_merge.rename(columns={'Unique.reads': 'Clonal_barcode',
                                 'Center': 'Clonal_barcode_center'}, inplace=True)
    
    # Drop 'Cluster.ID' and merge with bartender input data.
    temp_merge = temp_merge.drop(columns=['Cluster.ID']).merge(temp_df3, on='Clonal_barcode', how='right')
    
    return temp_merge

def main():
    """
    Combines sgRNA and clonal barcode information from the input sample folder 
    and saves the results to specified output CSV files.

    Steps:
    1. Parse command line arguments for input and output directories.
    2. Validate the input folder and create the output folder if necessary.
    3. Extract the sample ID from the input folder.
    4. Process data by calling Combine_sgRNA_barcode_from_the_Same_mouse.
    5. Save the detailed and deduplicated outputs as CSV files.
    """
    # Set up argument parser.
    parser = argparse.ArgumentParser(description='Combine sgRNA and clonal barcode information.')
    parser.add_argument("--a", required=True, help="Input sample folder containing data to process.")
    parser.add_argument("--o", required=True, help="Output folder to save combined results.")
    args = parser.parse_args()

    input_folder = args.a
    output_prefix = args.o

    # Validate input folder.
    if not os.path.isdir(input_folder):
        raise FileNotFoundError(f"Error: The input folder '{input_folder}' does not exist.")

    # Create the output directory if it doesn't exist.
    os.makedirs(output_prefix, exist_ok=True)
    
    # Extract sample ID from the input folder.
    sample_id = os.path.basename(os.path.normpath(input_folder))
    print(f"Processing sample: {sample_id}")

    # Define output file paths.
    combined_nd_file = os.path.join(output_prefix, 'Combined_ND_df.csv')
    combined_deduplexed_file = os.path.join(output_prefix, 'Combined_deduplexed_df.csv')

    # Process data and combine sgRNA and barcode information.
    try:
        df1, df2 = Combine_sgRNA_barcode_from_the_Same_mouse(input_folder)
    except Exception as e:
        raise RuntimeError(f"Error during data processing for sample {sample_id}: {e}")

    # Save results to CSV files.
    try:
        df1.to_csv(combined_nd_file, index=False)
        print(f"Saved Combined ND DataFrame to: {combined_nd_file}")
        df2.to_csv(combined_deduplexed_file, index=False)
        print(f"Saved Combined Deduplexed DataFrame to: {combined_deduplexed_file}")
    except Exception as e:
        raise RuntimeError(f"Error while saving output files: {e}")

if __name__ == "__main__":
    main()
