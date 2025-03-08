#!/usr/bin/env python
# coding: utf-8
"""
README / Documentation
----------------------
This script consolidates sgRNA and clonal barcode information from multiple sample subfolders into a single CSV file.
This script generates a merged DataFrame for all data, including combinations that were not designed.
It is intended to be used for quality control (QC) to assess the probability of template switching between guides.

Usage:
------
Run the script from the command line with the following required arguments:
    --a : Prefix (directory) of input files containing sample subfolders.
    --o : Prefix (directory) of output files where the combined result will be saved.

Example:
--------
    python combine_full.py --a path/to/input_prefix/ --o path/to/output_prefix/

Outputs:
--------
- gRNA_clonalbarcode_combined_full.csv: A CSV file that merges the data from all sample subfolders.

Processing Overview:
--------------------
1. Identify all subdirectories in the input prefix directory.
2. For each subdirectory, read the file "Combined_deduplexed_df_full.csv".
3. Concatenate all these CSV files into a single DataFrame.
4. Save the final combined DataFrame as "gRNA_clonalbarcode_combined_full.csv" in the output prefix.

"""

import pandas as pd
import argparse
import os

def main():
    """
    Consolidates sgRNA and clonal barcode data from multiple sample subfolders into a single CSV file.

    Processing Steps:
    -----------------
    1. Parse command line arguments to obtain the input and output prefixes.
    2. Identify all subdirectories within the input prefix directory.
    3. For each subdirectory:
         - Extract the sample ID (the folder name).
         - Construct the path to "Combined_deduplexed_df_full.csv" within the corresponding output folder.
         - Read the CSV file and append the DataFrame to a list.
    4. Concatenate all DataFrames into one final DataFrame.
    5. Save the combined DataFrame as "gRNA_clonalbarcode_combined_full.csv" in the output prefix.

    Important:
    ----------
    This process generates a merged DataFrame that includes all combinations of sgRNA and clonal barcodes,
    even those combinations that were not originally designed. This comprehensive dataset can be used for QC,
    specifically to evaluate the probability of template switching between guides.
    """
    parser = argparse.ArgumentParser(description='Combine sgRNA and clonal barcode information from multiple samples.')
    parser.add_argument("--a", required=True, help="Prefix of input files (directory containing sample subfolders).")
    parser.add_argument("--o", required=True, help="Prefix of output files (directory where combined output will be saved).")
    args = parser.parse_args()
    
    temp_input_prefix = args.a
    temp_output_prefix = args.o

    # List all subdirectories in the input prefix.
    folder_paths = []
    for entry_name in os.listdir(temp_input_prefix):
        entry_path = os.path.join(temp_input_prefix, entry_name)
        if os.path.isdir(entry_path):
            folder_paths.append(entry_path)
    
    # Initialize list to store DataFrames from each sample.
    temp_final_df_list = []
    for temp_folder in folder_paths:
        # Extract the sample ID from the folder name.
        temp_sampleID = os.path.basename(temp_folder)
        # Construct the path to the sample's combined deduplicated CSV file.
        temp_o2 = os.path.join(temp_input_prefix, temp_sampleID, 'Combined_deduplexed_df_full.csv')
        # Read the CSV file.
        temp_df2 = pd.read_csv(temp_o2)
        temp_final_df_list.append(temp_df2)
    
    # Concatenate all DataFrames into one.
    temp_final_df = pd.concat(temp_final_df_list, ignore_index=True)
    # Construct the output file path.
    output_file_path = os.path.join(temp_output_prefix, 'gRNA_clonalbarcode_combined_full.csv')
    # Save the final combined DataFrame.
    temp_final_df.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    main()
