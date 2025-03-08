#!/usr/bin/env python
# coding: utf-8
"""
README / Documentation
----------------------
This script consolidates sgRNA and clonal barcode information from multiple sample subfolders into a single CSV file.
It searches through a specified output folder for subdirectories containing deduplicated combined files,
merges these files, and produces a single consolidated output file.

Usage:
------
Run the script from the command line with the following required argument:
    --o : Path to the output folder containing subfolders for each sample.

Example:
--------
    python Combine_All_Samples.py --o path/to/output_folder

Outputs:
--------
- gRNA_clonalbarcode_combined.csv: A CSV file combining the deduplicated sgRNA and clonal barcode information 
  from all sample subfolders.

Processing Overview:
--------------------
1. Verify that the provided output folder exists.
2. Identify all subdirectories within the output folder.
3. For each subdirectory, locate the file "Combined_deduplexed_df.csv" that contains the processed data for that sample.
4. Read and collect each CSV file into a list of DataFrames.
5. Concatenate all DataFrames into a single DataFrame.
6. Save the final combined DataFrame as "gRNA_clonalbarcode_combined.csv" in the specified output folder.

File Dependencies:
------------------
- Each sample subfolder within the output folder must contain a file named "Combined_deduplexed_df.csv", 
  which is generated from prior processing steps.
"""

import pandas as pd
import argparse
import os

def main():
    """
    Consolidates sgRNA and clonal barcode information into a single CSV file.
    
    Processing Steps:
    -----------------
    1. Parse the command line argument to obtain the output folder path.
    2. Validate that the output folder exists.
    3. Identify all subdirectories (each representing a sample) within the output folder.
    4. For each subdirectory:
         - Construct the path to the 'Combined_deduplexed_df.csv' file.
         - Read the CSV file if it exists; otherwise, log a warning and skip the sample.
    5. Concatenate all valid DataFrames into a final combined DataFrame.
    6. Save the consolidated DataFrame as 'gRNA_clonalbarcode_combined.csv' in the output folder.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Combine sgRNA and clonal barcode information from multiple samples.')
    parser.add_argument("--o", required=True, help="Output folder prefix containing subfolders with sample data.")
    args = parser.parse_args()

    # Define the output directory and verify its existence
    output_prefix = args.o
    if not os.path.isdir(output_prefix):
        raise FileNotFoundError(f"Error: The specified output directory '{output_prefix}' does not exist.")

    # Identify all subdirectories (samples) within the output folder
    folder_paths = [
        os.path.join(output_prefix, entry_name)
        for entry_name in os.listdir(output_prefix)
        if os.path.isdir(os.path.join(output_prefix, entry_name))
    ]

    combined_dataframes = []

    # Process each sample subfolder and read its respective CSV file
    for folder_path in folder_paths:
        sample_id = os.path.basename(folder_path)
        combined_csv_path = os.path.join(folder_path, 'Combined_deduplexed_df.csv')

        if not os.path.isfile(combined_csv_path):
            print(f"Warning: File not found for sample {sample_id}: {combined_csv_path}. Skipping...")
            continue

        try:
            df = pd.read_csv(combined_csv_path)
            combined_dataframes.append(df)
        except Exception as e:
            print(f"Error reading file {combined_csv_path} for sample {sample_id}: {e}. Skipping...")
            continue

    # Combine all DataFrames into a single one, if any were successfully loaded
    if combined_dataframes:
        final_combined_df = pd.concat(combined_dataframes, ignore_index=True)
        output_file_path = os.path.join(output_prefix, 'gRNA_clonalbarcode_combined.csv')
        final_combined_df.to_csv(output_file_path, index=False)
        print(f"Combined data saved to: {output_file_path}")
    else:
        print("No valid dataframes found to combine. Exiting...")

if __name__ == "__main__":
    main()
