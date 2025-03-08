#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import argparse
import os

def main():
    """
    Combines sgRNA and clonal barcode information into a single CSV file.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Combine sgRNA and clonal barcode information.')
    parser.add_argument("--o", required=True, help="Output folder prefix for combined file.")
    args = parser.parse_args()

    # Define output directory and initialize variables
    output_prefix = args.o
    if not os.path.isdir(output_prefix):
        raise FileNotFoundError(f"Error: The specified output directory '{output_prefix}' does not exist.")

    folder_paths = [
        os.path.join(output_prefix, entry_name)
        for entry_name in os.listdir(output_prefix)
        if os.path.isdir(os.path.join(output_prefix, entry_name))
    ]

    combined_dataframes = []

    # Process each folder and read its respective CSV file
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

    # Combine all dataframes into a single one
    if combined_dataframes:
        final_combined_df = pd.concat(combined_dataframes, ignore_index=True)
        output_file_path = os.path.join(output_prefix, 'gRNA_clonalbarcode_combined.csv')
        final_combined_df.to_csv(output_file_path, index=False)
        print(f"Combined data saved to: {output_file_path}")
    else:
        print("No valid dataframes found to combine. Exiting...")

if __name__ == "__main__":
    main()
