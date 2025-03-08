#!/usr/bin/env python
# coding: utf-8
"""
README / Documentation
----------------------
This script extracts guide RNA (gRNA) array and clonal barcode information from a paired-end (PE) sequencing FASTQ.gz file.
It uses a gRNA array reference CSV file to validate and annotate the extracted gRNA sequences.

Usage:
------
Run the script from the command line with the following required arguments:
    --a : Path to the input FASTQ.gz file containing the sequencing data.
    --b : Path to the gRNA array reference CSV file.
    --o : Path to the output directory where the results will be saved.

Example:
--------
    python UltraSeq_Step1_Cas12a_TripleArray.py --a path/to/input.fastq.gz \
                                                 --b path/to/sgRNA_Cas12a_combine.csv \
                                                 --o path/to/output_directory

Outputs:
--------
- Intermediate_df.csv: A CSV file containing detailed extracted information for each sequencing read, including gRNA sequences, clonal barcodes, read IDs, sample ID, and classification.
  * Each row corresponds to one sequencing read.
  * Only reads with gRNA sequences matching the reference at the corresponding positions are included.
  
- Combined_deduplexed_df_full.csv: A summarized CSV file that groups records by unique gRNA combinations and provides their counts.
  * This file represents a deduplicated version of Intermediate_df.csv.
  * This file is useful for QC templated switching event
  
- Bartender input files: The script generates .bartender files (one per unique gRNA combination) that contain clonal barcode and read ID pairs for downstream processing.
  * A text file listing the paths to these .bartender files is also created.
  * Only reads with gRNA combinations that match the reference (both individual gRNA sequences and their combination) are included.

gRNA Array Reference File Requirements:
-----------------------------------------
The gRNA array reference CSV file (e.g., sgRNA_Cas12a_combine.csv) is essential for accurate annotation of the gRNA sequences. This file must include the following columns:
    - DR1, DR2, DR3, DR6: Direct repeat (DR) patterns used to build the regular expression for sequence extraction.
      * Each DR column must contain only one unique value; the script will exit with an error if multiple patterns are found.
    - Guide1_sequence, Guide2_sequence, Guide3_sequence: The expected guide RNA sequences corresponding to each of the three guides.
    - gRNA_combination: A concatenated string representing the combined gRNA sequence.
      
Note:
-----
- The regular expression pattern used for extraction is constructed as follows:
      TAGTT (16bp clonal barcode) TATGG DR1 (23bp) DR2 (23bp) DR3 (23bp) DR6
  where DR1, DR2, DR3, and DR6 are dynamically inserted from the reference file.
- Ensure that the reference file is correctly formatted to avoid errors during sequence extraction or annotation.

Processing Overview:
--------------------
1. Load and validate the gRNA array reference CSV file.
2. Extract the unique DR patterns (DR1, DR2, DR3, DR6) and construct the regular expression.
3. Open and read the input FASTQ.gz file to process each sequencing read.
4. Extract a 16-bp clonal barcode and three 23-bp gRNA sequences using the compiled pattern.
5. Verify the extracted gRNA sequences against the expected sequences in the reference file.
6. Generate detailed and summarized output CSV files.
7. Create .bartender files for downstream processing by grouping reads based on unique gRNA combinations.
"""


import gzip
import regex
import argparse
import sys
import pandas as pd

def find_reverse_complementary(input_string):
    temp_dic = {'A':'T','G':'C','T':'A','C':'G','N':'N',
                'a':'t','g':'c','t':'a','c':'g','n':'n'}
    return(''.join(tuple([temp_dic.get(x) for x in input_string[::-1]])))

def main():
    parser = argparse.ArgumentParser(description='Extract gRNA and clonal barcode from merged FASTQ.gz file')
    parser.add_argument("--a", required=True, help="Path to the input FASTQ.gz file")
    parser.add_argument("--b", required=True, help="Path to the gRNA array reference CSV file (e.g., sgRNA_Cas12a_combine.csv)")
    parser.add_argument("--o", required=True, help="Directory for output files")
    args = parser.parse_args()
    
    fastqgz_input_address = args.a
    ref_address = args.b
    output_dir = args.o

    # Load the gRNA array reference file
    ref_sgRNA_df = pd.read_csv(ref_address)  # reference DataFrame

    # Validate that each DR pattern column has a single unique value
    if len(ref_sgRNA_df.DR1.unique()) > 1:
        sys.exit("Error! More than one DR1 pattern found in the reference file.")
    else:
        DR1 = ref_sgRNA_df.DR1.unique()[0]

    if len(ref_sgRNA_df.DR2.unique()) > 1:
        sys.exit("Error! More than one DR2 pattern found in the reference file.")
    else:
        DR2 = ref_sgRNA_df.DR2.unique()[0]

    if len(ref_sgRNA_df.DR3.unique()) > 1:
        sys.exit("Error! More than one DR3 pattern found in the reference file.")
    else:
        DR3 = ref_sgRNA_df.DR3.unique()[0]

    if len(ref_sgRNA_df.DR6.unique()) > 1:
        sys.exit("Error! More than one DR6 pattern found in the reference file.")
    else:
        DR6 = ref_sgRNA_df.DR6.unique()[0]

    # Construct the regular expression pattern for extraction
    # The pattern is:
    # TAGTT (16bp barcode) TATGG DR1 (23bp) DR2 (23bp) DR3 (23bp) DR6
    temp_pattern1 = regex.compile('TAGTT' + '(.{16})' + 'TATGG' + DR1 + '(.{23})' + DR2 + '(.{23})' + DR3 + '(.{23})' + DR6)

    # Initialize counters and lists for collected data
    temp_total_read, temp_extracted_read, temp_matched_read = 0, 0, 0
    temp_sample_ID = output_dir.split('/')[-1]

    temp_gRNA1_list, temp_gRNA2_list, temp_gRNA3_list = [], [], []
    temp_read_ID_list, temp_Clonal_barcode_list, temp_label_list = [], [], []

    with gzip.open(fastqgz_input_address, 'rt') as handler1:
        # Read the first FASTQ record
        temp_readID = handler1.readline().rstrip()  # read ID
        temp_sequence1 = handler1.readline().rstrip()
        handler1.readline()  # skip +
        handler1.readline()  # skip quality scores

        while temp_readID:
            temp_total_read += 1
            # Apply the regular expression pattern to the sequence
            temp_search_result1 = temp_pattern1.search(temp_sequence1)
            if temp_search_result1:
                temp_clonal_barcode_r1 = temp_search_result1.group(1)
                temp_gRNA1_r1 = temp_search_result1.group(2)
                temp_gRNA2_r1 = temp_search_result1.group(3)
                temp_gRNA3_r1 = temp_search_result1.group(4)
                temp_extracted_read += 1
                temp_gRNA1_list.append(temp_gRNA1_r1)
                temp_gRNA2_list.append(temp_gRNA2_r1)
                temp_gRNA3_list.append(temp_gRNA3_r1)
                temp_read_ID_list.append(temp_readID)
                temp_Clonal_barcode_list.append(temp_clonal_barcode_r1)
                
                # Check if the extracted gRNAs match the expected sequences in the reference
                if ((temp_gRNA1_r1 in ref_sgRNA_df.Guide1_sequence.values) and
                    (temp_gRNA2_r1 in ref_sgRNA_df.Guide2_sequence.values) and
                    (temp_gRNA3_r1 in ref_sgRNA_df.Guide3_sequence.values)):
                    temp_matched_read += 1
                    temp_label_list.append('Expected')
                else:
                    temp_label_list.append('Unexpected')
            
            # Read the next FASTQ record (4 lines per record)
            temp_readID = handler1.readline().rstrip()  # read ID
            temp_sequence1 = handler1.readline().rstrip()
            handler1.readline()  # skip +
            handler1.readline()  # skip quality scores

    # Print a summary of the extraction process
    print(f"Sample {temp_sample_ID:s} has totally {temp_total_read:d} reads. Among them {temp_extracted_read:d} have barcode and sgRNA, which is about {temp_extracted_read/temp_total_read:.3f}. {temp_matched_read:d} have expected sgRNA, which is about {temp_matched_read/temp_total_read:.3f}")

    # Create a DataFrame with the extracted data
    Final_df = pd.DataFrame({
        'gRNA1': temp_gRNA1_list, 
        'gRNA2': temp_gRNA2_list,
        'gRNA3': temp_gRNA3_list,
        'Clonal_barcode': temp_Clonal_barcode_list, 
        'Read_ID': temp_read_ID_list,
        'Sample_ID': temp_sample_ID,
        'Class': temp_label_list
    })
    Final_df['gRNA_combination'] = Final_df['gRNA1'] + '_' + Final_df['gRNA2'] + '_' + Final_df['gRNA3']

    # Save detailed and summarized outputs
    Final_df.to_csv(f"{output_dir}/Intermediate_df.csv", index=False)
    Final_df_S = Final_df.groupby([col for col in Final_df.columns if col != 'Read_ID']).size().reset_index(name='Count')
    Final_df_S.to_csv(f"{output_dir}/Combined_deduplexed_df_full.csv", index=False)

    # Filter and generate .bartender input files for downstream processing
    Final_filtered_df = Final_df[Final_df.gRNA_combination.isin(ref_sgRNA_df.gRNA_combination)]
    sgRNA_groups = Final_filtered_df.groupby("gRNA_combination")
    temp_address_file = output_dir + '/Bartender_input_address'
    with open(temp_address_file, 'w') as file_a:
        for groups in sgRNA_groups:
            g, value = groups
            temp_name = output_dir + '/Clonal_barcode/' + g + '.bartender'
            file_a.write(temp_name + '\n')
            value[['Clonal_barcode', 'Read_ID']].to_csv(temp_name, sep=',', header=False, index=False)

if __name__ == "__main__":
    main()
