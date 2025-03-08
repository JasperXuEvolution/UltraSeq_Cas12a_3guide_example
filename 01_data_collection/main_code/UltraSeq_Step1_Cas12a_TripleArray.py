#!/usr/bin/env python
# coding: utf-8
#function
import gzip
import regex
import argparse
import sys
import pandas as pd

def find_reverse_complementary(input_string):
    temp_dic = {'A':'T','G':'C','T':'A','C':'G','N':'N',
                'a':'t','g':'c','t':'a','g':'c','n':'n'}
    return(''.join(tuple([temp_dic.get(x) for x in input_string[::-1]])))

def main():
    parser = argparse.ArgumentParser(description='A function to extract gRNA and clonal barcode from merged fastq gz file')
    parser.add_argument("--a", required=True, help="This is the input fastq gz file")
    parser.add_argument("--b", required=True, help="This is the input file of reference sgRNA")
    parser.add_argument("--o", required=True, help="This is the dir of output file")
    args = parser.parse_args()
    fastqgz_input_address = args.a
    ref_address = args.b
    output_dir = args.o

    ref_sgRNA_df = pd.read_csv(ref_address) # reference df

    if len(ref_sgRNA_df.DR1.unique())>1:
        sys.exit("Errors! More than one DR1 pattern")
    else:
        DR1 = ref_sgRNA_df.DR1.unique()[0]


    if len(ref_sgRNA_df.DR2.unique())>1:
        sys.exit("Errors! More than one DR2 pattern")
    else:
        DR2 = ref_sgRNA_df.DR2.unique()[0]

    if len(ref_sgRNA_df.DR3.unique())>1:
        sys.exit("Errors! More than one DR3 pattern")
    else:
        DR3 = ref_sgRNA_df.DR3.unique()[0]

    if len(ref_sgRNA_df.DR6.unique())>1:
        sys.exit("Errors! More than one DR6 pattern")
    else:
        DR6 = ref_sgRNA_df.DR6.unique()[0]

    temp_pattern1 = regex.compile('TAGTT' + '(.{16})' + 'TATGG' + DR1 + '(.{23})' + DR2 + '(.{23})' + DR3 + '(.{23})' + DR6)



    temp_total_read, temp_extracted_read, temp_matched_read = 0, 0, 0
    temp_sample_ID = output_dir.split('/')[-1]

    temp_gRNA1_list, temp_gRNA2_list, temp_gRNA3_list, temp_read_ID_list, temp_Clonal_barcode_list, temp_label_list = [], [], [], [], [], []
    with gzip.open(fastqgz_input_address,'rt') as handler1:
        temp_readID = handler1.readline().rstrip() # read ID
        temp_sequence1 = handler1.readline().rstrip()
        handler1.readline() # skip two lines
        handler1.readline()

        
        while temp_readID:
            # print("seq2 is "+temp_sequence2 +"/n")
            temp_total_read+=1
            # This is the regular expression pattern
            # 16 bp for barcode
            # 16-20 for sgRNA, some of the control sgRNA are shorter
            temp_search_result1 = temp_pattern1.search(temp_sequence1)
            if temp_search_result1:
                temp_clonal_barcode_r1 = temp_search_result1.group(1)
                temp_gRNA1_r1 = temp_search_result1.group(2)
                temp_gRNA2_r1 = temp_search_result1.group(3)
                temp_gRNA3_r1 = temp_search_result1.group(4)
                temp_extracted_read+=1
                temp_gRNA1_list.append(temp_gRNA1_r1)
                temp_gRNA2_list.append(temp_gRNA2_r1)
                temp_gRNA3_list.append(temp_gRNA3_r1)
                temp_read_ID_list.append(temp_readID)
                temp_Clonal_barcode_list.append(temp_clonal_barcode_r1)
                if (temp_gRNA1_r1 in ref_sgRNA_df.Guide1_sequence.values)&(temp_gRNA2_r1 in ref_sgRNA_df.Guide2_sequence.values)&(temp_gRNA3_r1 in ref_sgRNA_df.Guide3_sequence.values):
                    temp_matched_read+=1
                    temp_label_list.append('Expected')
                else:
                    temp_label_list.append('Unexpected')
            temp_readID = handler1.readline().rstrip() # read ID
            temp_sequence1 = handler1.readline().rstrip()
            handler1.readline() # skip two lines
            handler1.readline()
            
    print(f"Sample {temp_sample_ID:s} has totally {temp_total_read:d} reads. Among them {temp_extracted_read:d} has barcode and sgRNA, which is about {temp_extracted_read/temp_total_read:.3f}. {temp_matched_read:d} has expected sgRNA, which is about {temp_matched_read/temp_total_read:.3f}")
    Final_df = pd.DataFrame({'gRNA1': temp_gRNA1_list, 'gRNA2': temp_gRNA2_list,'gRNA3': temp_gRNA3_list,
                                   'Clonal_barcode':temp_Clonal_barcode_list, 'Read_ID': temp_read_ID_list,
                                   'Sample_ID':temp_sample_ID,'Class':temp_label_list})
    Final_df['gRNA_combination'] = Final_df['gRNA1'] + '_' + Final_df['gRNA2'] + '_' + Final_df['gRNA3']

    Final_df_S = Final_df.groupby([col for col in Final_df.columns if col != 'Read_ID']).size().reset_index(name='Count')
    Final_df.to_csv(f"{output_dir}/Intermediate_df.csv", index=False)
    Final_df_S.to_csv(f"{output_dir}/Combined_deduplexed_df_full.csv", index=False)

    Final_filtered_df = Final_df[Final_df.gRNA_combination.isin(ref_sgRNA_df.gRNA_combination)] # I only take those array exist.
    sgRNA_groups = Final_filtered_df.groupby("gRNA_combination")
    temp_address_file = output_dir+'/Bartender_input_address'
    file_a = open(temp_address_file, 'w')
    for groups in sgRNA_groups:
        g, value = groups
        temp_name = output_dir + '/Clonal_barcode/' + g + '.bartender'
        file_a.write(temp_name + '\n')
        value[['Clonal_barcode', 'Read_ID']].to_csv(temp_name, sep=',', header=False, index=False)
    file_a.close()

if __name__ == "__main__":
    main()  


