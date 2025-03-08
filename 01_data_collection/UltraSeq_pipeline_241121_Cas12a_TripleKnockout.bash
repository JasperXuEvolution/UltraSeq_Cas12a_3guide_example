#!/bin/bash
# The interpreter used to execute the script

# “#SBATCH” directives convey submission options:
#SBATCH --job-name=UltraSeq_pipeline_241121_Cas12a_TripleKnockout
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100g 
#SBATCH --time=24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch

# Instructions and checks:
# - Ensure the fastq file address file ends with a newline.
# - Verify the sgRNA reference file is correct.
# - Decide if sgRNA reference distance needs checking.
# - Check UltraSeq_Step3.py arguments (--a4 default is 4).

# Overall directory
LP="/labs/mwinslow/Haiqing/Ultra_seq/"

# Experiment-specific directories and variables
Input_experiment_ID="Analysis_241121_JH_HX_Cas12a_TripleKnockout"
Project_directory="${LP}${Input_experiment_ID}"
Python_script_address="${LP}${Input_experiment_ID}"

# Directories for input and output
Input_data_info_address="${Project_directory}/NGS_address"
Step1_address="${Project_directory}/Merging"
Step2_address="${Project_directory}/Bartender"
Step3_address="${Project_directory}/Processed_data"

# Load necessary modules
module load adapterremoval/2.3.1
source ~/miniconda3/etc/profile.d/conda.sh
conda activate UltraSeq

# Create directories for pipeline steps
mkdir -p "${Step1_address}" "${Step2_address}" "${Step3_address}"

# Main pipeline
while read -r line; do
   # Extract fields from the input file
   r1=$(echo "$line" | cut -d',' -f1)
   r2=$(echo "$line" | cut -d',' -f2)
   sampleID=$(echo "$line" | cut -d',' -f3)
   
   # Step 1: Sequence Merging
   temp_folder1="${Step1_address}/${sampleID}"
   mkdir -p "${temp_folder1}"
   AdapterRemoval --file1 "${r1}" --file2 "${r2}" \
       --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
       --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
       --basename "${temp_folder1}/Merged" --collapse --gzip
   echo "For sample ${sampleID}, sequence merging is complete."

   # Step 2: Extract gRNA array and barcode
   sample_folder="${Step2_address}/${sampleID}"
   clonal_barcode_folder="${sample_folder}/Clonal_barcode"
   mkdir -p "${sample_folder}" "${clonal_barcode_folder}"

   python3 "${Python_script_address}/UltraSeq_Step1_Cas12a_TripleArray.py" \
       --a "${Step1_address}/${sampleID}/Merged.collapsed.gz" \
       --b "${Project_directory}/sgRNA_Cas12a_combine.csv" \
       --o "${sample_folder}"

   # Step 3: Clonal barcode clustering
   bartender_input="${sample_folder}/Bartender_input_address"
   if [[ -f "${bartender_input}" ]]; then
       while read -r line2; do
           new_name="${line2/.bartender/}"
           bartender_single_com -z -1 -d 1 -l 3 -f "${line2}" -o "${new_name}" || {
               echo "Error: Bartender clustering failed for ${line2}."
               exit 1
           }
       done < "${bartender_input}"
   else
       echo "Warning: Bartender input file not found for sample ${sampleID}. Skipping clustering..."
   fi

   # Step 4: Process clustered data
   processed_sample_folder="${Step3_address}/${sampleID}"
   mkdir -p "${processed_sample_folder}"
   python3 "${Python_script_address}/UltraSeq_Step2_Cas12a_TripleArray.py" \
       --a "${sample_folder}" \
       --o "${processed_sample_folder}/" || {
           echo "Error: Step 2 failed for sample ${sampleID}."
           exit 1
       }

   # Clean up intermediate files
   rm -rf "${clonal_barcode_folder}" || {
       echo "Error: Failed to clean up ${clonal_barcode_folder}."
       exit 1
   }
done < "${Input_data_info_address}"

# Final Steps: Combine and process all data
python3 "${Python_script_address}/UltraSeq_Step3_Cas12a_TripleArray.py" \
    --o "${Step3_address}/" || {
        echo "Error: Step 3 failed."
        exit 1
    }

python3 "${Python_script_address}/UltraSeq_Step4_Cas12a_TripleArray.py" \
    --o "${Step2_address}/" || {
        echo "Error: Step 4 failed."
        exit 1
    }

# Job summary
sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j "$SLURM_JOBID"
