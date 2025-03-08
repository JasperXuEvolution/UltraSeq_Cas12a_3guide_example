#!/bin/bash
# 02-info_extraction_cas12a_triple_guide.bash
# This shell script is part of the UltraSeq pipeline for processing Cas12a triple guide data.
# It performs sequence merging, gRNA/barcode extraction, clustering, and data aggregation.
# Additionally, it cleans up intermediate files and outputs a summary log for quality control.

# ---------------------------------------------------------------------
# SBATCH directives: Job submission options for the SLURM scheduler.
# ---------------------------------------------------------------------
#SBATCH --job-name=UltraSeq_pipeline_cas12a_triple_guide
#SBATCH --mail-user=xhq@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100g 
#SBATCH --time=24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch

# ---------------------------------------------------------------------
# Environment Setup: Load configuration, modules, and activate Conda environment.
# ---------------------------------------------------------------------
source ../config.sh

# Load required modules and activate the Conda environment.
module load adapterremoval/2.3.1
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate UltraSeq

# ---------------------------------------------------------------------
# Directories and Input Files Setup
# ---------------------------------------------------------------------
working_dir="$PROJECT_DIR/01_data_collection"
input_data_info_address="${working_dir}/data/NGS_address"
# Use the appropriate guide reference file (uncomment the desired file).
guide_ref="$working_dir/data/guide_reference-cas12a_triple_KO.csv"

# Directory containing the Python scripts for the pipeline.
python_script_dir="$working_dir/main_code"

# Define directories for each pipeline step.
Step1_address="${working_dir}/data/Merging"
Step2_address="${working_dir}/data/Bartender"
Step3_address="${working_dir}/data/Processed_data"

# Create directories for pipeline outputs if they don't exist.
mkdir -p "${Step1_address}" "${Step2_address}" "${Step3_address}"

# ---------------------------------------------------------------------
# Main Pipeline: Process each sample listed in the input data info file.
# ---------------------------------------------------------------------
while read -r line; do
    # Extract fields from the input CSV file (assumed to be comma-separated):
    # r1: path to Read 1 FASTQ file
    # r2: path to Read 2 FASTQ file
    # sampleID: unique sample identifier
    r1=$(echo "$line" | cut -d',' -f1)
    r2=$(echo "$line" | cut -d',' -f2)
    sampleID=$(echo "$line" | cut -d',' -f3)
    
    echo "Processing sample: ${sampleID}"
    
    # ---------------------------
    # Step 1: Sequence Merging
    # ---------------------------
    temp_folder1="${Step1_address}/${sampleID}"
    mkdir -p "${temp_folder1}"
    AdapterRemoval --file1 "${r1}" --file2 "${r2}" \
       --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
       --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
       --basename "${temp_folder1}/Merged" --collapse --gzip
    echo "For sample ${sampleID}, sequence merging is complete."

    # ---------------------------
    # Step 2: gRNA Array and Barcode Extraction
    # ---------------------------
    sample_folder="${Step2_address}/${sampleID}"
    clonal_barcode_folder="${sample_folder}/Clonal_barcode"
    mkdir -p "${sample_folder}" "${clonal_barcode_folder}"

    # Execute the Python script to extract gRNA array and barcode information.
    python3 "${python_script_dir}/cas12a_triple_guide_parsing.py" \
       --a "${Step1_address}/${sampleID}/Merged.collapsed.gz" \
       --b "$guide_ref" \
       --o "${sample_folder}"

    # ---------------------------
    # Step 3: Clonal Barcode Clustering
    # ---------------------------
    bartender_input="${sample_folder}/Bartender_input_address"
    if [[ -f "${bartender_input}" ]]; then
        while read -r line2; do
            # Remove the ".bartender" extension to form the output file name.
            new_name="${line2/.bartender/}"
            # Run the bartender clustering command on the input file.
            bartender_single_com -z -1 -d 1 -l 5 -f "${line2}" -o "${new_name}" || {
                echo "Error: Bartender clustering failed for ${line2}."
                exit 1
            }
        done < "${bartender_input}"
    else
        echo "Warning: Bartender input file not found for sample ${sampleID}. Skipping clustering..."
    fi

    # ---------------------------
    # Step 4: Process Clustered Data
    # ---------------------------
    processed_sample_folder="${Step3_address}/${sampleID}"
    mkdir -p "${processed_sample_folder}"
    python3 "${python_script_dir}/cas12a_triple_guide_aggregate_barcode.py" \
       --a "${sample_folder}" \
       --o "${processed_sample_folder}/" || {
           echo "Error: Data aggregation failed for sample ${sampleID}."
           exit 1
       }

    # Clean up intermediate files (remove clonal barcode folder)
    rm -rf "${clonal_barcode_folder}" || {
        echo "Error: Failed to clean up ${clonal_barcode_folder}."
        exit 1
    }
done < "${input_data_info_address}"

# ---------------------------------------------------------------------
# Final Steps: Combine and Process All Data Across Samples
# ---------------------------------------------------------------------
python3 "${python_script_dir}/cas12a_triple_guide_aggregate_sample.py" \
    --o "${Step3_address}/" || {
        echo "Error: Final sample aggregation failed."
        exit 1
    }

python3 "${python_script_dir}/cas12a_triple_guide_aggregate_sample_for_QC.py" \
    --a "${Step2_address}/" --o "${Step3_address}/" || {
        echo "Error: QC aggregation failed."
        exit 1
    }

# ---------------------------------------------------------------------
# Job Summary: Display SLURM job statistics.
# ---------------------------------------------------------------------
sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j "$SLURM_JOBID"

# ---------------------------------------------------------------------
# Clean SLURM Output Log:
# Remove specific sections from the SLURM output log to produce a cleaner log file.
# The sed commands below remove lines between specific patterns.
# ---------------------------------------------------------------------
sed '/Loading barcodes from the file/,/The default or specified step is larger than the seed length../d' slurm-${SLURM_JOBID}.out | \
sed '/Trimming paired end reads/,/Processed a total ../d' > slurm-${SLURM_JOBID}_clean.out
