#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=create_sample_list_from_fastp
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=end,fail

#script will be run in /data/users/aikiror/rnaseq/outputReports

#this script should loop over each *1.fastq.gz file in the raw data folder; extract the name without the numbers and file name
#and echo the extracted name, the file name it looped over ending in 1.* and ending in 2.*

# Get the FASTP folder with the trimmed fastq files from the first argument
# /data/users/aikiror/rnaseq/outputFiles/fastpFiles_num3
FASTP_FOLDER=$1 

# Define the output file where the list will be saved
OUTPUT_FILE="/data/users/aikiror/rnaseq/rnaseqScripts/fastqc_num4/sampleListfromFastp_num4.tsv"


for FILE in "$FASTP_FOLDER"/*_fastp_R1.fastq.gz
do 
    SAMPLE="${FILE%_fastp_R1.fastq.gz}"  # Remove the _R1.fastq.gz part to get the base name
    READ1="${SAMPLE}_fastp_R1.fastq.gz"  # Construct full path for _R1 file
    READ2="${SAMPLE}_fastp_R2.fastq.gz"  # Construct full path for _R2 file

    echo -e "$(basename "$SAMPLE")\t$READ1\t$READ2" >> "$OUTPUT_FILE"
done

#to run: sbatch ../../rnaseqScripts/fastqc_num4/sampleList_num4.sh /data/users/aikiror/rnaseq/outputFiles/fastpFiles_num3
