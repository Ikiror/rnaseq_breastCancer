#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=create_sample_list
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=end,fail

#script will be run in /data/users/aikiror/rnaseq/outputReports

#this script should loop over each *1.fastq.gz file in the raw data folder; extract the name without the numbers and file name
#and echo the extracted name, the file name it looped over ending in 1.* and ending in 2.*


####mostly original code copied from /data/courses/rnaseq_course/tools #####
# FASTQ_FOLDER=$1 #file where the raw data is - will be first argument
# OUTPUT_FILE=/data/users/aikiror/rnaseq/rnaseqScripts/sampleList.tsv #where I want the output to be directed to 

# for FILE in $FASTQ_FOLDER/*_*1.fastq.gz
# do 
#     PREFIX="${FILE%_*.fastq.gz}"
#     SAMPLE=`basename $PREFIX`
#     echo -e "${SAMPLE}\t$FILE\t${FILE%?.fastq.gz}2.fastq.gz" >> "$OUTPUT_FILE"
# done



# FASTQ_FOLDER=$1 #file where the raw data is - will be first argument
# OUTPUT_FILE=/data/users/aikiror/rnaseq/rnaseqScripts/sampleList.tsv #where I want the output to be directed to 

# for FILE in $FASTQ_FOLDER/*_*1.fastq.gz
# do 
#     PREFIX="${FILE%_R1.fastq.gz}"
#     SAMPLE=`basename $PREFIX`
#     echo -e "${SAMPLE}\t$FILE\t${PREFIX}_R2.fastq.gz" >> "$OUTPUT_FILE"
# done



# Get the FASTQ folder from the first argument
FASTQ_FOLDER=$1 

# Define the output file where the list will be saved
OUTPUT_FILE="/data/users/aikiror/rnaseq/rnaseqScripts/sampleList.tsv"

# # Loop over files with _1.fastq.gz suffix (R1)
# for FILE in $FASTQ_FOLDER/*_*1.fastq.gz
# do 
#     # Ensure the filename is processed correctly
#     PREFIX="${FILE%_1.fastq.gz}"

#     # Extract the sample name from the prefix
#     SAMPLE=$(basename $PREFIX)

#     # Get the R2 file by replacing _1.fastq.gz with _2.fastq.gz
#     R1_FILE=$(basename $FILE)
#     R2_FILE=$(basename ${FILE%_1.fastq.gz}_2.fastq.gz)

#     # Write the sample name and corresponding R1 and R2 files to the output file
#     echo -e "${SAMPLE}\t$R1_FILE\t$R2_FILE" >> "$OUTPUT_FILE"
# done

for FILE in "$FASTQ_FOLDER"/*_R1.fastq.gz
do 
    SAMPLE="${FILE%_R1.fastq.gz}"  # Remove the _R1.fastq.gz part to get the base name
    READ1="${SAMPLE}_R1.fastq.gz"  # Construct full path for _R1 file
    READ2="${SAMPLE}_R2.fastq.gz"  # Construct full path for _R2 file

    echo -e "$(basename "$SAMPLE")\t$READ1\t$READ2" >> "$OUTPUT_FILE"
done
