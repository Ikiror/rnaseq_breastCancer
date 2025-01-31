#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=index_sorted_bam_sample_list
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=end,fail

#script will be run in /data/users/aikiror/rnaseq/outputReports/indexSortedBam_num9/sampleListReport_num9

#this script should loop over each *.bam file in the sorted folder; extract the name without the _sorted.bam
#and echo the extracted name, and the file path 

#folder with .bam files
#/data/users/aikiror/rnaseq/outputFiles/sortedBamFiles_num8
BAM_FOLDER=$1 

#Output file where the list will be saved
#empty file created
OUTPUT_FILE="/data/users/aikiror/rnaseq/rnaseqScripts/indexSortedBam_num9/sampleList_num9.tsv"

for FILE in "$BAM_FOLDER"/*.bam
do 
    SAMPLE="${FILE%_sorted.bam}"  # Remove the _mapped.bam part to get the base name
    #READ1="${SAMPLE}_mapped.bam"  # Construct full path for _R1 file - redundant

    echo -e "$(basename "$SAMPLE")\t$FILE" >> "$OUTPUT_FILE"
done

#to run: sbatch /data/users/aikiror/rnaseq/rnaseqScripts/indexSortedBam_num9/sampleListScript_num9.sh /data/users/aikiror/rnaseq/outputFiles/sortedBamFiles_num8

