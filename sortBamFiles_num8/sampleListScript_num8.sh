#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bam_sort_sample_list
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=end,fail

#script will be run in /data/users/aikiror/rnaseq/outputReports/sortBamReports_num8/sampleListReport_num8

#this script should loop over each *.bam file in the raw data folder; extract the name without the _mapped.bam
#and echo the extracted name, and the file path 

#folder with .bam files
#/data/users/aikiror/rnaseq/outputFiles/mappedReadsHisat_num7
BAM_FOLDER=$1 

#Output file where the list will be saved
#empty file created
OUTPUT_FILE="/data/users/aikiror/rnaseq/rnaseqScripts/sortBamFiles_num8/sampleList_num8.sh"

for FILE in "$BAM_FOLDER"/*.bam
do 
    SAMPLE="${FILE%_mapped.bam}"  # Remove the _mapped.bam part to get the base name
    #READ1="${SAMPLE}_mapped.bam"  # Construct full path for _R1 file - redundant

    echo -e "$(basename "$SAMPLE")\t$FILE" >> "$OUTPUT_FILE"
done

#to run: sbatch /data/users/aikiror/rnaseq/rnaseqScripts/sortBamFiles_num8/sampleListScript_num8.sh /data/users/aikiror/rnaseq/outputFiles/mappedReadsHisat_num7
#../../../rnaseqScripts/sortBamFiles_num8/sampleListScript_num8.sh
#../../../outputFiles/mappedReadsHisat_num7
#/data/users/aikiror/rnaseq/outputFiles/mappedReadsHisat_num7
