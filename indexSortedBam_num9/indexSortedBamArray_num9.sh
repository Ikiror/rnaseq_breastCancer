#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=10:00:00
#SBATCH --mem=25gb
#SBATCH --cpus-per-task=16
#SBATCH --job-name=index_sorted_bam_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail


#steps -> index the coordinate sorted bam files using Samtools 
#      -> index the sorted bam files from the last step with samtools index

#run script in outputReports/indexSortedBam_num9 - to store generated error and output files
#sbatch ../../rnaseqScripts/indexSortedBam_num9/indexSortedBamArray_num9.sh

WORKDIR="/data/users/${USER}/rnaseq"
OUTPUTFILESDIR="$WORKDIR/outputFiles/indexedSortedBamFiles_num9"
LOGFILES="$WORKDIR/logfiles/indexSortedBamFiles_num9"

#make directories if they dont already exist
mkdir -p $LOGFILES
mkdir -p $OUTPUTFILESDIR 

#not used
#location of sorted bam files - *.bam
BAMFILES="${WORKDIR}/outputFiles/sortedBamFiles_num8"
# $BAMFILES/*.bam

SAMPLELIST="$WORKDIR/rnaseqScripts/indexSortedBam_num9/sampleList_num9.tsv"
#absolute path sample list to the .bam files: /data/users/aikiror/rnaseq/rnaseqScripts/indexSortedBam_num9/sampleList_num9.tsv

#Hisat samtool container
HISATCONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Extract sample information
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST)


apptainer exec --bind /data/ $HISATCONTAINER samtools index $READ1 $OUTPUTFILESDIR/index_${SAMPLE}.bai 2> $LOGFILES/${SAMPLE}_samtool_index_summary.log
# apptainer exec --bind /data/ $HISATCONTAINER samtools index $READ1 2> $LOGFILES/${SAMPLE}_samtool_index_summary.log

