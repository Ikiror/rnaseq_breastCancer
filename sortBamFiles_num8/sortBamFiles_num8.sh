#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=10:00:00
#SBATCH --mem=25gb
#SBATCH --cpus-per-task=16
#SBATCH --job-name=sort_bam_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail


#steps -> sort the bam files by genomic coordinates using Samtools

#run script in outputReports - to store generated error and output files
#sbatch ../../rnaseqScripts/sortBamFiles_num8/sortBamFiles_num8.sh

WORKDIR="/data/users/${USER}/rnaseq"
OUTPUTFILESDIR="$WORKDIR/outputFiles/sortedBamFiles_num8"
LOGFILES="$WORKDIR/logfiles/sortedBamFiles_num8"
TEMPFILES="$OUTPUTFILESDIR/tempFiles"

#make directories if they dont already exist
mkdir -p $LOGFILES
mkdir -p $OUTPUTFILESDIR
mkdir -p $TEMPFILES

#not used
#location of bam files - *.bam
BAMFILES="${WORKDIR}/outputFiles/mappedReadsHisat_num7"
# $BAMFILES/*.bam

SAMPLELIST="$WORKDIR/rnaseqScripts/sortBamFiles_num8/sampleList_num8.sh"
#absolute path to the .bam files: /data/users/aikiror/rnaseq/rnaseqScripts/sortBamFiles_num8/sampleList_num8.sh

#Hisat samtook container
HISATCONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Extract sample information
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST)


#apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    #samtools sort -m <memory> -@ <threads> -o <sorted.bam> -T temp <mappedReads.bam>

#memory is less (16gb) than specified above (20gb) in the #SBATCH specifications to allow for overhead
apptainer exec --bind /data/ $HISATCONTAINER samtools sort -m 2G -@ 4 -o $OUTPUTFILESDIR/${SAMPLE}_sorted.bam -T $TEMPFILES/temp_${SAMPLE}_${SLURM_JOB_ID} $READ1 2> $LOGFILES/${SAMPLE}_samtool_sort_summary.log



# apptainer exec --bind /data /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "hisat2 -x \
#     $GENOME_FILESDIR/reference_file -1 $READ1 -2 \
#     $READ2 -S $OUTPUTFILESDIR/${SAMPLE}_mappedReads.sam -p 16 2> $LOGFILES/${SAMPLE}_hisat2_summary.log; \
#     samtools sort -S -b $OUTPUTFILESDIR/${SAMPLE}_mappedReads.sam > $OUTPUTFILESDIR/${SAMPLE}_mapped.bam"

#samtools

#sort 


