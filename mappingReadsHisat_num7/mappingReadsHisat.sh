#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=10:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=16
#SBATCH --job-name=Hisat_mapping_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

#from lecture html
##steps -> For each sample separately, map the reads to the reference genome using Hisat2.
#       -> Convert the resulting sam files to bam format using Samtools


#run script in outputReports - to store generated error and output files
#sbatch ../../rnaseqScripts/mappingReadsHisat_num7/mappingReadsHisat.sh


WORKDIR="/data/users/${USER}/rnaseq"
OUTPUTFILESDIR="$WORKDIR/outputFiles/mappedReadsHisat_num7"
LOGFILES="$WORKDIR/logfiles/mappedReadsHisatLogfiles_num7"

SAMPLELIST="$WORKDIR/rnaseqScripts/fastqc_num4/sampleListfromFastp_num4.tsv"
#absolute path to the fastqc.gz files: /data/users/aikiror/rnaseq/rnaseqScripts/fastqc_num4/sampleListfromFastp_num4.tsv
GENOME_FILESDIR="$WORKDIR/outputFiles/genome_index" #reference genomes

HISATCONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
#REFERENCEGENOME="$WORKDIR/referenceGenomes/reference_file.fa" 

mkdir -p $OUTPUTFILESDIR
mkdir -p $LOGFILES

# Extract sample information
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST)
READ2=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST)

apptainer exec --bind /data /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "hisat2 -x \
    $GENOME_FILESDIR/reference_file -1 $READ1 -2 \
    $READ2 -S $OUTPUTFILESDIR/${SAMPLE}_mappedReads.sam -p 16 2> $LOGFILES/${SAMPLE}_hisat2_summary.log; \
    samtools view -S -b $OUTPUTFILESDIR/${SAMPLE}_mappedReads.sam > $OUTPUTFILESDIR/${SAMPLE}_mapped.bam"

#samtools

#view 
#With no options or regions specified, prints all alignments in the specified input alignment file 
#(in SAM, BAM, or CRAM format)to standard output in SAM format (with no header by default). You may specify one or 
#more space-separated region specifications after the input filename to restrict output to only those alignments 
#which overlap the specified region(s). Use of region specifications requires a coordinate-sorted and indexed input 
#file. Options exist to change the output format from SAM to BAM or CRAM, so this command also acts as a file format 
#conversion utility.

#-S -> not 100% necessary here. previous samtool versions needed this if the input was a sam file. now it's automatically detected
#-b -> converts to sam to bam format

