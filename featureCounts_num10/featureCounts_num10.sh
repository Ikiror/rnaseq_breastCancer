#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=30g
#SBATCH --cpus-per-task=6
#SBATCH --job-name=feature_counts
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

#loaded subread image  
#contains featureCounts
SUBREAD_IMAGE="/containers/apptainer/subread_2.0.1--hed695b0_0.sif"

WORKDIR="/data/users/${USER}/rnaseq"
OUTPUTFILESDIR="$WORKDIR/outputFiles/featureCounts_num10"
LOGFILES="$WORKDIR/logfiles/featureCounts_num10"

#make directories if they dont already exist
mkdir -p $LOGFILES
mkdir -p $OUTPUTFILESDIR 

#location of sorted bam files - *.bam
BAMFILES="${WORKDIR}/outputFiles/sortedBamFiles_num8"
# $BAMFILES/*.bam  - will iterate over all bam files

#reference_genomes
#annotation file
ANNOTATION_FILE="/data/users/aikiror/rnaseq/referenceGenomes/annotation_file.gtf"
#reference file
REFERENCE_FILE="/data/users/aikiror/rnaseq/referenceGenomes/reference_file.fa"

apptainer exec --bind /data/ $SUBREAD_IMAGE featureCounts -Q 10 -C -p -t exon -g gene_id -T 6 -a $ANNOTATION_FILE -G $REFERENCE_FILE -o $OUTPUTFILESDIR/featureCounts_num10_${SLURM_JOB_ID}.txt ${BAMFILES}/*.bam 2> $LOGFILES/featureCounts_num10_${SLURM_JOB_ID}_summary.log
# apptainer exec --bind /data/ $HISATCONTAINER samtools index $READ1 2> $LOGFILES/${SAMPLE}_samtool_index_summary.log
# apptainer exec --bind /data/ $HISATCONTAINER samtools index $READ1 2> $LOGFILES/${SAMPLE}_samtool_index_summary.log

#-Q = The minimum mapping quality score a read must satisfy in order to be counted
#-p = indicates that paired reads will be used
#-t = feature type to be counted (used exon)
#-g = category to use to group (used gene_id)
#-T = threads to be used (using 5)
#-a = annotation file .gtf format
#-G = reference file .fa format
#-o = output file name chosen
#-C = removes multi-mapped reads (align to multiple regions) and only leaves uniquely aligned reads; used to help eliminate noise

#resources:
#https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html
#https://subread.sourceforge.net/SubreadUsersGuide.pdf

#What proportion of reads overlaps with annotated genes in each sample?
#How many reads, on average, are unassigned due to ambiguity? 
#Can you think of situation when it may not be possible to assign a read unambiguously to a particular gene?