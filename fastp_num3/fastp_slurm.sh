#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=fastp_slurm
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

# Define variables
WORKDIR="/data/users/aikiror/rnaseq"
OUTDIR="${WORKDIR}/outputFiles/fastpFiles_num3"
SAMPLELIST="${WORKDIR}/rnaseqScripts/absolutePathSampleList.tsv"

mkdir -p $OUTDIR

# Extract sample information
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST)
READ2=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST)

#for the analysis fastp version 0.23.2 was used
FASTPCONTAINER="/containers/apptainer/fastp_0.23.2--h5f740d0_3.sif"


#load the fastp module
#bind your specific folder for rerun
apptainer exec --bind /data ${FASTPCONTAINER} fastp \
    --dont_overwrite \
    -h "${OUTDIR}/${SAMPLE}_fastp.html" \
    --detect_adapter_for_pe \
    -i "${READ1}" \
    -I "${READ2}" \
    -o "${OUTDIR}/${SAMPLE}_fastp_R1.fastq.gz" \
    -O "${OUTDIR}/${SAMPLE}_fastp_R2.fastq.gz"

#to run: sbatch ../../rnaseqScripts/fastp_num3/fastp_slurm.sh