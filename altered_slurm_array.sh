#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=16:00:00
#SBATCH --mem=16g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=slurm_array
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail


#added sbatch scripts to email me of beginning of execution, failure, and the end of the script.
#will run sbatch script in /data/users/aikiror/rnaseq/outputReports/fastqcOutputReports  - any error or slurm output files will be here

# define variables
WORKDIR="/data/users/aikiror/rnaseq"
REPORTDIR="$WORKDIR/outputReports/fastqcOutputReports"
SAMPLELIST="$WORKDIR/rnaseqScripts/sampleList.tsv"
FASTQCCONTAINER="/containers/apptainer/fastqc-0.12.1.sif"
OUTDIR="$WORKDIR/outputFiles/fastqcOutputFiles"

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

OUTFILE="$REPORTDIR/${SAMPLE}.txt"

############################


mkdir -p $REPORTDIR
mkdir -p $OUTDIR

echo "Starting FastQC for $SAMPLE at $(date)" > $OUTFILE

apptainer exec "$FASTQCCONTAINER" fastqc --outdir "$OUTDIR" "$READ1" "$READ2"

echo "Completed fastqc for $SAMPLE: outputs saved to $OUTDIR" >> "$OUTFILE"

sleep 30

# FASTQCPATH="/containers/apptainer/fastqc-0.12.1.sif"

# apptainer exec $FASTQCPATH fastqc --outdir . "$READS_DIR"/*.fastq

 #######

# #####inspiration from SLURM ibu documentation
# READS_DIR="$WORKDIR/breastCancerRawData/breastcancer_de/reads"

# apptainer exec fastqc_0.11.9--hdfd78af_1.sif fastqc --outdir . "$READS_DIR"/ecoli_*.fastq

# ###end inspiration

