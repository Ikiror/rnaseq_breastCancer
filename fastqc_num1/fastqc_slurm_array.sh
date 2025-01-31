#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=10:00:00
#SBATCH --mem=16g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=slurm_array
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

#added sbatch scripts to email me of beginning of execution, failure, and the end of the script.

# define variables
# FASTQCCONTAINER="/containers/apptainer/fastqc-0.12.1.sif"

WORKINGDIR="/data/users/aikiror/rnaseq"
OUTPUTFILESDIR="$WORKINGDIR/outputFiles/fastqcFiles_num1"
SAMPLELIST="$WORKINGDIR/rnaseqScripts/absolutePathSampleList.tsv"


SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

OUTFILE="$OUTPUTFILESDIR/${SAMPLE}"

############################

mkdir -p $OUTPUTFILESDIR
mkdir -p $OUTFILE

#apptainer not working *should work now - wasn't working because original apptainer code did  not bind /data
# apptainer exec --bind /data "$FASTQCCONTAINER" fastqc --outdir "$OUTFILE" "$READ1" "$READ2"
#--bind makes the folder /data available to the fastqc container

# #echo "Run task for $SAMPLE with $READ1 and $READ2 ..." > $OUTFILE

module load FastQC/0.11.9-Java-11  

fastqc $READ1 $READ2 -o $OUTFILE 
sleep 30
