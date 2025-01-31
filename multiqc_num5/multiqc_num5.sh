#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=16g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiqc_slurm
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail


#added sbatch scripts to email me of beginning of execution, failure, and the end of the script.

# define variables
MULTIQCCONTAINER="/containers/apptainer/multiqc-1.19.sif"

#location of fastqc files
#fastqc on data that we did fastp on
WORKINGDIR="/data/users/aikiror/rnaseq/outputFiles/fastqcFiles_num4"


#storing the output
OUTPUTFILESDIR="/data/users/aikiror/rnaseq/outputFiles/multiqcFiles_num5"




#SAMPLELIST="$WORKINGDIR/rnaseqScripts/relativePathSampleList.tsv"
#SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
#READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
#READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

#OUTFILE="$OUTPUTFILESDIR/${SAMPLE}"

############################

mkdir -p $OUTPUTFILESDIR
#mkdir -p $OUTFILE

#apptainer not working *should work now
# apptainer exec --bind /data "$FASTQCCONTAINER" fastqc --outdir "$OUTFILE" "$READ1" "$READ2"
#--bind makes the folder /data available to the fastqc container

# #echo "Run task for $SAMPLE with $READ1 and $READ2 ..." > $OUTFILE
 
apptainer exec --bind /data "$MULTIQCCONTAINER" multiqc -o "$OUTPUTFILESDIR" "$WORKINGDIR"
# module load MultiQC/1.11-foss-2021a 
# multiqc $READ1 $READ2 -o $OUTFILE 
# sleep 30

#to run: sbatch ../../rnaseqScripts/multiqc_num5/multiqc_num5.sh