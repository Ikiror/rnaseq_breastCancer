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

#location of mappedReadsHisatLogFiles
WORKINGDIR="/data/users/aikiror/rnaseq/logfiles/mappedReadsHisatLogfiles_num7"

#storing the output
OUTPUTFILESDIR="/data/users/aikiror/rnaseq/outputFiles/multiqcFiles_num7_5"

mkdir -p $OUTPUTFILESDIR


apptainer exec --bind /data "$MULTIQCCONTAINER" multiqc -o "$OUTPUTFILESDIR" "$WORKINGDIR"
# module load MultiQC/1.11-foss-2021a 
# multiqc $READ1 $READ2 -o $OUTFILE 
# sleep 30

#to run: sbatch ../../rnaseqScripts/multiqc_num7_5/multiqc_num_7_5_slurm_array.sh