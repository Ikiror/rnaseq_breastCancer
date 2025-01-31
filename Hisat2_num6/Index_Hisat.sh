#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=16
#SBATCH --job-name=Hisat_index_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

#run script in outputReports - to store generated error and output files


# get ref genome and rename
# wget ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#used: curl https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > reference_file.fa.gz


# check sum of acquired file (Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)
# sum reference_file.fa.gz

# get annotation genome and rename
# wget ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
#used: curl https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz > annotation_file.gtf.gz

# check sum (Homo_sapiens.GRCh38.113.gtf.gz)
# sum annotation_file.gtf.gz

# unzip 
#gzip -d annotation_file.gtf.gz 
#gzip -d reference_file.fa.gz 

WORKDIR="/data/users/${USER}/rnaseq/"
OUTPUTFILESDIR="$WORKDIR/outputFiles/genome_index"
REFERENCEGENOME="$WORKDIR/referenceGenomes/reference_file.fa"
HISATCONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

mkdir -p $OUTPUTFILESDIR


# apptainer exec --bind /data HISATCONTAINER hisat2-build -p 16 genome.fa genome

# apptainer exec --bind /data HISATCONTAINER -p 16 /data/users/aikiror/rnaseq/referenceGenomes/sequence_file.fa genome_index/genome_index


#The --bind /data/users/${USER}/rnaseq/referenceGenomes:/data part of the command binds host machine’s 
#/data/users/${USER}/rnaseq/referenceGenomes directory to /data inside the container.

#The --bind /data/users/${USER}/rnaseq/outputFiles/genome_index:/index part binds the host directory 
#/data/users/${USER}/rnaseq/outputFiles/genome_index to /index inside the container.

apptainer exec --bind $(dirname ${REFERENCEGENOME}):/data,${OUTPUTFILESDIR}:/index ${HISATCONTAINER} \
    hisat2-build -p 16 /data/reference_file.fa /index/reference_file

#to run: sbatch ../../rnaseqScripts/Hisat2_num6/Index_Hisat.sh