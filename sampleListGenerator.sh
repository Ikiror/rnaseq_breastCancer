#!/bin/bash

FASTQ_FOLDER=$1
#abs path - /data/users/aikiror/rnaseq/breastCancerRawData/breastcancer_de/reads
#relative path - rnaseq/breastCancerRawData/breastcancer_de/reads


for FILE in $FASTQ_FOLDER/*_*1.fastq.gz
do 
    PREFIX="${FILE%_*.fastq.gz}"
    SAMPLE=`basename $PREFIX`
    echo -e "${SAMPLE}\t$FILE\t${FILE%?.fastq.gz}2.fastq.gz" 
done

#code used in terminal:
#./sampleListGenerator.sh ../breastCancerRawData/breastcancer_de/reads > relativePathSampleList.tsv
#resulting file: /data/users/aikiror/rnaseq/rnaseqScripts/relativePathSampleList.tsv