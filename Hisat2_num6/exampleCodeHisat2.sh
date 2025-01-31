#code from lecture html 

#builds index files (as done in /data/users/aikiror/rnaseq/rnaseqScripts/Hisat2_num6/Index_Hisat.sh and as seen in /data/users/aikiror/rnaseq/outputFiles/genome_index)
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build \
    <reference_sequence_fasta> <index_basename>

where <reference_sequence_fasta> is data/users/${USER}/rnaseq/referenceGenomes/reference_file.fa
where <index_basename> is /data/users/${USER}/rnaseq/outputFiles/genome_index/reference_file

#code used :
apptainer exec --bind $(dirname ${REFERENCEGENOME}):/data,${OUTPUTFILESDIR}:/index ${HISATCONTAINER} \
    hisat2-build -p 16 /data/reference_file.fa /index/reference_file




#Maps Reads to the Reference Genome (Paired-End Reads)
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2 -x <index_basename> -1 <mate1.fastq.gz> -2 <mate2.fastq.gz> -S <mappedReads.sam> -p <threads>

apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2 -x /data/users/${USER}/rnaseq/outputFiles/genome_index/reference_file -1 <mate1.fastq.gz> -2 <mate2.fastq.gz> -S <mappedReads.sam> -p 16





apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools view -hbS <mappedReads.sam> > <mappedReads.bam>
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools sort -m <memory> -@ <threads> -o <sorted.bam> -T temp <mappedReads.bam>
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools index <sorted.bam>