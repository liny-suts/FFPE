#!/bin/bash
index_dir=grch38_snp_tran
index_name=genome_snp_tran
trim_dir='fastq file directory'
sam_dir='directory to output sam file'
bam_dir='directory to output bam file'
count_dir='directory to output featurecounts results'
annotation_f=Homo_sapiens.GRCh38.104.coding.gtf
samplelist='sample id list'
bam_dir='directory to output bam file'
read_distribution_dir='directory to output read_distribution results'
##change your options above
mkdir -p ${sam_dir}
mkdir -p ${bam_dir}
mkdir -p ${count_dir}
mkdir -p ${read_distribution_dir}

while read sample_id
do
##hisat2 - featurecounts
hisat2 -p 20 -x ${index_dir}/${index_name} -1 ${trim_dir}/${sample_id}_R1_trim.fastq.gz -2 ${trim_dir}/${sample_id}_R2_trim.fastq.gz -S ${sam_dir}/${sample_id}.sam 2>${sam_dir}/${sample_id}.maplog
featureCounts -T 16 -p --countReadPairs -t exon -g gene_id -a ${annotation_f} -o ${count_dir}/${sample_id}.counts ${sam_dir}/${sample_id}.sam 2>${count_dir}/${sample_id}.featureCountlog
samtools view -S -b ${sam_dir}/${sample_id}.sam -@ 8 | samtools sort -o ${sample_id}_sort.bam -@ 8
samtools index -b ${sample_id}_sort.bam

##stringtie
stringtie -p 16 -eB -G ${annotation_f} -A ./abundance/${bam}.gene_abund.tab ${prefix}/${bam}.bam

##get read distribution
read_distribution.py -i ${bamdir}/${sample_id}_sort.bam -r Homo_sapiens.GRCh38.104.chr.bed>${read_distribution_dir}/dist_${sample_id}

done<${samplelist}