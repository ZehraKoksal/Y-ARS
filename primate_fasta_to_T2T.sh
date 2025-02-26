#!/bin/bash

#Example PanPan1:

#Bash code for aligning primate fasta sequence to T2T reference using minimap2
minimap2 -a T2T_Reference.fa mPanPan1-v2.0_bonoboY.fasta > mPanPan1-v2.0_t2t_mm2.sam
#Then compress the sam file to a bam file
samtools view -bS mPanPan1-v2.0_t2t_mm2.sam > mPanPan1-v2.0_t2t_mm2.bam
#Sort reads
samtools sort mPanPan1-v2.0_t2t_mm2.bam  -o mPanPan1-v2.0_t2t_mm2_sorted.bam
#Index bam file
samtools index mPanPan1-v2.0_t2t_mm2_sorted.bam

#Variant calling using quality thresholds and specifying sites in an additional .csv file after -R parameter
bcftools mpileup -Ou -q 50 -Q 20 -R 2520_sites.csv -f T2T_Reference.fa mPanPan1-v2.0_t2t_mm2_sorted.bam | \
bcftools call --ploidy 1 -O v -o mPanPan1-v2.0_t2t_mm2_sorted_t2t_sites.vcf -c
#Exclude indels from resulting vcf file
bcftools filter -e 'TYPE="INDEL"' mPanPan1-v2.0_t2t_mm2_sorted_t2t_sites.vcf -o mPanPan1-v2.0_t2t_mm2_sorted_t2t_sites_noindels.vcf
