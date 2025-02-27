### Processing short read sequencing data: From fastq files to vcf files
#EXAMPLE: GRCh37 (repeat for GRCh38, T2T, Y-ARS)

##TRIMMING STEP
java -jar trimmomatic-0.39.jar PE -threads 6 ./sample_name_1.fastq.gz ./sample_name_2.fastq.gz ./sample_name_1_Paired.fastq.gz  ./sample_name_1_UnPaired.fastq.gz ./sample_name_2_Paired.fastq.gz ./sample_name_2_UnPaired.fastq.gz ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36


##ALIGNING TO REFERENCE 
bwa_v07.7.17 mem -t 24 ./GRCh37.fna ./sample_name_1_Paired.fastq.gz ./sample_name_2_Paired.fastq.gz > ./sample_name_paired_grch37.bam
bwa_v07.7.17 mem -t 24 ./GRCh37.fna ./sample_name_1_UnPaired.fastq.gz > ./sample_name_1_unpaired_grch37.bam
bwa_v07.7.17 mem -t 24 ./GRCh37.fna ./sample_name_2_UnPaired.fastq.gz > ./sample_name_2_unpaired_grch37.bam


##SORTING READS
samtools sort ./sample_name_paired_grch37.bam -o ./sample_name_paired_grch37_sorted.bam
samtools sort ./sample_name_1_unpaired_grch37.bam -o ./sample_name_1_unpaired_grch37_sorted.bam
samtools sort ./sample_name_2_unpaired_grch37.bam -o ./sample_name_2_unpaired_grch37_sorted.bam


##MERGE AND SORT
samtools merge -@24 -r ./sample_merged_grch37.bam ./sample_name_paired_grch37_sorted.bam ./sample_name_1_unpaired_grch37_sorted.bam ./sample_name_2_unpaired_grch37_sorted.bam
samtools sort ./sample_merged_grch37.bam -o ./sample_merged_sorted_grch37.bam


##ADD READ GROUP
samtools view -H ./sample_merged_sorted_grch37.bam | sed '/^@RG/s/$/\tSM:1KGP\tLB:None\tPL:Illumina/' | samtools reheader - ./sample_merged_sorted_grch37.bam > ./sample_merged_sorted_RG_grch37.bam


##MARK DUPLICATES
java -jar ./picard.jar MarkDuplicates -I ./sample_merged_sorted_RG_grch37.bam -O ./sample_merged_sorted_RG_duprmv_grch37.bam -M File_to_write_duplication_metrics_to.txt --REMOVE_DUPLICATES true


##INDEX BAM FILE
samtools index ./sample_merged_sorted_RG_duprmv_grch37.bam


##BASE QUALITY SCORE RECALIBRATION
gatk BaseRecalibrator -I ./GRCh37.fna --known-sites ./00-All-SNPs_Yfiltered.recode.chrnames.vcf.gz  -O ./sample_grch37.recaldata.table
gatk ApplyBQSR -R ./GRCh37.fna -I ./sample_merged_sorted_RG_duprmv_grch37.bam --bqsr-recal-file ./sample_grch37.recaldata.table -O ./sample_merged_sorted_RG_duprmv_BQSR_grch37.bam


##VARIANT CALLING
bcftools mpileup -Ou --skip-indels -q 50 -Q 20 -r NC_000024.9 -f ./GRCh37.fna ./sample_merged_sorted_RG_duprmv_BQSR_grch37.bam | bcftools call -v --ploidy 1 -O v -o ./sample_grch37.vcf -c
bcftools filter -i 'INFO/DP >= 2' ./sample_grch37.vcf -O v -o ./sample_grch37_rd2.vcf
