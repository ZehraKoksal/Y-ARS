#

#Example: GRCh37 aligned to GRCh38
#alignment
minimap2 -a ./GCF_000001405.25_GRCh37.p13_genomic.fna ./GCF_000001405.40_GRCh38.p14_genomic.fna > ./grch38_to_grch37_mm2.sam
#compress sam to binary bam format
samtools view -bS ./grch38_to_grch37_mm2.sam > ./grch38_to_grch37_mm2.bam
#sort reads
samtools sort ./grch38_to_grch37_mm2.bam  -o ./grch38_to_grch37_mm2_sorted.bam
#indem bam file (creating bam.bai)
samtools index ./grch38_to_grch37_mm2_sorted.bam
#call variants on Y chromosome of ploidy 1
bcftools mpileup -Ou -q 50 -Q 20 -r NC_000024.9 -f ./GCF_000001405.25_GRCh37.p13_genomic.fna ./grch38_to_grch37_mm2_sorted.bam | \
bcftools call --ploidy 1 -O v -o ./grch38_to_grch37.vcf -c
#remove indels
bcftools filter -e 'TYPE="INDEL"' ./grch38_to_grch37.vcf -o ./grch38_to_grch37_noindels.vcf

#Remaining fasta files GRCh38, T2T, Y-ARS were aligned to GRCh37
