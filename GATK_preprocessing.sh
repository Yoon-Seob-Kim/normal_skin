#!/bin/bash
#SBATCH -J NS
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o %x%j.o
#SBATCH -e %x%j.e
#SBATCH -p mrcs
#SBATCH -w Zeta

#tool_path
BWA=/data/Delta_data4/kysbbubbu/tools/bwa-0.7.17
GATK=/data/Delta_data4/kysbbubbu/tools/gatk-4.1.9.0
QUALIMAP=/data/Delta_data4/kysbbubbu/tools/qualimap_v2.2.1
DB=/data/Delta_data4/kysbbubbu/genomicdb
FASTQC=/data/Delta_data4/kysbbubbu/tools/FastQC_v0.11.9
FASTQ_DIR=/data/Delta_data4/kysbbubbu/Project_Rawdata/NormalSkin
##OUTPUT
FASTQC_DIR=./01_FastQC
BAM_DIR=./02_BAM
BAMQC_DIR=./03_BAMQC
MT_DIR=./04_MT
for i in NS49N NS49E NS51N NS51E
do
#$FASTQC/fastqc -o $FASTQC_DIR -t 4 $FASTQ_DIR/$i\_1.fq.gz $FASTQ_DIR/$i\_2.fq.gz
#$BWA/bwa mem -t 4 $DB/human_g1k_v37.fasta $FASTQ_DIR/$i\_1.fq.gz $FASTQ_DIR/$i\_2.fq.gz > $i\.sam
#$GATK/gatk --java-options "-Xmx4g" AddOrReplaceReadGroups -I $i\.sam -O $i\_RG.bam -SO coordinate --CREATE_INDEX -ID $i -LB $i -PU $i -PL ILLUMINA -SM $i
#rm $i\.sam
#$GATK/gatk --java-options "-Xmx4g" MarkDuplicates -I $i\_RG.bam -O $i\.bam -M $BAMQC_DIR/$i\_dup_metrics.txt
#rm $i\_RG.bam
#$GATK/gatk --java-options "-Xmx4g" BaseRecalibrator -R $DB/human_g1k_v37.fasta -I $i\.bam --known-sites $DB/dbsnp_138.b37.vcf --known-sites $DB/1000G_phase1.indels.b37.vcf -O $i\_recal_data.table
#$GATK/gatk --java-options "-Xmx4g" ApplyBQSR -R $DB/human_g1k_v37.fasta -I $i\.bam -bqsr $i\_recal_data.table -O $BAM_DIR/$i\_b37.bam
#$GATK/gatk --java-options "-Xmx4g" AnalyzeCovariates -bqsr $i\_recal_data.table -plots $BAMQC_DIR/$i\_.pdf
#rm $i\.bam $i\.bam.bai $i\.bam.sbi $i\_recal_data.table
#$GATK/gatk --java-options "-Xmx4g" GetPileupSummaries -R $DB/human_g1k_v37.fasta -I $BAM_DIR/$i\_b37.bam -O $MT_DIR/$i\_pileups.table -V $DB/af-only-gnomad.raw.sites.b37.vcf.gz --intervals $DB/SureselectV4.list
$QUALIMAP/qualimap bamqc -bam $BAM_DIR/$i\_b37.bam -gff $DB/SureselectV4.bed -outdir $BAMQC_DIR/$i -c --java-mem-size=4G
#$GATK/gatk --java-options "-Xmx4g" HaplotypeCaller -R $DB/human_g1k_v37.fasta -I $BAM_DIR/$i\_b37.bam --intervals $DB/SureselectV4.list -ERC GVCF -O $MT_DIR/$i\_germline.vcf.gz
done

