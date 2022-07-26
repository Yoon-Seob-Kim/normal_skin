#!/bin/bash
#SBATCH -J NSvcall
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p mrcs
#SBATCH -w Zeta
#tool_path
BWA=/data/MRC1_data4/kysbbubbu/tools/bwa-0.7.17
GATK=/data/MRC1_data4/kysbbubbu/tools/gatk-4.1.9.0
DB=/data/MRC1_data4/kysbbubbu/genomicdb

##OUTPUT
BAM_DIR=./02_BAM
MT_DIR=./04_MT
MUT_DIR=./08_mutsigcv
item2=filtered_mutation
$GATK/gatk Funcotator -R $DB/hg19/ucsc.hg19.fasta -V $MUT_DIR/$item2\.vcf  -O $MUT_DIR/$item2\.maf --output-file-format MAF --data-sources-path $DB/funcotator/funcotator_dataSources.v1.7.20200521s --ref-version hg19
