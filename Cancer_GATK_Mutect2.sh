#!/bin/bash
#SBATCH -J NSvcall
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p mrcs
#SBATCH -w Zeta
#tool_path
GATK=/data/Delta_data4/kysbbubbu/tools/gatk-4.1.9.0
DB=/data/Delta_data4/kysbbubbu/genomicdb
FUNCTO=/data/Delta_data4/kysbbubbu/genomicdb/funcotator/funcotator_dataSources.v1.6.20190124s

##OUTPUT
BAM_DIR=./02_BAM
MT_DIR=./04_MT
ANNO_DIR=./05_Annovar
FUNCO_DIR=./05_FUNCTO

list1="NS09N NS18N NS28N NS32N NS39N NS41N NS42N NS43N NS45N NS47N"
list2="NS09T NS18T NS28T NS32T NS39T NS41T NS42T NS43T NS45T NS47T"
echo $list1 | sed 's/ /\n/g' > /tmp/a.$$
echo $list2 | sed 's/ /\n/g' > /tmp/b.$$
paste /tmp/a.$$ /tmp/b.$$ | while read item1 item2; do
$GATK/gatk --java-options "-Xmx4g" Mutect2 -R $DB/human_g1k_v37.fasta -I $BAM_DIR/$item1\_b37.bam -I $BAM_DIR/$item2\_b37.bam -normal $item1 -tumor $item2 --intervals $DB/SureselectV4.list --germline-resource $DB/af-only-gnomad.raw.sites.b37.vcf.gz --f1r2-tar-gz $MT_DIR/$item2\.tar.gz -O $MT_DIR/$item2\.vcf.gz
$GATK/gatk --java-options "-Xmx4g" LearnReadOrientationModel -I $MT_DIR/$item2\.tar.gz -O $MT_DIR/$item2\_cal.tar.gz
$GATK/gatk --java-options "-Xmx4g" CalculateContamination -I $MT_DIR/$item2\_pileups.table -matched $MT_DIR/$item1\_pileups.table -O $MT_DIR/$item2\_contamination.table
$GATK/gatk --java-options "-Xmx4g" FilterMutectCalls --intervals $DB/SureselectV4.list -V $MT_DIR/$item2\.vcf.gz --contamination-table $MT_DIR/$item2\_contamination.table -ob-priors $MT_DIR/$item2\_cal.tar.gz -R $DB/human_g1k_v37.fasta -O $MT_DIR/$item2\_filt.vcf
awk -F '\t' '{if($7== NULL) print; else if($7 == "FILTER") print ; else if($7 == "PASS") print}' $MT_DIR/$item2\_filt.vcf > $MT_DIR/$item2\_filt2.vcf
$GATK/gatk Funcotator --data-sources-path $FUNCTO -V $MT_DIR/$item2\_filt2.vcf --output $FUNCO_DIR/$item2\.txt --output-file-format VCF --ref-version hg19 -R $DB/human_g1k_v37.fasta --force-b37-to-hg19-reference-contig-conversion
done
rm /tmp/a.$$
rm /tmp/b.$$


