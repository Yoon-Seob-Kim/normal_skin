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

#for epidermis
list1="NS49N NS51N"
list2="NS49E NS51E"
#list1="NS02N NS04N NS05N NS06N NS07N NS08N NS09N NS10N NS11N NS12N NS13N NS14N NS16N NS17N NS18N NS19N NS20N NS21N NS22N NS25N NS26N NS27N NS28N NS29N NS30N NS31N NS32N NS33N NS37N NS38N NS39N NS40N NS41N NS42N NS43N NS45N NS47N NS48N NS49N NS51N"
#list2="NS02E NS04E NS05E NS06E NS07E NS08E NS09E NS10E NS11E NS12E NS13E NS14E NS16E NS17E NS18E NS19E NS20E NS21E NS22E NS25E NS26E NS27E NS28E NS29E NS30E NS31E NS32E NS33E NS37E NS38E NS39E NS40E NS41E NS42E NS43E NS45E NS47E NS48E NS49E NS51E"
echo $list1 | sed 's/ /\n/g' > /tmp/a.$$
echo $list2 | sed 's/ /\n/g' > /tmp/b.$$
paste /tmp/a.$$ /tmp/b.$$ | while read item1 item2; do
$GATK/gatk --java-options "-Xmx4g" Mutect2 -R $DB/human_g1k_v37.fasta -I $BAM_DIR/$item1\_b37.bam -I $BAM_DIR/$item2\_b37.bam -normal $item1 -tumor $item2 --intervals $DB/SureselectV4.list --germline-resource $DB/af-only-gnomad.raw.sites.b37.vcf.gz --f1r2-tar-gz $MT_DIR/$item2\.tar.gz -O $MT_DIR/$item2\.vcf.gz
$GATK/gatk --java-options "-Xmx4g" LearnReadOrientationModel -I $MT_DIR/$item2\.tar.gz -O $MT_DIR/$item2\_cal.tar.gz
$GATK/gatk --java-options "-Xmx4g" FilterMutectCalls -V $MT_DIR/$item2\.vcf.gz -ob-priors $MT_DIR/$item2\_cal.tar.gz -R $DB/human_g1k_v37.fasta --intervals $DB/SureselectV4.list  -O $MT_DIR/$item2\_filt.vcf
awk -F '\t' '{if($7== NULL) print; else if($7 == "FILTER") print ; else if($7 == "PASS") print}' $MT_DIR/$item2\_filt.vcf > $MT_DIR/$item2\_filt2.vcf
$GATK/gatk Funcotator --data-sources-path $FUNCTO -V $MT_DIR/$item2\_filt2.vcf --output $FUNCO_DIR/$item2\.txt --output-file-format VCF --ref-version hg19 -R $DB/human_g1k_v37.fasta --force-b37-to-hg19-reference-contig-conversion
done
