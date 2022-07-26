#!/bin/bash
#SBATCH -J FACET
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p mrcs
#SBATCH -w Zeta
#tool_path
DB=/data/Delta_data4/kysbbubbu/genomicdb

##OUTPUT
BAM_DIR=./02_BAM
FACET_DIR=./06_FACET
list1="NS02N NS04N NS05N NS06N NS07N NS08N NS09N NS10N NS11N NS12N NS13N NS14N NS16N NS17N NS18N NS19N NS20N NS21N NS22N NS25N NS26N NS27N NS28N NS29N NS30N NS31N NS32N NS33N NS37N NS38N NS39N NS40N NS41N NS42N NS43N NS45N NS47N NS49N NS51N"
list2="NS02E NS04E NS05E NS06E NS07E NS08E NS09E NS10E NS11E NS12E NS13E NS14E NS16E NS17E NS18E NS19E NS20E NS21E NS22E NS25E NS26E NS27E NS28E NS29E NS30E NS31E NS32E NS33E NS37E NS38E NS39E NS40E NS41E NS42E NS43E NS45E NS47E NS49E NS51E"
echo $list1 | sed 's/ /\n/g' > /tmp/e.$$
echo $list2 | sed 's/ /\n/g' > /tmp/f.$$
paste /tmp/e.$$ /tmp/f.$$ | while read item1 item2; do
cnv_facets.R -t $BAM_DIR/$item2\_b37.bam -n $BAM_DIR/$item1\_b37.bam -vcf $DB/00-common_all.vcf.gz -bq 20 -mq 30 -T $DB/SureselectV4.bed -N 8 -g hg19 -a $DB/TCGA_drivers.bed -o $FACET_DIR/$item2
done
rm /tmp/e.$$
rm /tmp/f.$$


