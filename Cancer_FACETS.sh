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
list1="NS09N NS18N NS28N NS32N NS39N NS41N NS42N NS43N NS45N NS47N"
list2="NS09T NS18T NS28T NS32T NS39T NS41T NS42T NS43T NS45T NS47T"
echo $list1 | sed 's/ /\n/g' > /tmp/e.$$
echo $list2 | sed 's/ /\n/g' > /tmp/f.$$
paste /tmp/e.$$ /tmp/f.$$ | while read item1 item2; do
cnv_facets.R -t $BAM_DIR/$item2\_b37.bam -n $BAM_DIR/$item1\_b37.bam -vcf $DB/00-common_all.vcf.gz -bq 20 -mq 30 --cval 100 1000 -T $DB/SureselectV4.bed -N 8 -g hg19 -a $DB/TCGA_drivers.bed -o $FACET_DIR/$item2
done
rm /tmp/e.$$
rm /tmp/f.$$


