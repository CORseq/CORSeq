#!/bin/bash

#PBS -N gene_quantification
#PBS -q smp
#PBS -P HEAL1390
#PBS -l select=1:ncpus=24
#PBS -l walltime=24:00:00


# SCRIPT: 3a.gene_transcript_quantification.sh
# AUTHOR: Denis Awany
# DATE: 11/02/2021
# REV: 1.0
# PURPOSE: This is a script to do gene quantification on the aligned reads. 
# 	   The HTSeq tool is utilized for the quantification. 
# PLATFORM: Linux



. /home/adenis/lustre3p/SATVI/RNA_Seq/config.txt # load the config file



mkdir -p ${wd}out_gene_transcript_quantification # directory to store output from gene counting.


for file in ${wd}out_alignment/*_Aligned.sortedByCoord.marked_duplicates.bam
do
	prefix=$(basename ${file} _Aligned.sortedByCoord.marked_duplicates.bam)

	htseq-count -f bam -i gene_id -s no \
	-m intersection-strict \
	-r name \
	-a 10 \
	$file ${wd}ref_data/gencode.v36.annotation.gtf > ${wd}out_gene_transcript_quantification/${prefix}_counts.txt
done


# We will have a text file with counts for every BAM we provided. Next, we merge these individual files to 
# make a single file. We use awk to do so.


awk '{arr[$1]=arr[$1]"\t"$2}END{for(i in arr)print i,arr[i]}' ${wd}out_gene_transcript_quantification/*_counts.txt >> ${wd}out_gene_transcript_quantification/merged_htseq_counts.tsv



# End of 3a.gene_transcript_quantification.sh
