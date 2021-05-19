#!/bin/bash

#PBS -N alignment
#PBS -q smp
#PBS -P HEAL1390
#PBS -l select=1:ncpus=24
#PBS -l walltime=12:00:00
#PBS -m bae
#PBS -M awanydenis@gmail.com


# SCRIPT: 2a.alignment_using_GSNAP.sh
# AUTHOR: Denis Awany
# DATE: 09/02/2021
# REV: 1.0
# PLATFORM: Linux
# PURPOSE: This is a script to perform  reads mapping. The script performs alignment using 'GSNAP'. 
#          Mapping reads using GSNAP is a two step process:(1) creating a genome index, and 
#          (2) mapping reads to the genome.





. /home/adenis/lustre3p/SATVI/RNA_Seq/config.txt # load the config file



Filt_Fastq=${wd}out_quality_filtering/trimmomatic # directory to the filtered fastq files.



# Step 1: generating genome indexes.
# ---------------------------------------------
mkdir -p ${wd}ref_data/genome_index_dir # directory to store genome indices.


gmap_build -d ${wd}ref_data/genome_index_dir ${wd}ref_data/GRCh38.primary_assembly.genome.fasta




# Step 2 (final): mapping reads to the genome
# --------------------------------------------


mkdir -p ${wd}out_alignment


# For all read pairs, we can use a loop.

for file in $(ls ${Filt_Fastq}/*.trimmed.fq.gz | sed -r 's/_[12].trimmed.fq.gz//' | uniq)
do
	filename=$(basename -- $file)

	gsnap -d grch38_chr1 -D ${wd}ref_data/genome_index_dir \
	--gunzip \
	-t 6 -M 2 -n 10 -N 1 \
	--quality-protocol=sanger \
	-w 200000 --pairmax-rna=200000 \
	-E 1 -B 2 \
	-A sam ${file}_1.trimmed.fq.gz ${file}_2.trimmed.fq.gz | \
	samtools view -bS - | samtools sort - \
	> ${wd}out_alignment/${filename}_Aligned.sortedByCoord.out.bam
done





# End of 2a.alignment.sh
