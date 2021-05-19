#!/bin/bash

#PBS -N globin_align
#PBS -q smp
#PBS -P HEAL1390
#PBS -l select=1:ncpus=24
#PBS -l walltime=24:00:00


# SCRIPT: 2d.globin_removal.sh
# AUTHOR: Denis Awany
# DATE: 09/02/2021
# REV: 1.0
# PLATFORM: Linux
# PURPOSE: This is a script to to estimate proportion of globin reads for RNA-seq data. This allow us to then bioinformatically remove the globin reads.
#	   We achieve this by performing alignment to hgb reads.


. /home/adenis/lustre3p/SATVI/RNA_Seq/config.txt # load the config file



# Estimate globin rate by aligning to globin redas.

# create a directory to store the output.
mkdir -p ${wd}out_globin_align

# step 1: index the reference genome
bwa index -p ${wd}ref_data/globin_gene_sequences -a bwtsw ${wd}ref_data/globin_gene_sequences.fasta



# step 2: perform alignment
for file in $(ls ${wd}out_quality_filtering/trimmomatic/*.trimmed.fq.gz | sed -r 's/_[12].trimmed.fq.gz//' | uniq)
do
        filename=$(basename -- $file)

	bwa mem -M -t 32 ${wd}ref_data/globin_gene_sequences \
	${file}_1.trimmed.fq.gz ${file}_2.trimmed.fq.gz \
	2> ${wd}out_globin_align/${filename}_bwa.err \
	> ${wd}out_globin_align/${filename}.sam
done



# step 3: Alignment metrics
for file in ${wd}out_globin_align/*.sam
do
  	sampleName=$(basename ${file} .sam)
	
	samtools flagstat $file > ${wd}out_globin_align/${sampleName}_alignment_metrics.txt
done





# End of 2d.globin_removal.sh
