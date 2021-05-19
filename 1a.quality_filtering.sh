#!/bin/bash

#PBS -N quality_filtering
#PBS -q smp
#PBS -P HEAL1390
#PBS -l select=1:ncpus=24
#PBS -l walltime=12:00:00
##PBS -m bae
##PBS -M awanydenis@gmail.com



# SCRIPT: 1a.quality_filtering.sh
# AUTHOR: Denis Awany
# DATE: 09/02/2021
# REV: 1.0
# PLATFORM: Linux
# PURPOSE: This is a script to do quality (i) inspection, and (ii) filtering on the raw reads.  




# 1. QUALITY INSPECTION
#=====================================================================
# We shall use FastQC, MultiQC to inspect the quality of the raw reads.


. /home/adenis/lustre3p/SATVI/RNA_Seq/config.txt  # Load the config.txt file


fastq_files=$fastq_data  # the fastq files, passed from the config file.



# Create a directories to store QC results from FastQC and MultiQC tools.
mkdir -p ${wd}out_quality_filtering/FastQC

mkdir -p ${wk}out_quality_filtering/MultiQC


# 1(a) Quality inspection using FastQC
# ---------------------------------------------------------------------

fastqc=${wk}softwares/FastQC/fastqc # the fastq software

${fastqc} -t 12 ${fastq_files}*.fq.gz -o ${wk}out_quality_filtering/FastQC #


# 2(b) Quality inspection using MultiQC
# ---------------------------------------------------------------------
# For this, we need to pass the .zip output from FastQC as input to MultiQC
n
multiqc -o ${wk}out_quality_filtering/MultiQC/ -n multiqc_report ${wk}out_quality_filtering/fqC/*.zip






# 2. QUALITY FILTERING
#========================================================================

#  Quality filtering using Trimmomatic
# ------------------------------------------------------------------------

mkdir -p ${wd}out_quality_filtering/trimmomatic


# Use a loop to operate on all files:

for infile in ${fastq_files}*_1.fq.gz
do
	base=$(basename ${infile} _1.fq.gz)
	trimmomatic PE -phred33 \
	-threads 32 \
	${fastq_files}${base}_1.fq.gz ${fastq_files}${base}_2.fq.gz \
	${base}_1.trimmed.fq.gz ${base}_1.un.trimmed.fq.gz \
	${base}_2.trimmed.fq.gz ${base}_2.un.trimmed.fq.gz \
	ILLUMINACLIP:${wd}softwares/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done



# move the trimmd files to the out_quality_filtering directory
mv *trimmed.fq.gz ${wd}out_quality_filtering/trimmomatic


# create a director to move to the undesired 'orphaned' reads

mkdir -p ${wd}out_quality_filtering/trimmomatic/orphaned_reads

mv ${wd}out_quality_filtering/trimmomatic/*un.trimmed.fq.gz ${wd}out_quality_filtering/trimmomatic/orphaned_reads/





# End of 1a.quality_filtering.sh
