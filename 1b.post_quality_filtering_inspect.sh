#!/bin/bash

#PBS -N post_quality_filtering_inspect
#PBS -q smp
#PBS -P HEAL1390
#PBS -l select=1:ncpus=24
#PBS -l walltime=12:00:00
##PBS -m bae
##PBS -M awanydenis@gmail.com


# SCRIPT: 1b.post_quality_filtering_inspect.sh
# AUTHOR: Denis Awany
# DATE: 09/02/2021
# REV: 1.0
# PLATFORM: Linux
# PURPOSE: This is a script to do quality inspection on the quality filtered (trimmed) reads.
#          This is an optional step to assess the performance of quality filtering (trimming).



. /home/adenis/lustre3p/SATVI/RNA_Seq/config.txt  # Load the config.txt file




# Create a directories to store QC results from FastQC and MultiQC tools.
mkdir -p ${wd}out_post_quality_filtering_inspect/FastQC

mkdir -p ${wk}out_post_quality_filtering_inspect/MultiQC


# (a) Quality inspection using FastQC
# ------------------------------------------------------------------------

fastqc=${wk}softwares/FastQC/fastqc # the fastq software

${fastqc} -t 12 ${wd}out_quality_filtering/trimmomatic/*.trimmed.fq.gz -o ${wk}out_post_quality_filtering_inspect/FastQC #



# (b) Quality inspection using MultiQC
# -------------------------------------------------------------------------
# For this, we need to pass the .zip output from FastQC as input to MultiQC

multiqc -o ${wk}out_post_quality_filtering_inspect/MultiQC/ -n multiqc_report ${wk}out_post_quality_filtering_inspect/FastQC/*.zip





# End of 1b.post_quality_filtering_inspect.sh
