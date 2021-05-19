#!/bin/bash

#PBS -N index.bam
#PBS -q smp
#PBS -P HEAL1390
#PBS -l select=1:ncpus=24
#PBS -l walltime=48:00:00


# SCRIPT: 2b.visualize_align_with_igv.sh
# AUTHOR: Denis Awany
# DATE: 09/02/2021
# REV: 1.0
# PURPOSE: This is a script to create indexes on the bam files (from the alignnment step).
#	   The bam-indexed files can then be visualized on IGV.
#          This ia an optional step.
# PLATFORM: Linux



. /home/adenis/lustre3p/SATVI/RNA_Seq/config.txt



# The IGV software will be used to visualize the BAM files. For IGV to read the
# BAM files, the .bam files need to be indexed. We will use the samtools software:



for file in ${wd}out_alignment/*_Aligned.sortedByCoord.out.bam
do

	samtools index ${file}

done



# OTHER NOTES

: '
This creates a .bai file for each .bam file. Then download the *.bam, *.bai, .fa and .gtf files 
to the local computer. In the IGV (opened on a local computer), create your own genome database.
Click "Genomes"->"Create .genome" file. Fill out the following fields:
Unique identifier: 'SATVI' or other of your choice.
Descript name: 'SATVI' or other name.
Fasta: use the "Browse" button to find the .fa file
Gene file: use the "Browse" button to find the .gtf file
Then save the genome database on your computer. 

Finally, to inspect reghions using the IGV, from menu "File" -> "Load file", open the "filename_Aligned.sortedByCoord.out.bam".

'




# End of 2b.visualize_align_with_igv.sh
