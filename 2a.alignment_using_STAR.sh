#!/bin/bash

#PBS -N alignment
#PBS -q smp
#PBS -P HEAL1390
#PBS -l select=1:ncpus=24
#PBS -l walltime=96:00:00
#PBS -m bae
#PBS -M awanydenis@gmail.com


# SCRIPT: 2a.alignment.sh
# AUTHOR: Denis Awany
# DATE: 09/02/2021
# REV: 1.0
# PLATFORM: Linux
# PURPOSE: This is a script to perform  reads mapping. The  'STAR' alignment tool is chosen. Mapping reads using STAR 
#	   is a two step process:(1) creating a genome index, and (2) mapping reads to the genome.





. /home/adenis/lustre3p/SATVI/RNA_Seq/config.txt # load the config file




STAR_soft=${wd}softwares/star/STAR-2.7.7a/bin/Linux_x86_64_static/STAR #the STAR tool
Filt_Fastq=${wd}out_quality_filtering/trimmomatic # directory to the filtered fastq files.



# Step 1: generating genome indexes files (done once)
# ---------------------------------------------
# (a) Download (and unzip) reference genome sequences (FASTA files) - the Genome sequence, primary assembly (GRCh38) fasta file
#wget -P ${wd}softwares/star/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz
#gunzip ${wd}ref_data/GRCh38.primary_assembly.genome.fa.gz

# (b) And download the corresponding annotations (GTF file), and unxip
#wget -P ${wd}softwares/star/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
#gunzip ${wd}ref_data/gencode.v36.annotation.gtf.gz


: '
mkdir -p ${wd}ref_data/genome_index_dir # directory to store genome indices.

${STAR_soft} --runMode genomeGenerate \
--genomeDir ${wd}ref_data/genome_index_dir \
--genomeFastaFiles ${wd}ref_data/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile ${wd}ref_data/gencode.v36.annotation.gtf \
--runThreadN 12 \
--genomeChrBinNbits 14 \
--genomeSAsparseD 2 \
--genomeSAindexNbases 12 



'

# Step 2 (final): mapping reads to the genome
# ------------------------------------------------

mkdir -p ${wd}out_alignment


# For all read pairs, we can use a loop.

for file in $(ls ${Filt_Fastq}/*.trimmed.fq.gz | sed -r 's/_[12].trimmed.fq.gz//' | uniq)
do
	filename=$(basename -- $file) 
 
	${STAR_soft} --runThreadN 8 \
	--genomeDir ${wd}ref_data/genome_index_dir \
	--sjdbGTFfile ${wd}ref_data/gencode.v36.annotation.gtf \
	--outFilterType BySJout \
	--readFilesIn ${file}_1.trimmed.fq.gz ${file}_2.trimmed.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix ${wd}out_alignment/${filename}_ \
	--outSAMtype BAM SortedByCoordinate Unsorted \
	--quantMode GeneCounts \
	--twopassMode Basic
done





# End of 2a.alignment.sh
