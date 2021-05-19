#!/bin/bash

#PBS -N post_alignment_qc
#PBS -q smp
#PBS -P HEAL1390
#PBS -l select=1:ncpus=24
#PBS -l walltime=48:00:00


# SCRIPT: 2c.post_alignment.sh
# AUTHOR: Denis Awany
# DATE: 09/02/2021
# REV: 1.0
# PLATFORM: Linux
# PURPOSE: This is a script to to calculate post-alignment QC metrics. The RNA-SeQC (v2.4.0) tool is used,
#	   allowing us to determine QC metrics, and estimate the proportion of reads belonging to rRNA as 
#	   well as calculate proportion of reads belonging to hgbRNAs.


. /home/adenis/lustre3p/SATVI/RNA_Seq/config.txt # load the config file

mkdir -p ${wd}out_post_alignment_qc # directory to store output from the RNA-SeQC run.


# Current version of RNA-SeQC does not support full annotation gtfs. So we first use the python script
# 'collapse_annotation.py' to collapse our annotation, which effectively produces the union of transcripts for each gene.

python collapse_annotation.py ${wd}ref_data/gencode.v36.annotation.gtf ${wd}ref_data/gencode.v36.annotation_collapsed.gtf

# Index the reference fasta files (done once) to create .fai file
samtools faidx ${wd}ref_data/GRCh38.primary_assembly.genome.fasta


# Create sequence dictionary for the reference sequence
java -jar ${wd}softwares/picard.jar CreateSequenceDictionary.jar \
R=${wd}ref_data/GRCh38.primary_assembly.genome.fasta \
O=${wd}ref_data/GRCh38.primary_assembly.genome.dict



# Deduplicate pre-processed reads for RNASeqQC analysis. The Picard MarkDuplicates program identifies duplicates reads 
# in each dataset and retains only one of the duplicated reads for RNASeqQC.

for file in ${wd}out_alignment/*_Aligned.sortedByCoord.out.bam
do
	sampleName=$(basename ${file} _Aligned.sortedByCoord.out.bam)
	
	java -jar ${wd}softwares/picard.jar MarkDuplicates \
	I=$file \
	O=${wd}out_alignment/${sampleName}_Aligned.sortedByCoord.marked_duplicates.bam \
	M=${wd}out_alignment/${sampleName}.marked_dup_metrics.txt VALIDATION_STRINGENCY=SILENT
done



# Add read groups information to the BAM file (required by RNA-SeQC tool).

for file in ${wd}out_alignment/*_Aligned.sortedByCoord.marked_duplicates.bam
do
	sampleName=$(basename ${file} _Aligned.sortedByCoord.marked_duplicates.bam)

	java -jar  ${wd}softwares/picard.jar AddOrReplaceReadGroups \
      		I=$file \
      		O=${wd}out_alignment/${sampleName}_Aligned.sortedByCoord.marked_duplicates_withReadgroup.bam \
      		RGID=4 \
     		RGLB=lib1 \
      		RGPL=illumina \
      		RGPU=unit1 \
      		RGSM=20
done



# Index the bam files
for file in ${wd}out_alignment/*_Aligned.sortedByCoord.marked_duplicates_withReadgroup.bam
do

        samtools index ${file}

done



# Now perform post-alignment QC.

for file in ${wd}out_alignment/*_Aligned.sortedByCoord.marked_duplicates_withReadgroup.bam
do
  	sampleName=$(basename ${file} _Aligned.sortedByCoord.marked_duplicates_withReadgroup.bam)

        mkdir -p ${wd}out_post_alignment_qc/$sampleName

        echo -e "Sample ID\tBam File\tNotes \n$sampleName\t$file\t " >> ${wd}out_post_alignment_qc/${sampleName}/rnaseq_qc_sample.txt

        java -jar ${wd}softwares/RNA-SeQC_v1.1.8.jar -n 1000 \
                -t ${wd}ref_data/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf \
                -s "$sampleName|$file|NA" \
                -r ${wd}ref_data/GRCh38.primary_assembly.genome.fasta \
                -o ${wd}out_post_alignment_qc/$sampleName \
		-BWArRNA ${wd}ref_data/human_all_rRNA.fasta \
		-bwa /apps/chpc/bio/bwa/0.7.17/bwakit/bwa \
                -strat gc \
		-transcriptDetails \
                -gc ${wd}ref_data/gencode.v7.gc.txt
done





# End of 2c.post_alignment.sh
