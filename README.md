
**CORseq/CORSeq** respository contains a workflow created by the SATVI (South African Tuberculosis Vaccine Initiative) Bioanalytical Node, for the analysis of data in the CORTIS study. 


**Analysis workflow:**
--------------------

**Scripts**. This contains all the scripts used in the analysis of the RNA-seq data. The .sh files contain PBS bash scripts that can be run directly on the terminal or, in the perferred case, submitted to the cluster. The scripts are labeled 1.,2.,3.,...according to the order in which they ought to be run, to proceed from quality inspection/control through reads alignment to gene/transcript quantification. While each script can be run individually, the 'chain_all_pbs_jobs.sh' scrip may be used to sequentially run all the scripts in a PBS job submission.

The 'config.txt' file is a reference file for all .sh scripts. This is the file where all required modules are loaded, paths to the datasets and files are specified, and parameter values are defined.

The the python script 'collapse_annotation.py' is utilized in the post alignment quality control step with the RNA-SeQC tool, where the script is used to collapse the full GTF annotation, which effectively produces the union of transcripts for each gene. This step is invoked prior to running RNA-SeQC tool because the current version of RNA-SeQC does not support full annotation GTFs.

**Data.** 

**Results.**
