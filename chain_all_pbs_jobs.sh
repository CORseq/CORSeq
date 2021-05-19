#!/bin/bash
first=$(qsub 1a.quality_filtering.pbs)
echo $first
second=$(qsub -W depend=afterok:$first 1b.post_quality_filtering_inspect.pbs)
echo $second
third=$(qsub -W depend=afterok:$second 2a.alignment_using_GSNAP.pbs)
echo $third
#fourth=$(qsub -W depend=afterok:$third 2b.visualize_align_with_igv.pbs)
#echo $fourth
#fifth=$(qsub -W depend=afterok:$fourth 2c.post_alignment_qc.pbs)
#echo $fifth
#sixth=$(qsub -W depend=afterok:$fifth 2d.globin_removal.pbs)
#echo $sixth
#seventh=$(qsub -W depend=afterok:$sixth 3a.gene_transcript_quantification.pbs)
#echo $seventh
