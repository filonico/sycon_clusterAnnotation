#!/bin/bash

while read sample; do

    cellranger count --id="$sample"_count \
        --sample=$sample \
        --fastqs=02c_cellranger_count/01_input/01_fastqs \
        --transcriptome=02c_cellranger_count/02_slac_ref_genome \
        --nosecondary --create-bam=true
        
done <00_input/SRA_toDownload.ls
