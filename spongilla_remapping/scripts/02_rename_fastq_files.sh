#!/bin/bash

for i in 01_raw_reads/SRR*/*fastq.gz; do

    FULLPATH=$(realpath $i)
    SAMPLE=$(basename $i | cut -d_ -f1)
    SAMPLE_NUM=$(grep $SAMPLE 00_input/SRA_toDownload.tsv | cut -f3)
    LANE_NUM=$(grep $SAMPLE 00_input/SRA_toDownload.tsv | cut -f4)
    READ=$(basename $i | cut -d_ -f2 | sed -E 's/\..+$//')

    # files with "_1" specification are the index files, which are not required for CellRanger
    if [ $READ == "1" ]; then
        continue
    elif [ $READ == "2" ]; then
        READ_TYPE="R1"
    else
        READ_TYPE="R2"
    fi

    FILENAME=$(echo $SAMPLE"_"$SAMPLE_NUM"_"$LANE_NUM"_"$READ_TYPE"_001.fastq.gz")

    ln -s $FULLPATH 02a_cellranger_count/01_input/01_fastqs/$FILENAME

    #cp $FULLPATH 02a_cellranger_count/01_input/01_fastqs/$FILENAME

done