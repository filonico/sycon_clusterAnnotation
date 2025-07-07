#!/bin/bash

source ~/miniforge3/bin/activate busco_env

for i in 01*/*ol.faa; do

	spID="$(basename $i | awk -F "_" '{print $1}')" &&
	LINEAGE="metazoa_odb10" &&
	MODE="" &&
		
	if [ "${i#*.}" == "faa" ]; then
		MODE="proteins"
	else
		MODE="transcriptome"
	fi

	busco -i $i -l $LINEAGE -o 02_busco/$(basename $i).busco -m "$MODE" -c 12
done

rm -rf busco_downloads/ busco*log
