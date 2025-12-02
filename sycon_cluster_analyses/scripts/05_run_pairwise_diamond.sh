#!/bin/bash

source ~/miniforge3/bin/activate diamond_env

for i in 01_proteomes/*ol.faa; do

	database_SP="$(basename ${i::-4} | sed -E 's/_.+$//')"

	diamond makedb --in $i --db ${i::-4} > /dev/null 2>&1

	for j in 01_proteomes/*ol.faa; do

		if [ $i != $j ]; then
			
			query_SP="$(basename ${j::-4} | sed -E 's/_.+$//')"
			
			if [[ ! -d 03_pairwise_diamond/"$query_SP""$database_SP" ]]; then
				
				OUTFILE="$(echo 03_pairwise_diamond/"$database_SP""$query_SP"/"$query_SP"_to_"$database_SP".txt)"
				mkdir 03_pairwise_diamond/"$database_SP""$query_SP"

			else
				OUTFILE="$(echo 03_pairwise_diamond/"$query_SP""$database_SP"/"$query_SP"_to_"$database_SP".txt)"
			fi

			echo "Running DIAMOND on $database_SP vs $query_SP"
			echo -e "\tCOMMAND: diamond blastp --query $j --db ${i::-4}.dmnd --threads 8 --evalue 1e-6 --ultra-sensitive --outfmt 6 --out $OUTFILE\n"


			diamond blastp --query $j --db ${i::-4}.dmnd --threads 8 --evalue 1e-6 --ultra-sensitive --outfmt 6 --out $OUTFILE > /dev/null 2>&1
		fi
	done
done
