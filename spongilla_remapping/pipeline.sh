# Genome annotation downloaded from:    
# https://ftp.ebi.ac.uk/pub/ensemblorganisms/Spongilla_lacustris/GCA_949361645.1/ensembl/geneset/2024_03/

# COMMANDS TO RECOVER (DELETED ON 21/11/2025)

# read QC
multiqc --outdir 01_raw_reads/ 01_raw_reads/01_fastqc/

# prepare files for CellRanger
zcat 00_input/slac_genome/GCA.949361645.1/GCA_949361645.1_odSpoLacu1.1_genomic.fna.gz | sed -E 's/^>.+chromosome: />/; s/ .+$//' > 02a_cellranger_count/01_input/02_genome/GCA_949361645.1_odSpoLacu1.1_genomic_edited.fna

#########################
#     RUN PEAKS2UTR     #
#########################

cat 02a_cellranger_count/01_input/01_fastqs/*R2*gz > 02b_peaks2utr/R2.fastq.gz

STAR --runMode genomeGenerate --runThreadN 10 --genomeDir 02b_peaks2utr/01_genome_directory --genomeFastaFiles 02a_cellranger_count/02_slac_ref_genome/fasta/genome.fa --genomeSAindexNbases 12

STAR --genomeDir 02b_peaks2utr/01_genome_directory/ --outFilterMultimapNmax 10 --runThreadN 12 --readFilesIn 02b_peaks2utr/R2.fastq.gz --outFileNamePrefix 02b_peaks2utr/02_STAR_mapping/slac_R2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard

peaks2utr 02a_cellranger_count/02_slac_ref_genome/genes/genes.gtf slac_R2.Aligned.sortedByCoord.out.bam -o slac_R2_peaks2utr.gff -p 20


################################
#     RUN CELLRANGER AGAIN     #
################################

mkdir -p 02c_cellranger_count/01_input/02_genome
ln -s "$(realpath 02a_cellranger_count/01_input/01_fastqs)" 02c_cellranger_count/01_input/
ln -s "$(realpath 02a_cellranger_count/01_input/02_genome/GCA_949361645.1_odSpoLacu1.1_genomic_edited.fna)" 02c_cellranger_count/01_input/02_genome/
ln -s "$(realpath 02b_peaks2utr/slac_R2_peaks2utr_new.gff)" 02c_cellranger_count/01_input/02_genome/
ll 02c_cellranger_count/01_input/02_genome/

agat_convert_sp_gff2gtf.pl --gff 02b_peaks2utr/slac_R2_peaks2utr_new.gff -o 02b_peaks2utr/slac_R2_peaks2utr_new.gtf

cellranger mkgtf 02b_peaks2utr/slac_R2_peaks2utr_new.gtf 02c_cellranger_count/01_input/02_genome/slac_R2_peaks2utr_cellranger_format.gtf --attribute=gene_biotype:protein_coding
cellranger mkref --genome=02_slac_ref_genome --fasta=02c_cellranger_count/01_input/02_genome/GCA_949361645.1_odSpoLacu1.1_genomic_edited.fna --genes=02c_cellranger_count/01_input/02_genome/slac_R2_peaks2utr_cellranger_format.gtf

mv 02_slac_ref_genome/ 02c_cellranger_count/

sbatch --job-name=cellranger --cpus-per-task=20 --time=23:00:00 --mem=0 --export=NONE --account=evassvis --mail-type=BEGIN,END,FAIL --mail-user=fn76@le.ac.uk scripts/11_run_cellRanger_count_again.sh

mv SRR98410*count 02c_cellranger_count/

mkdir 02c_cellranger_count/03_aggregated/

while read j; do echo $j","$(realpath 02c_cellranger_count/"$j"_count/outs/molecule_info.h5); done <00_input/SRA_toDownload.ls | sed '1i sample_id,molecule_h5' > 02c_cellranger_count/03_aggregated/slac_3ext_aggr.csv

cellranger aggr --id=slac_3ext_aggr --csv=02c_cellranger_count/03_aggregated/slac_3ext_aggr.csv

mv slac_3ext_aggr/ 02c_cellranger_count/03_aggregated/


############################
#     PAIRWISE DIAMOND     #
############################

# to run pairwise blastp properly, you'd need to blast the two proteomes one against the other
# however, cellranger uses the gene ID (and not the fucking transcript ID) to assign UMIs
# and of course, the genome annotation nomenclature is not consistent
# (i.e., ENSLPGP00000015110.1 is not a transcript/protein derived from ENSLPGG00000015110.1, as it would be natural to do,
# but rather from ENSLPGG00000010709.1)
# To avoid mis-assignation, we keep both the transcript ID and the gene ID for the diamond search,
# then remove the transcript ID later on for SAMap, so that IDs in results match those computed by cellranger

# c95711-g1       ENSLPGP00000012527    ENSLPGG00000008927

sed -E 's/\..+gene:/_/; s/\..+$//' 00_input/slac_genome/pep.fa > 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genome_proteins.faa
sed -E 's/_i.+$//; s/_/-/' 00_input/slac_transcriptome/slac_pep.faa > 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_transcriptome_proteins.faa

diamond makedb --in 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genome_proteins.faa --db 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genome_proteins
diamond makedb --in 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_transcriptome_proteins.faa --db 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_transcriptome_proteins

diamond blastp --query 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genome_proteins.faa --db 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_transcriptome_proteins.dmnd --threads 8 --evalue 1e-6 --ultra-sensitive --outfmt 6 --out 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genomeVStranscriptome.tsv
diamond blastp --query 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_transcriptome_proteins.faa --db 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genome_proteins.dmnd --threads 8 --evalue 1e-6 --ultra-sensitive --outfmt 6 --out 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_transcriptomeVSgenome.tsv


##########################################
#     SAMap BETWEEN THE TWO DATASETS     #
##########################################

mkdir -p 04_SAMap/{01_pairwise_diamond/slacOriginalslacRemapped/,02_prepared_h5ads}

sed -E 's/^.+_//' 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_genomeVStranscriptome.tsv
sed -E 's/\t.+_/\t/' 03_slac_remapped_clustering/01_diamond_transcriptome_genome/slac_transcriptomeVSgenome.tsv > 04_SAMap/01_pairwise_diamond/slacOriginalslacRemapped/slacOriginal_to_slacRemapped.txt

python scripts/07_run_SAMap.py | tee -a 04_SAMap/slac_integrated_originalVSremapped.log
python scripts/08_get_SAMap_genePairs.py -p 04_SAMap/slac_integrated_originalVSremapped.pkl -o 04_SAMap/03_SAMap_statistics/ -t 0.2
python scripts/09_get_SAMap_mappingTables.py -p 04_SAMap/slac_integrated_originalVSremapped.pkl -o 04_SAMap/03_SAMap_statistics/ -n 100
