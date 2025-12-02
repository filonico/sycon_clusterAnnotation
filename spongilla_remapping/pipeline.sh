# COMMANDS TO RECOVER (DELETED ON 21/11/2025)

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

cellranger mkgtf 02b_peaks2utr/slac_R2_peaks2utr_new.gtf 02c_cellranger_count/01_input/02_genome/slac_R2_peaks2utr_cellranger_format.gtf --attribute=gene_biotype:protein_codingcellranger mkref --genome=02_slac_ref_genome --fasta=02c_cellranger_count/01_input/02_genome/GCA_949361645.1_odSpoLacu1.1_genomic_edited.fna --genes=02c_cellranger_count/01_input/02_genome/slac_R2_peaks2utr_cellranger_format.gtf
cellranger mkref --genome=02_slac_ref_genome --fasta=02c_cellranger_count/01_input/02_genome/GCA_949361645.1_odSpoLacu1.1_genomic_edited.fna --genes=02c_cellranger_count/01_input/02_genome/slac_R2_peaks2utr_cellranger_format.gtf

mv 02_slac_ref_genome/ 02c_cellranger_count/

sbatch --job-name=cellranger --cpus-per-task=20 --time=23:00:00 --mem=0 --export=NONE --account=evassvis --mail-type=BEGIN,END,FAIL --mail-user=fn76@le.ac.uk scripts/11_run_cellRanger_count_again.sh

mv SRR98410*count 02c_cellranger_count/

mkdir 02c_cellranger_count/03_aggregated/

while read j; do echo $j","$(realpath 02c_cellranger_count/"$j"_count/outs/molecule_info.h5); done <00_input/SRA_toDownload.ls | sed '1i sample_id,molecule_h5' > 02c_cellranger_count/03_aggregated/slac_3ext_aggr.csv

cellranger aggr --id=slac_3ext_aggr --csv=02c_cellranger_count/03_aggregated/slac_3ext_aggr.csv

mv slac_3ext_aggr/ 02c_cellranger_count/03_aggregated/