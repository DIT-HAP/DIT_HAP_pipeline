#!/bin/bash

# This script is used to run the workflow for the analysis of the genome region classification

# set the input files
pombe_gff="/data/c/yangyusheng_optimized/resource/Pombase_FTP/release/20230701/Schizosaccharomyces_pombe_all_chromosomes.gff3"
pombe_fai="/data/c/yangyusheng_optimized/resource/Pombase_FTP/release/20230701/Schizosaccharomyces_pombe_all_chromosomes.fa.fai"
gene_product="/data/c/yangyusheng_optimized/resource/Pombase_FTP/release/20230701/gene_IDs_names_products.tsv"
gene_essentiality="../references/Hayles_2013_OB_merged_categories.xlsx"
release_version="20230701"

# set the output directory
output_dir="."/$release_version
mkdir -p $output_dir

# find the coding genes
python extract_coding_region_from_gff3.py \
    -g $pombe_gff \
    -p "/data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/20230701/PeptideStats.tsv" \
    -o $output_dir"/Coding_genes.bed"

# find the overlapped gene regions
python find_overlapped_gene_regions.py \
    -i $output_dir"/Coding_genes.bed" \
    -o $output_dir"/Overlapped_gene_regions.bed"

# annotate the intergenic regions
python annotate_intergenic_regions.py \
    -ic $output_dir"/Coding_genes.bed" \
    -iol $output_dir"/Overlapped_gene_regions.bed" \
    -f $pombe_fai \
    -o $output_dir"/Intergenic_regions.bed"

# merge these regions
python merge_genome_regions.py \
    -ic $output_dir"/Coding_genes.bed" \
    -iol $output_dir"/Overlapped_gene_regions.bed" \
    -iIGR $output_dir"/Intergenic_regions.bed" \
    -o1 $output_dir"/Genome_regions_overlapped_regions_not_separated.bed" \
    -o2 $output_dir"/Genome_regions_overlapped_regions_separated.bed"

# extract the CDS and intron regions
python extract_CDS_intron_from_coding_region_gff3.py \
    -g $pombe_gff \
    -o $output_dir"/Coding_gene_CDS_and_intron.bed"

# merge CDS, intron and intergenic regions
python merge_CDS_intron_intergenic.py \
    -ci $output_dir"/Coding_gene_CDS_and_intron.bed" \
    -i $output_dir"/Intergenic_regions.bed" \
    -o $output_dir"/Genome_regions_CDS_intron_IGR.bed" 

python annotate_the_genome_regions.py \
    -i $output_dir"/Genome_regions_CDS_intron_IGR.bed" \
    -n $gene_product \
    -e $gene_essentiality \
    -o $output_dir"/Genome_regions_CDS_intron_IGR_annotated.bed"

# python create_bed_for_each_nucleotide.py \
#     -g $pombe_gff \
#     -c 56 \
#     -o $output_dir"/nucletides_for_each_coding_region.bed"