# download pombase data from ftp server
# -----------------------------------------------------
rule download_pombase_data:
    output:
        fasta = "resources/pombase_data/{release_version}/genome_sequence_and_features/Schizosaccharomyces_pombe_all_chromosomes.fa",
        gff = "resources/pombase_data/{release_version}/genome_sequence_and_features/Schizosaccharomyces_pombe_all_chromosomes.gff3",
        peptide_stats = "resources/pombase_data/{release_version}/Protein_features/peptide_stats.tsv",
        gene_IDs_names_products = "resources/pombase_data/{release_version}/Gene_metadata/gene_IDs_names_products.tsv",
        gene_viability = "resources/pombase_data/{release_version}/Gene_metadata/gene_viability.tsv"
    params:
        release_version="{release_version}",
        download_dir="resources/pombase_data/{release_version}"
    log:
        f"logs/{project_name}/preparation/download_pombase_data_{{release_version}}.log"
    message:
        "*** Downloading pombase data from ftp server"
    shell:
        "workflow/scripts/preparation/fetch_pombase_datasets.sh {params.release_version} {params.download_dir} &> {log}"

# index the genome fasta file with samtools faidx
# -----------------------------------------------------
rule samtools_faidx:
    input:
        rules.download_pombase_data.output.fasta
    output:
        "resources/pombase_data/{release_version}/genome_sequence_and_features/Schizosaccharomyces_pombe_all_chromosomes.fa.fai"
    log:
        f"logs/{project_name}/preparation/samtools_faidx_{{release_version}}.log"
    message:
        "*** Indexing genome fasta file with samtools faidx"
    wrapper:
        "v6.2.0/bio/samtools/faidx"

# index the genome fasta file with bwa
# -----------------------------------------------------
rule bwa_index:
    input:
        rules.download_pombase_data.output.fasta
    output:
        multiext(rules.download_pombase_data.output.fasta, ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        f"logs/{project_name}/preparation/bwa_index_{{release_version}}.log"
    message:
        "*** Indexing genome fasta file with bwa"
    wrapper:
        "v6.2.0/bio/bwa/index"

# extract genome region from gff3 file
# -----------------------------------------------------
rule extract_genome_region:
    input:
        rules.download_pombase_data.output.gff,
        rules.samtools_faidx.output,
        rules.download_pombase_data.output.peptide_stats,
        rules.download_pombase_data.output.gene_IDs_names_products,
        rules.download_pombase_data.output.gene_viability,
        rules.bwa_index.output
    output:
        primary_transcripts_bed = "resources/pombase_data/{release_version}/genome_region/coding_gene_primary_transcripts.bed",
        intergenic_regions_bed = "resources/pombase_data/{release_version}/genome_region/intergenic_regions.bed",
        non_coding_rna_bed = "resources/pombase_data/{release_version}/genome_region/non_coding_rna.bed",
        genome_intervals_bed = "resources/pombase_data/{release_version}/genome_region/genome_intervals.bed",
        overlapped_region_bed = "resources/pombase_data/{release_version}/genome_region/overlapped_region.bed"
    log:
        notebook=f"logs/{project_name}/preparation/gff_processing_and_annotation_{{release_version}}.ipynb"
    params:
        hayles_viability_path = "resources/Hayles_2013_OB_merged_categories.xlsx"
    message:
        "*** Extracting genome region from gff3 file"
    conda:
        "../envs/pybedtools.yml"
    notebook:
        "../notebooks/gff_processing_and_annotation.ipynb"