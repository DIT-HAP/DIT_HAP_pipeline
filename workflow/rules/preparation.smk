# rule all:
#     input:
#         "resources/pombase_data/2025-05-01/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.fai",
#         "resources/pombase_data/2025-05-01/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.amb"

# download pombase data from ftp server
# -----------------------------------------------------
checkpoint download_pombase_data:
    output:
        fasta = "resources/pombase_data/{release_version}/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa"
    params:
        release_version="{release_version}",
        download_dir="resources/pombase_data/{release_version}"
    log:
        "logs/preparation/download_pombase_data_{release_version}.log"
    message:
        "Downloading pombase data from ftp server"
    shell:
        "workflow/scripts/preparation/download_file_from_pombaseFTP.sh {params.release_version} {params.download_dir} &> {log}"

# index the genome fasta file with samtools faidx
# -----------------------------------------------------
rule samtools_faidx:
    input:
        rules.download_pombase_data.output.fasta
    output:
        f"{rules.download_pombase_data.output.fasta}.fai"
    log:
        "logs/preparation/samtools_faidx_{release_version}.log"
    message:
        "Indexing genome fasta file with samtools faidx"
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
        "logs/preparation/bwa_index_{release_version}.log"
    message:
        "Indexing genome fasta file with bwa"
    wrapper:
        "v6.2.0/bio/bwa/index"

# extract coding region from gff3 file
# -----------------------------------------------------
