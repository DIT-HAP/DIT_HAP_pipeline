#!/usr/bin/env bash

# ---
# Script to download various data files from PomBase FTP for a specific release.
#
# Usage:
#   ./download_file_from_pombaseFTP.sh <release_version> [download_dir]
#   ./download_file_from_pombaseFTP.sh -h|--help
#
# Parameters:
#   release_version  : PomBase release version (required, e.g., "2025-05-01")
#   download_dir     : Target directory for downloads (optional, defaults to "./release/<release_version>")
#
# Examples:
#   ./download_file_from_pombaseFTP.sh "2025-05-01"
#   ./download_file_from_pombaseFTP.sh "2025-05-01" "/path/to/custom/download/dir"
#
# The script will download files into the specified directory structure.
# It uses `wget -nc` to avoid re-downloading files if they already exist and
# `gunzip -f` to avoid overwriting uncompressed files if they exist.
# ---

set -e
set -o pipefail

# --- Help function ---
show_help() {
    cat << EOF
Script to download various data files from PomBase FTP for a specific release.

Usage:
    $0 <release_version> [download_dir]
    $0 -h|--help

Parameters:
    release_version  : PomBase release version (required, e.g., "2025-05-01")
    download_dir     : Target directory for downloads (optional, defaults to "./release/<release_version>")

Examples:
    $0 "2025-05-01"
    $0 "2025-05-01" "/path/to/custom/download/dir"

The script downloads various PomBase data files including:
- FASTA files (chromosomes, peptides)
- GFF annotation files
- Miscellaneous metadata (gene information, expression, protein features)
- Ontology files (GO, FYPO, Mondo)
- Export files (interactions, annotations)

Files are organized by type and use wget -nc (no-clobber) to avoid re-downloading.
EOF
}

# --- Parameter parsing ---
# Check for help option
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

# Check if release_version is provided
if [[ $# -lt 1 ]]; then
    echo "Error: release_version parameter is required."
    echo ""
    show_help
    exit 1
fi

# Set parameters
release_version="$1"
year=$(echo "${release_version}" | cut -d'-' -f1)
download_dir="${2:-./release/${release_version}}"

# Validate release_version is not empty
if [[ -z "$release_version" ]]; then
    echo "Error: release_version cannot be empty."
    exit 1
fi

# --- Main Script ---
echo "--- Starting PomBase data download for release: ${release_version} ---"
echo "Download directory: ${download_dir}"

# Create the target directory if it doesn't exist
mkdir -p "${download_dir}"

# --- Download FASTA files ---
echo "[INFO] Downloading FASTA files..."
# Create FASTA subfolder
mkdir -p "${download_dir}/genome_sequence_and_features"
# Chromosomes
wget -nc -P "${download_dir}/genome_sequence_and_features" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/genome_sequence_and_features/fasta_format/chromosomes/Schizosaccharomyces_pombe_all_chromosomes.fa"
# Feature sequences (peptides)
wget -nc -P "${download_dir}/genome_sequence_and_features" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/genome_sequence_and_features/fasta_format/feature_sequences/peptide.fa"

# --- Download yGFF file ---
echo "[INFO] Downloading GFF file..."
wget -nc -P "${download_dir}/genome_sequence_and_features" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/genome_sequence_and_features/gff_format/Schizosaccharomyces_pombe_all_chromosomes.gff3"

# --- Download metadata files ---
echo "[INFO] Downloading metadata files..."
## Basic gene information
mkdir -p "${download_dir}/Gene_metadata"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/gene_names_and_identifiers/gene_IDs_names_products.tsv"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/phenotypes_and_genotypes/gene_viability.tsv"

## Gene expression
mkdir -p "${download_dir}/RNA_metadata"
wget -nc -P "${download_dir}/RNA_metadata" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/gene_expression/qualitative_gene_expression.tsv"
wget -nc -P "${download_dir}/RNA_metadata" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/gene_expression/quantitative_gene_expression.tsv"

## Protein features
mkdir -p "${download_dir}/Protein_features"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/protein_features/peptide_stats.tsv"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/protein_features/protein_families_and_domains.tsv"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/protein_features/disordered_regions.tsv"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/protein_features/protein_modifications.tsv"

# --- Download Ontology files ---
echo "[INFO] Downloading ontology files..."
# Create ontologies subfolder
mkdir -p "${download_dir}/ontologies_and_associations"

# Download GO basic ontology
wget -nc -P "${download_dir}/ontologies_and_associations" "https://purl.obolibrary.org/obo/go/go-basic.obo"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://purl.obolibrary.org/obo/mondo.obo"
# Get latest FYPO version from GitHub API
echo "[INFO] Getting latest FYPO version from GitHub..."
latest_fypo_version=$(curl -s https://api.github.com/repos/pombase/fypo/releases/latest | grep '"tag_name"' | cut -d'"' -f4)
if [[ -z "$latest_fypo_version" ]]; then
    echo "Warning: Could not detect latest FYPO version, falling back to v2025-05-09"
    latest_fypo_version="v2025-05-09"
fi
echo "Latest FYPO version detected: ${latest_fypo_version}"
# Download FYPO files using detected version
wget -nc -P "${download_dir}/ontologies_and_associations" "https://github.com/pombase/fypo/releases/download/${latest_fypo_version}/fypo-base.obo"

## Slim ontologies
wget -nc -P "${download_dir}/ontologies_and_associations" "https://current.geneontology.org/ontology/subsets/goslim_pombe.obo" # Note: This downloads current GO slim

# Download slim ontology metadata
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/gene_ontology/bp_go_slim_terms.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/gene_ontology/cc_go_slim_terms.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/gene_ontology/mf_go_slim_terms.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/phenotypes_and_genotypes/fypo_slim_ids_and_names.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/human_disease_annotation/pombe_mondo_disease_slim_terms.tsv"

## Term associations
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/gene_ontology/gene_ontology_annotation.gaf.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/macromolecular_complexes/macromolecular_complex_annotation.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/phenotypes_and_genotypes/pombase_phenotype_annotation.phaf.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/phenotypes_and_genotypes/pombase_phenotype_annotation.eco.phaf.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/monthly_releases/${year}/pombase-${release_version}/human_disease_annotation/human_disease_association.tsv"

echo "--- PomBase data download finished successfully for release: ${release_version} ---"
echo "Files are organized in subfolders under: ${download_dir}"
echo "  - fasta/       : Genome and protein sequences"
echo "  - gff/         : Genome annotation files"  
echo "  - Gene_metadata/        : Basic gene information"
echo "  - RNA_metadata/        : Gene expression"
echo "  - Protein_features/        : Protein features"
echo "  - ontologies_and_associations/  : Ontology files (GO, FYPO v${latest_fypo_version}, Mondo) and associations"
echo "  - exports/     : Interaction and annotation export files"






