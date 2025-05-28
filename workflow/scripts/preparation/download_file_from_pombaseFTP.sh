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
# `gunzip -n` to avoid overwriting uncompressed files if they exist.
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
mkdir -p "${download_dir}/fasta"

# Chromosomes
wget -nc -P "${download_dir}/fasta" "https://www.pombase.org/data/releases/pombase-${release_version}/fasta/chromosomes/Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
gunzip -n "${download_dir}/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa.gz"

# Feature sequences (peptides)
wget -nc -P "${download_dir}/fasta" "https://www.pombase.org/data/releases/pombase-${release_version}/fasta/feature_sequences/peptide.fa.gz"
gunzip -n "${download_dir}/fasta/peptide.fa.gz"

# --- Download GFF file ---
echo "[INFO] Downloading GFF file..."
# Create GFF subfolder
mkdir -p "${download_dir}/gff"
wget -nc -P "${download_dir}/gff" "https://www.pombase.org/data/releases/pombase-${release_version}/gff/Schizosaccharomyces_pombe_all_chromosomes.gff3"

# --- Download Miscellaneous metadata files ---
echo "[INFO] Downloading miscellaneous metadata files..."

## Basic gene information
mkdir -p "${download_dir}/Gene_metadata"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/gene_IDs_names.tsv"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/gene_IDs_names_products.tsv"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/sysID2product.tsv"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/sysID2product.rna.tsv"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/uniprot_id_mapping.tsv"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/pseudogeneIDs.tsv"
wget -nc -P "${download_dir}/Gene_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/FYPOviability.tsv"

## Gene expression
mkdir -p "${download_dir}/RNA_metadata"
wget -nc -P "${download_dir}/RNA_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/full_gene_expression_table.tsv"
wget -nc -P "${download_dir}/RNA_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/qualitative_gene_expression.tsv"
wget -nc -P "${download_dir}/RNA_metadata" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/htp_gene_expression_table.tsv"

## Protein features
mkdir -p "${download_dir}/Protein_features"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/PeptideStats.tsv"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/ProteinFeatures.tsv"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/disordered_regions.tsv"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/interactions.tsv"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/modifications.tsv"
wget -nc -P "${download_dir}/Protein_features" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/transmembrane_domain_coords_and_seqs.tsv"

# --- Download Ontology files ---
echo "[INFO] Downloading ontology files..."
# Create ontologies subfolder
mkdir -p "${download_dir}/ontologies_and_associations"

# Download PomBase ontology files (tied to release version)
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/ontologies/go-basic.obo"

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
wget -nc -P "${download_dir}/ontologies_and_associations" "https://github.com/pombase/fypo/releases/download/${latest_fypo_version}/fypo-full.obo"

# Download other ontology files (not version-specific)
wget -nc -P "${download_dir}/ontologies_and_associations" "http://purl.obolibrary.org/obo/mondo.obo"

## Slim ontologies
wget -nc -P "${download_dir}/ontologies_and_associations" "https://current.geneontology.org/ontology/subsets/goslim_pombe.obo" # Note: This downloads current GO slim
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/ontologies/fypo-simple.obo"

# Download slim ontology metadata to misc folder
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/bp_goslim_pombe_ids_and_names.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/cc_goslim_pombe_ids_and_names.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/mf_goslim_pombe_ids_and_names.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/fypo_slim_ids_and_names.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/pombe_mondo_slim_ids_and_names.tsv"

## Term associations
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/go_style_gaf.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/pombase_style_gaf.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/extended_pombase_style_gaf.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/Complex_annotation.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/disease_association.tsv"
wget -nc -P "${download_dir}/ontologies_and_associations" "https://www.pombase.org/data/releases/pombase-${release_version}/misc/single_locus_phenotype_annotations_taxon_4896.phaf"

# --- Download Export files ---
echo "[INFO] Downloading export files..."
# Create exports subfolder
mkdir -p "${download_dir}/exports"

wget -nc -P "${download_dir}/exports" "https://www.pombase.org/data/releases/pombase-${release_version}/exports/pombase-go-physical-interactions.tsv.gz"
gunzip -n "${download_dir}/exports/pombase-go-physical-interactions.tsv.gz"

wget -nc -P "${download_dir}/exports" "https://www.pombase.org/data/releases/pombase-${release_version}/exports/pombase-interactions-since-v62-2017-01-30.gz"
gunzip -n "${download_dir}/exports/pombase-interactions-since-v62-2017-01-30.gz" # Note: Filename indicates 'since-v62', might not align with release_version

wget -nc -P "${download_dir}/exports" "https://www.pombase.org/data/releases/pombase-${release_version}/pombase-${release_version}.gaf.gz"
gunzip -n "${download_dir}/exports/pombase-${release_version}.gaf.gz"

wget -nc -P "${download_dir}/exports" "https://www.pombase.org/data/releases/pombase-${release_version}/pombase-${release_version}.phaf.gz"
gunzip -n "${download_dir}/exports/pombase-${release_version}.phaf.gz"

echo "--- PomBase data download finished successfully for release: ${release_version} ---"
echo "Files are organized in subfolders under: ${download_dir}"
echo "  - fasta/       : Genome and protein sequences"
echo "  - gff/         : Genome annotation files"  
echo "  - Gene_metadata/        : Basic gene information"
echo "  - RNA_metadata/        : Gene expression"
echo "  - Protein_features/        : Protein features"
echo "  - ontologies_and_associations/  : Ontology files (GO, FYPO v${latest_fypo_version}, Mondo) and associations"
echo "  - exports/     : Interaction and annotation export files"






