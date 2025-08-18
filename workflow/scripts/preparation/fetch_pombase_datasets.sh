#!/usr/bin/env bash

#===============================================================================
# PomBase Data Download Script
#===============================================================================
#
# DESCRIPTION:
#   Downloads comprehensive datasets from PomBase (Schizosaccharomyces pombe 
#   database) for a specified release version. This script organizes data into 
#   structured directories and includes genome sequences, annotations, metadata,
#   and ontology files.
#
# USAGE:
#   ./download_file_from_pombaseFTP.sh <release_version> [download_directory]
#
# ARGUMENTS:
#   release_version    - PomBase release version (e.g., "2024-06-05")
#   download_directory - Target directory (optional, default: ./release/<version>)
#
# EXAMPLES:
#   ./download_file_from_pombaseFTP.sh "2024-06-05"
#   ./download_file_from_pombaseFTP.sh "2024-06-05" "/data/pombase"
#
# OUTPUT STRUCTURE:
#   <download_dir>/
#   ├── genome_sequence_and_features/    # FASTA and GFF files
#   ├── Gene_metadata/                   # Gene information and viability
#   ├── RNA_metadata/                    # Expression data
#   ├── Protein_features/                # Protein annotations
#   └── ontologies_and_associations/     # Ontologies and term associations
#
# DEPENDENCIES:
#   - wget (for file downloads)
#   - curl (for API calls)
#   - Internet connection
#
# AUTHOR: Bioinformatics Pipeline
# VERSION: 2.0
# UPDATED: $(date +%Y-%m-%d)
#===============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

#===============================================================================
# CONFIGURATION AND CONSTANTS
#===============================================================================

readonly SCRIPT_NAME="$(basename "${0}")"
readonly DEFAULT_FALLBACK_FYPO_VERSION="v2025-08-13"

# Base URLs for different data sources
readonly POMBASE_BASE_URL="https://www.pombase.org/monthly_releases"
readonly FYPO_GITHUB_API="https://api.github.com/repos/pombase/fypo/releases/latest"
readonly FYPO_RELEASE_URL="https://github.com/pombase/fypo/releases/download"
readonly GO_OBO_URL="https://purl.obolibrary.org/obo/go/go-basic.obo"
readonly MONDO_OBO_URL="https://purl.obolibrary.org/obo/mondo.obo"
readonly GO_SLIM_URL="https://current.geneontology.org/ontology/subsets/goslim_pombe.obo"

# Wget common options
readonly WGET_OPTS="-nc --progress=bar:force --timeout=30 --tries=3"

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

# Print formatted log messages with timestamps
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*" >&2
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*" >&2
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

# Display usage information
show_usage() {
    cat << EOF
Usage: ${SCRIPT_NAME} <release_version> [download_directory]

Download PomBase datasets for a specific release version.

Arguments:
    release_version     PomBase release version (format: YYYY-MM-DD)
    download_directory  Target directory (optional, default: ./release/<version>)

Examples:
    ${SCRIPT_NAME} "2024-06-05"
    ${SCRIPT_NAME} "2024-06-05" "/data/pombase"

For more information, see the script header documentation.
EOF
}

# Validate release version format
validate_release_version() {
    local version="$1"
    
    if [[ ! "${version}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
        log_error "Invalid release version format: '${version}'"
        log_error "Expected format: YYYY-MM-DD (e.g., 2025-08-01)"
        return 1
    fi
    
    return 0
}

# Create directory with error handling
create_directory() {
    local dir_path="$1"
    local description="$2"
    
    if ! mkdir -p "${dir_path}"; then
        log_error "Failed to create ${description} directory: ${dir_path}"
        return 1
    fi
    
    log_info "Created directory: ${dir_path}"
    return 0
}

# Download file with enhanced error handling and validation
download_file() {
    local url="$1"
    local target_dir="$2"
    local description="${3:-file}"
    
    log_info "Downloading ${description}: $(basename "${url}")"
    
    if ! wget ${WGET_OPTS} -P "${target_dir}" "${url}"; then
        log_error "Failed to download ${description}: ${url}"
        return 1
    fi
    
    # Verify file was downloaded and has non-zero size
    local filename
    filename="$(basename "${url}")"
    local filepath="${target_dir}/${filename}"
    
    if [[ ! -f "${filepath}" ]]; then
        log_error "Downloaded file not found: ${filepath}"
        return 1
    fi
    
    if [[ ! -s "${filepath}" ]]; then
        log_error "Downloaded file is empty: ${filepath}"
        return 1
    fi
    
    log_info "Successfully downloaded: ${filename} ($(du -h "${filepath}" | cut -f1))"
    return 0
}

# Get latest FYPO version from GitHub API
get_latest_fypo_version() {
    log_info "Fetching latest FYPO version from GitHub API..."
    
    local fypo_version
    fypo_version=$(curl -s --max-time 10 "${FYPO_GITHUB_API}" | grep '"tag_name"' | cut -d'"' -f4 2>/dev/null)
    
    if [[ -z "${fypo_version}" ]]; then
        log_warn "Could not detect latest FYPO version from API"
        log_warn "Falling back to default version: ${DEFAULT_FALLBACK_FYPO_VERSION}"
        fypo_version="${DEFAULT_FALLBACK_FYPO_VERSION}"
    else
        log_info "Latest FYPO version detected: ${fypo_version}"
    fi
    
    echo "${fypo_version}"
}

#===============================================================================
# DOWNLOAD FUNCTIONS
#===============================================================================

# Download genome sequence and annotation files
download_genome_files() {
    local base_url="$1"
    local target_dir="$2"
    
    log_info "Starting genome sequence and annotation downloads..."
    
    local genome_dir="${target_dir}/genome_sequence_and_features"
    create_directory "${genome_dir}" "genome sequence and features" || return 1
    
    # Chromosome sequences
    download_file \
        "${base_url}/genome_sequence_and_features/fasta_format/chromosomes/Schizosaccharomyces_pombe_all_chromosomes.fa" \
        "${genome_dir}" \
        "chromosome sequences" || return 1
    
    # Peptide sequences
    download_file \
        "${base_url}/genome_sequence_and_features/fasta_format/feature_sequences/peptide.fa" \
        "${genome_dir}" \
        "peptide sequences" || return 1
    
    # GFF annotation file
    download_file \
        "${base_url}/genome_sequence_and_features/gff_format/Schizosaccharomyces_pombe_all_chromosomes.gff3" \
        "${genome_dir}" \
        "genome annotation (GFF3)" || return 1
    
    log_info "Genome files download completed"
    return 0
}

# Download gene metadata files
download_gene_metadata() {
    local base_url="$1"
    local target_dir="$2"
    
    log_info "Starting gene metadata downloads..."
    
    local metadata_dir="${target_dir}/Gene_metadata"
    create_directory "${metadata_dir}" "gene metadata" || return 1
    
    # Gene IDs, names, and products
    download_file \
        "${base_url}/gene_names_and_identifiers/gene_IDs_names_products.tsv" \
        "${metadata_dir}" \
        "gene IDs and names" || return 1
    
    # Gene viability data
    download_file \
        "${base_url}/phenotypes_and_genotypes/gene_viability.tsv" \
        "${metadata_dir}" \
        "gene viability data" || return 1
    
    log_info "Gene metadata download completed"
    return 0
}

# Download RNA/expression data
download_rna_metadata() {
    local base_url="$1"
    local target_dir="$2"
    
    log_info "Starting RNA metadata downloads..."
    
    local rna_dir="${target_dir}/RNA_metadata"
    create_directory "${rna_dir}" "RNA metadata" || return 1
    
    # Qualitative gene expression
    download_file \
        "${base_url}/gene_expression/qualitative_gene_expression.tsv" \
        "${rna_dir}" \
        "qualitative gene expression" || return 1
    
    # Quantitative gene expression
    download_file \
        "${base_url}/gene_expression/quantitative_gene_expression.tsv" \
        "${rna_dir}" \
        "quantitative gene expression" || return 1
    
    log_info "RNA metadata download completed"
    return 0
}

# Download protein feature files
download_protein_features() {
    local base_url="$1"
    local target_dir="$2"
    
    log_info "Starting protein features downloads..."
    
    local protein_dir="${target_dir}/Protein_features"
    create_directory "${protein_dir}" "protein features" || return 1
    
    # Peptide statistics
    download_file \
        "${base_url}/protein_features/peptide_stats.tsv" \
        "${protein_dir}" \
        "peptide statistics" || return 1
    
    # Protein families and domains
    download_file \
        "${base_url}/protein_features/protein_families_and_domains.tsv" \
        "${protein_dir}" \
        "protein families and domains" || return 1
    
    # Disordered regions
    download_file \
        "${base_url}/protein_features/disordered_regions.tsv" \
        "${protein_dir}" \
        "disordered regions" || return 1
    
    # Protein modifications
    download_file \
        "${base_url}/protein_features/protein_modifications.tsv" \
        "${protein_dir}" \
        "protein modifications" || return 1
    
    log_info "Protein features download completed"
    return 0
}

# Download ontology files and associations
download_ontologies() {
    local base_url="$1"
    local target_dir="$2"
    local fypo_version="$3"
    
    log_info "Starting ontology and association downloads..."
    
    local onto_dir="${target_dir}/ontologies_and_associations"
    create_directory "${onto_dir}" "ontologies and associations" || return 1
    
    # Core ontology files
    download_file "${GO_OBO_URL}" "${onto_dir}" "GO basic ontology" || return 1
    download_file "${MONDO_OBO_URL}" "${onto_dir}" "Mondo ontology" || return 1
    download_file "${FYPO_RELEASE_URL}/${fypo_version}/fypo-base.obo" "${onto_dir}" "FYPO ontology" || return 1
    download_file "${GO_SLIM_URL}" "${onto_dir}" "GO slim pombe" || return 1
    
    # Slim ontology metadata
    local slim_files=(
        "gene_ontology/bp_go_slim_terms.tsv"
        "gene_ontology/cc_go_slim_terms.tsv"
        "gene_ontology/mf_go_slim_terms.tsv"
        "phenotypes_and_genotypes/fypo_slim_ids_and_names.tsv"
        "human_disease_annotation/pombe_mondo_disease_slim_terms.tsv"
    )
    
    for file in "${slim_files[@]}"; do
        download_file \
            "${base_url}/${file}" \
            "${onto_dir}" \
            "$(basename "${file}" .tsv) metadata" || return 1
    done
    
    # Term associations
    local association_files=(
        "gene_ontology/gene_ontology_annotation.gaf.tsv"
        "macromolecular_complexes/macromolecular_complex_annotation.tsv"
        "phenotypes_and_genotypes/pombase_phenotype_annotation.phaf.tsv"
        "phenotypes_and_genotypes/pombase_phenotype_annotation.eco.phaf.tsv"
        "human_disease_annotation/human_disease_association.tsv"
    )
    
    for file in "${association_files[@]}"; do
        download_file \
            "${base_url}/${file}" \
            "${onto_dir}" \
            "$(basename "${file}" .tsv) associations" || return 1
    done
    
    log_info "Ontology and association downloads completed"
    return 0
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

main() {
    # Parse and validate arguments
    local release_version="$1"
    local download_dir="${2:-./release/${release_version}}"
    
    log_info "Starting PomBase data download for release: ${release_version}"
    log_info "Target directory: ${download_dir}"
    
    # Validate inputs
    validate_release_version "${release_version}" || {
        show_usage
        exit 1
    }
    
    # Extract year from release version for URL construction
    local year
    year=$(echo "${release_version}" | cut -d'-' -f1)
    local base_url="${POMBASE_BASE_URL}/${year}/pombase-${release_version}"
    
    # Create main download directory
    create_directory "${download_dir}" "main download" || exit 1
    
    # Get latest FYPO version
    local fypo_version
    fypo_version=$(get_latest_fypo_version)
    
    # Execute download phases
    log_info "Beginning multi-phase download process..."
    
    download_genome_files "${base_url}" "${download_dir}" || exit 1
    download_gene_metadata "${base_url}" "${download_dir}" || exit 1
    download_rna_metadata "${base_url}" "${download_dir}" || exit 1
    download_protein_features "${base_url}" "${download_dir}" || exit 1
    download_ontologies "${base_url}" "${download_dir}" "${fypo_version}" || exit 1
    
    # Generate completion summary
    log_info "All downloads completed successfully!"
    echo ""
    echo "==================== DOWNLOAD SUMMARY ===================="
    echo "Release Version: ${release_version}"
    echo "Download Location: ${download_dir}"
    echo "FYPO Version: ${fypo_version}"
    echo ""
    echo "Directory Structure:"
    echo "  ${download_dir}/"
    echo "  ├── genome_sequence_and_features/    # FASTA and GFF files"
    echo "  ├── Gene_metadata/                   # Gene information and viability"
    echo "  ├── RNA_metadata/                    # Expression data (qualitative/quantitative)"
    echo "  ├── Protein_features/                # Protein annotations and modifications"
    echo "  └── ontologies_and_associations/     # Ontologies (GO, FYPO, Mondo) and associations"
    echo ""
    echo "Total downloaded files: $(find "${download_dir}" -type f | wc -l)"
    echo "Total download size: $(du -sh "${download_dir}" | cut -f1)"
    echo "=========================================================="
}

# Script entry point with argument validation
if [[ $# -lt 1 || $# -gt 2 ]]; then
    log_error "Invalid number of arguments"
    show_usage
    exit 1
fi

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_usage
    exit 0
fi

# Validate required tools
for tool in wget curl; do
    if ! command -v "${tool}" &> /dev/null; then
        log_error "Required tool '${tool}' is not installed or not in PATH"
        exit 1
    fi
done

# Execute main function
main "$@"
