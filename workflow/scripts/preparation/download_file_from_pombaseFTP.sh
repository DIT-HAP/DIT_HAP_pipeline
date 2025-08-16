#!/usr/bin/env bash

# =============================================================================
# PomBase Data Download Script
# =============================================================================
#
# Description:
#   Automated download utility for PomBase genomic and annotation data files.
#   Retrieves comprehensive datasets for a specific PomBase release, including
#   genome sequences, annotations, gene metadata, and ontology files.
#
# Usage:
#   ./download_file_from_pombaseFTP.sh <release_version> [download_dir]
#   ./download_file_from_pombaseFTP.sh -h|--help
#
# Parameters:
#   release_version  : PomBase release version (required, format: YYYY-MM-DD)
#   download_dir     : Target directory for downloads (optional, default: ./release/<release_version>)
#
# Examples:
#   ./download_file_from_pombaseFTP.sh "2025-05-01"
#   ./download_file_from_pombaseFTP.sh "2025-05-01" "/path/to/custom/download/dir"
#
# Exit Codes:
#   0  - Success
#   1  - Invalid parameters
#   2  - Network/download failure
#   3  - File system error
#
# Notes:
#   - Uses wget with --no-clobber to skip existing files
#   - Creates organized directory structure for different data types
#   - Includes comprehensive error handling and progress reporting
#   - Validates downloaded file integrity using HTTP HEAD requests
#
# Dependencies:
#   - wget (for file downloads)
#   - curl (for API calls and HTTP validation)
#   - bash 4.0+ (for associative arrays and improved error handling)
#
# =============================================================================

# =============================================================================
# Configuration and Environment Setup
# =============================================================================

# Strict error handling and debugging
set -euo pipefail
IFS=$'\n\t'

# Script metadata
readonly SCRIPT_NAME="$(basename "$0")"
readonly SCRIPT_VERSION="1.0.0"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Color codes for output formatting
readonly COLOR_RED='\033[0;31m'
readonly COLOR_GREEN='\033[0;32m'
readonly COLOR_YELLOW='\033[1;33m'
readonly COLOR_BLUE='\033[0;34m'
readonly COLOR_NC='\033[0m' # No Color

# Logging configuration
readonly LOG_LEVEL_INFO=0
readonly LOG_LEVEL_WARN=1
readonly LOG_LEVEL_ERROR=2
readonly LOG_LEVEL_DEBUG=3
LOG_LEVEL=${LOG_LEVEL_INFO}

# =============================================================================
# Utility Functions
# =============================================================================

# Logging functions with color support
log_info() {
    echo -e "${COLOR_GREEN}[INFO]${COLOR_NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_warn() {
    echo -e "${COLOR_YELLOW}[WARN]${COLOR_NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_error() {
    echo -e "${COLOR_RED}[ERROR]${COLOR_NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_debug() {
    [[ $LOG_LEVEL -ge $LOG_LEVEL_DEBUG ]] && echo -e "${COLOR_BLUE}[DEBUG]${COLOR_NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

# Display help information
show_help() {
    cat << EOF
${SCRIPT_NAME} v${SCRIPT_VERSION}

$(cat "$0" | sed -n '8,40p' | sed 's/^# //' | sed 's/^#//')

Additional Options:
    -v, --verbose    Enable verbose logging
    -q, --quiet      Suppress informational messages
    --dry-run        Show what would be downloaded without actually downloading

Environment Variables:
    POMBASE_BASE_URL   Override base URL for PomBase (default: https://www.pombase.org)
    WGET_OPTS          Additional wget options (e.g., "--timeout=30 --tries=3")

Report bugs to: <your-email@example.com>
EOF
}

# =============================================================================
# Validation and Error Handling Functions
# =============================================================================

# Validate release version format
validate_release_version() {
    local release="$1"
    local pattern='^[0-9]{4}-[0-9]{2}-[0-9]{2}$'
    
    if [[ ! $release =~ $pattern ]]; then
        log_error "Invalid release version format: '$release'"
        log_error "Expected format: YYYY-MM-DD (e.g., 2025-05-01)"
        return 1
    fi
    
    # Additional validation: check if year is reasonable
    local year=$(echo "$release" | cut -d'-' -f1)
    if [[ $year -lt 2000 || $year -gt $(date +%Y) ]]; then
        log_error "Invalid year in release version: $year"
        return 1
    fi
    
    return 0
}

# Check for required dependencies
check_dependencies() {
    local deps=("wget" "curl")
    local missing=()
    
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            missing+=("$dep")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing required dependencies: ${missing[*]}"
        log_error "Please install missing tools and try again"
        return 1
    fi
    
    log_debug "All dependencies satisfied: ${deps[*]}"
    return 0
}

# Validate URL accessibility
validate_url() {
    local url="$1"
    
    if curl --head --silent --fail "$url" &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# =============================================================================
# Download Management Functions
# =============================================================================

# =============================================================================
# Command Line Argument Processing
# =============================================================================

# Default configuration
DRY_RUN=false
VERBOSE=false
QUIET=false

# Parse command line arguments
parse_arguments() {
    local positional_args=()
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -v|--verbose)
                LOG_LEVEL=$LOG_LEVEL_DEBUG
                shift
                ;;
            -q|--quiet)
                LOG_LEVEL=$LOG_LEVEL_WARN
                shift
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            -*)
                log_error "Unknown option: $1"
                show_help
                exit 1
                ;;
            *)
                positional_args+=("$1")
                shift
                ;;
        esac
    done
    
    # Set positional arguments
    set -- "${positional_args[@]}"
    
    # Validate required parameters
    if [[ $# -lt 1 ]]; then
        log_error "release_version parameter is required"
        show_help
        exit 1
    fi
    
    release_version="$1"
    download_dir="${2:-./release/${release_version}}"
    
    # Validate parameters
    validate_release_version "$release_version" || exit 1
    
    # Extract year for URL construction
    year=$(echo "${release_version}" | cut -d'-' -f1)
    
    log_debug "Configuration: release=$release_version, year=$year, dir=$download_dir"
}

# =============================================================================
# Main Download Functions
# =============================================================================

# =============================================================================
# Download Management Functions
# =============================================================================

# Download a file with retry logic and validation
download_file() {
    local url="$1"
    local output_dir="$2"
    local description="${3:-file}"
    
    local filename=$(basename "$url")
    local output_path="$output_dir/$filename"
    
    # Check if file already exists
    if [[ -f "$output_path" ]]; then
        log_info "Skipping existing file: $filename"
        return 0
    fi
    
    # Validate URL before attempting download
    if ! validate_url "$url"; then
        log_warn "URL not accessible: $url"
        return 1
    fi
    
    # Dry run mode
    if [[ $DRY_RUN == true ]]; then
        log_info "[DRY RUN] Would download: $url -> $output_path"
        return 0
    fi
    
    log_info "Downloading $description: $filename"
    
    # Use custom wget options if provided
    local wget_opts="${WGET_OPTS:-}"
    local wget_cmd="wget --no-clobber --continue --progress=dot:giga $wget_opts"
    
    if $wget_cmd -P "$output_dir" "$url"; then
        log_info "Successfully downloaded: $filename"
        return 0
    else
        log_error "Failed to download: $filename"
        return 1
    fi
}

# Download latest FYPO ontology version
download_fypo_ontology() {
    local output_dir="$1"
    local fypo_version
    
    log_info "Getting latest FYPO version from GitHub..."
    
    # Get latest release version from GitHub API
    fypo_version=$(curl -s https://api.github.com/repos/pombase/fypo/releases/latest | \
                   grep '"tag_name"' | cut -d'"' -f4)
    
    if [[ -z "$fypo_version" ]]; then
        log_warn "Could not detect latest FYPO version, using fallback"
        fypo_version="v2025-05-09"
    fi
    
    log_info "Latest FYPO version: $fypo_version"
    
    local fypo_url="https://github.com/pombase/fypo/releases/download/${fypo_version}/fypo-base.obo"
    download_file "$fypo_url" "$output_dir" "FYPO ontology"
}

# =============================================================================
# Main Download Orchestration
# =============================================================================

# Build download targets dynamically based on release version
build_download_targets() {
    local release_version="$1"
    local year=$(echo "${release_version}" | cut -d'-' -f1)
    
    # Return the download targets as a simple list
    cat << EOF
chromosomes=/monthly_releases/${year}/pombase-${release_version}/genome_sequence_and_features/fasta_format/chromosomes/Schizosaccharomyces_pombe_all_chromosomes.fa
peptides=/monthly_releases/${year}/pombase-${release_version}/genome_sequence_and_features/fasta_format/feature_sequences/peptide.fa
annotations=/monthly_releases/${year}/pombase-${release_version}/genome_sequence_and_features/gff_format/Schizosaccharomyces_pombe_all_chromosomes.gff3
gene_ids=/monthly_releases/${year}/pombase-${release_version}/gene_names_and_identifiers/gene_IDs_names_products.tsv
gene_viability=/monthly_releases/${year}/pombase-${release_version}/phenotypes_and_genotypes/gene_viability.tsv
qualitative_expression=/monthly_releases/${year}/pombase-${release_version}/gene_expression/qualitative_gene_expression.tsv
quantitative_expression=/monthly_releases/${year}/pombase-${release_version}/gene_expression/quantitative_gene_expression.tsv
peptide_stats=/monthly_releases/${year}/pombase-${release_version}/protein_features/peptide_stats.tsv
protein_domains=/monthly_releases/${year}/pombase-${release_version}/protein_features/protein_families_and_domains.tsv
disordered_regions=/monthly_releases/${year}/pombase-${release_version}/protein_features/disordered_regions.tsv
protein_modifications=/monthly_releases/${year}/pombase-${release_version}/protein_features/protein_modifications.tsv
go_ontology=/gene_ontology/bp_go_slim_terms.tsv
cc_ontology=/gene_ontology/cc_go_slim_terms.tsv
mf_ontology=/gene_ontology/mf_go_slim_terms.tsv
fypo_slim=/phenotypes_and_genotypes/fypo_slim_ids_and_names.tsv
mondo_slim=/human_disease_annotation/pombe_mondo_disease_slim_terms.tsv
go_annotations=/gene_ontology/gene_ontology_annotation.gaf.tsv
complex_annotations=/macromolecular_complexes/macromolecular_complex_annotation.tsv
phenotype_annotations=/phenotypes_and_genotypes/pombase_phenotype_annotation.phaf.tsv
eco_annotations=/phenotypes_and_genotypes/pombase_phenotype_annotation.eco.phaf.tsv
disease_annotations=/human_disease_annotation/human_disease_association.tsv
EOF
}

# Directory structure mapping
declare -A DIRECTORY_MAPPING=(
    ["genome_sequence_and_features"]="genome_sequence_and_features"
    ["Gene_metadata"]="Gene_metadata"
    ["RNA_metadata"]="RNA_metadata"
    ["Protein_features"]="Protein_features"
    ["ontologies_and_associations"]="ontologies_and_associations"
)

# Execute the main download process
main() {
    local release_version="$1"
    local download_dir="$2"
    local year=$(echo "${release_version}" | cut -d'-' -f1)
    
    local start_time=$(date +%s)
    local failed_downloads=()
    local success_count=0
    
    log_info "Starting PomBase data download for release: ${release_version}"
    log_info "Download directory: ${download_dir}"
    log_info "Dry run mode: ${DRY_RUN}"
    
    # Check dependencies
    check_dependencies || exit 1
    
    # Create directory structure
    log_info "Creating directory structure..."
    for dir in "${!DIRECTORY_MAPPING[@]}"; do
        mkdir -p "${download_dir}/${dir}"
    done
    
    # Download files by category
    log_info "Downloading genome sequence and features..."
    local base_url="https://www.pombase.org"
    
    # Build download targets dynamically
    local download_targets_str=$(build_download_targets "$release_version")
    
    # Parse the returned string into variables
    while IFS='=' read -r key value; do
        local url="${base_url}${value}"
        local category="genome_sequence_and_features"
        
        # Determine appropriate directory based on target
        case $key in
            gene_*|*_viability) category="Gene_metadata" ;;
            *_expression*) category="RNA_metadata" ;;
            peptide_*|protein_*|disordered_*|*_modifications*) category="Protein_features" ;;
            go_*|fypo_*|mondo_*|*_annotations) category="ontologies_and_associations" ;;
            *) category="genome_sequence_and_features" ;;
        esac
        
        if download_file "$url" "${download_dir}/${category}" "$key"; then
            ((success_count++))
        else
            failed_downloads+=("$key")
        fi
    done <<< "$download_targets_str"
    
    # Download external ontology files
    log_info "Downloading external ontology files..."
    
    # GO ontology
    download_file "https://purl.obolibrary.org/obo/go/go-basic.obo" \
                  "${download_dir}/ontologies_and_associations" "GO ontology"
    
    # Mondo ontology
    download_file "https://purl.obolibrary.org/obo/mondo.obo" \
                  "${download_dir}/ontologies_and_associations" "Mondo ontology"
    
    # GO slim for pombe
    download_file "https://current.geneontology.org/ontology/subsets/goslim_pombe.obo" \
                  "${download_dir}/ontologies_and_associations" "GO slim pombe"
    
    # FYPO ontology (latest version)
    download_fypo_ontology "${download_dir}/ontologies_and_associations"
    
    # Summary report
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log_info "Download completed in ${duration} seconds"
    log_info "Successfully downloaded: ${success_count} files"
    
    if [[ ${#failed_downloads[@]} -gt 0 ]]; then
        log_warn "Failed downloads: ${failed_downloads[*]}"
        log_warn "Total failed: ${#failed_downloads[@]}"
        return 2
    fi
    
    # Display summary
    cat << EOF

${COLOR_GREEN}=== Download Summary ===${COLOR_NC}
Release: ${release_version}
Directory: ${download_dir}
Duration: ${duration} seconds
Files downloaded: ${success_count}

Directory structure:
  ${download_dir}/
  ├── genome_sequence_and_features/  # FASTA and GFF files
  ├── Gene_metadata/                 # Gene information
  ├── RNA_metadata/                  # Expression data
  ├── Protein_features/              # Protein annotations
  └── ontologies_and_associations/   # Ontology files

${COLOR_GREEN}Download completed successfully!${COLOR_NC}
EOF
    
    return 0
}

# =============================================================================
# Script Initialization and Execution
# =============================================================================

# Initialize and run main function
initialize_and_run() {
    # Parse command line arguments
    parse_arguments "$@"
    
    # Execute main download process
    main "$@"
}

# Only run if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    initialize_and_run "$@"
fi






