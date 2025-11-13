#!/usr/bin/env bash

#===============================================================================
# STRING Database Data Fetching Script
#===============================================================================
#
# OVERVIEW:
#   Downloads comprehensive protein interaction datasets from STRING database
#   for Schizosaccharomyces pombe (taxon ID: 284812). This script automates
#   the acquisition of protein-protein interaction networks, protein information,
#   sequences, annotations, and orthology data with error recovery and validation.
#
# DATA SOURCES:
#   - STRING Database v12.0: Protein interaction networks and scores
#   - Protein Information: Basic protein metadata and descriptions
#   - Protein Sequences: Amino acid sequences in FASTA format
#   - Protein Aliases: Alternative protein identifiers and cross-references
#   - Homology/Orthology: Evolutionary relationships and orthologous proteins
#
# USAGE:
#   ./fetch_stringdb_datasets.sh [taxon_id] [download_directory]
#
# PARAMETERS:
#   taxon_id         - NCBI taxonomy ID for target organism
#                      (optional, default: 284812 for S. pombe)
#   download_directory - Target directory for downloaded files
#                       (optional, default: ./stringdb_datasets_<taxon_id>)
#
# EXAMPLES:
#   # Download S. pombe data to default location
#   ./fetch_stringdb_datasets.sh
#
#   # Download S. pombe data to specific directory
#   ./fetch_stringdb_datasets.sh "284812" "/data/stringdb"
#
#   # Download human protein data
#   ./fetch_stringdb_datasets.sh "9606"
#
#   # Download yeast protein data to custom directory
#   ./fetch_stringdb_datasets.sh "4932" "/data/yeast_string"
#
# OUTPUT STRUCTURE:
#   <download_dir>/
#   ├── protein_interactions/           # Protein-protein interaction networks
#   │   ├── protein.links.v12.0.txt.gz              # Combined interaction scores
#   │   ├── protein.links.detailed.v12.0.txt.gz     # Detailed interaction evidence
#   │   ├── protein.links.full.v12.0.txt.gz         # Complete interaction data
#   │   ├── protein.physical.links.v12.0.txt.gz     # Physical interactions only
#   │   ├── protein.physical.links.detailed.v12.0.txt.gz  # Detailed physical interactions
#   │   └── protein.physical.links.full.v12.0.txt.gz      # Complete physical interactions
#   ├── protein_information/           # Protein metadata and sequences
#   │   ├── protein.info.v12.0.txt.gz               # Basic protein information
#   │   ├── protein.sequences.v12.0.fa.gz           # Protein sequences (FASTA)
#   │   └── protein.aliases.v12.0.txt.gz            # Protein aliases and cross-refs
#   └── evolutionary_data/             # Homology and orthology information
#       ├── protein.homology.v12.0.txt.gz           # Homologous protein relationships
#       └── protein.orthology.v12.0.txt.gz          # Orthologous protein relationships
#
# ERROR HANDLING:
#   - Validates directory creation and permissions
#   - Verifies file downloads and non-zero file sizes
#   - Implements retry logic for failed downloads
#   - Creates directories with proper error checking
#   - Provides comprehensive logging and progress tracking
#
# REQUIREMENTS:
#   - wget: for file downloads
#   - Internet connection to STRING database servers
#   - Sufficient disk space for downloaded datasets (~500MB compressed)
#
# VERSION: 2.0
# UPDATED: $(date +%Y-%m-%d)
#===============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

#===============================================================================
# CONFIGURATION AND CONSTANTS
#===============================================================================

# Script identification
readonly SCRIPT_NAME="$(basename "${0}")"

# Default taxonomy ID for STRING database (can be overridden via command line)
# S. pombe: 284812, S. cerevisiae: 4932, Human: 9606, etc.
readonly DEFAULT_TAXON_ID="284812"

# STRING database version and base URL
readonly STRING_VERSION="v12.0"
readonly STRING_BASE_URL="https://stringdb-downloads.org/download"

# Wget configuration options
# -nc: no-clobber (don't overwrite existing files)
# --progress=bar:force: show progress bar always
# --timeout=30: timeout after 30 seconds
# --tries=3: retry up to 3 times
readonly WGET_OPTS="-nc --progress=bar:force --timeout=30 --tries=3"

#===============================================================================
# LOGGING AND OUTPUT FUNCTIONS
#===============================================================================

# Print formatted informational messages with timestamps
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*" >&2
}

# Print formatted warning messages with timestamps
log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*" >&2
}

# Print formatted error messages with timestamps
log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

#===============================================================================
# ARGUMENT PARSING AND VALIDATION
#===============================================================================

# Display usage information with detailed examples
show_usage() {
    cat << EOF
Usage: ${SCRIPT_NAME} [taxon_id] [download_directory]

Download comprehensive STRING database datasets for specified organism.

Arguments:
    taxon_id          NCBI taxonomy ID (optional, default: ${DEFAULT_TAXON_ID} for S. pombe)
                      Examples: 284812 (S. pombe), 4932 (S. cerevisiae), 9606 (Human)
    download_directory  Target directory (optional, default: ./stringdb_datasets_<taxon_id>)

Examples:
    # Download S. pombe data to default location
    ${SCRIPT_NAME}

    # Download human protein data
    ${SCRIPT_NAME} "9606"

    # Download to custom directory
    ${SCRIPT_NAME} "4932" "/data/yeast_string"

    # For help information
    ${SCRIPT_NAME} --help
    ${SCRIPT_NAME} -h

For detailed documentation about data sources and output structure,
see the script header documentation above.

Downloaded datasets include:
- Protein-protein interaction networks (multiple detail levels)
- Physical interaction networks (experimental evidence only)
- Protein information, sequences, and aliases
- Homology and orthology relationships
EOF
}

# Create directory with comprehensive error handling
create_directory() {
    local dir_path="$1"
    local description="$2"

    log_info "Creating ${description} directory: ${dir_path}"

    # Check if directory already exists
    if [[ -d "${dir_path}" ]]; then
        log_warn "Directory already exists: ${dir_path}"
        return 0
    fi

    # Create directory with parent directories if needed
    if ! mkdir -p "${dir_path}"; then
        log_error "Failed to create ${description} directory: ${dir_path}"
        log_error "Check permissions and available disk space"
        return 1
    fi

    log_info "Successfully created directory: ${dir_path}"
    return 0
}

# Download file with comprehensive error handling and validation
download_file() {
    local url="$1"
    local target_dir="$2"
    local description="${3:-file}"

    log_info "Downloading ${description}: $(basename "${url}")"

    # Download with retry logic and progress tracking
    if ! wget ${WGET_OPTS} -P "${target_dir}" "${url}"; then
        log_error "Failed to download ${description}: ${url}"
        log_error "Network error, timeout, or file not available"
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

    # Log success with file size
    local file_size
    file_size=$(du -h "${filepath}" | cut -f1)
    log_info "Successfully downloaded: ${filename} (${file_size})"
    return 0
}

#===============================================================================
# PROTEIN INTERACTION DOWNLOAD FUNCTIONS
#===============================================================================

# Download protein-protein interaction network files
download_protein_interactions() {
    local target_dir="$1"
    local base_url="$2"
    local taxon_id="$3"

    log_info "Starting protein-protein interaction downloads for taxon: ${taxon_id}..."

    local interactions_dir="${target_dir}/protein_interactions"
    create_directory "${interactions_dir}" "protein interactions" || return 1

    # Define interaction files to download
    local interaction_files=(
        "protein.links.${STRING_VERSION}/${taxon_id}.protein.links.${STRING_VERSION}.txt.gz"
        "protein.links.detailed.${STRING_VERSION}/${taxon_id}.protein.links.detailed.${STRING_VERSION}.txt.gz"
        "protein.links.full.${STRING_VERSION}/${taxon_id}.protein.links.full.${STRING_VERSION}.txt.gz"
        "protein.physical.links.${STRING_VERSION}/${taxon_id}.protein.physical.links.${STRING_VERSION}.txt.gz"
        "protein.physical.links.detailed.${STRING_VERSION}/${taxon_id}.protein.physical.links.detailed.${STRING_VERSION}.txt.gz"
        "protein.physical.links.full.${STRING_VERSION}/${taxon_id}.protein.physical.links.full.${STRING_VERSION}.txt.gz"
    )

    local descriptions=(
        "combined protein interaction scores"
        "detailed protein interaction evidence"
        "complete protein interaction data"
        "physical interaction scores"
        "detailed physical interaction evidence"
        "complete physical interaction data"
    )

    # Download each interaction file
    for i in "${!interaction_files[@]}"; do
        local file="${interaction_files[i]}"
        local desc="${descriptions[i]}"

        download_file \
            "${base_url}/${file}" \
            "${interactions_dir}" \
            "${desc}" || return 1
    done

    log_info "Protein interactions download completed successfully"
    return 0
}

#===============================================================================
# PROTEIN INFORMATION DOWNLOAD FUNCTIONS
#===============================================================================

# Download protein information, sequences, and aliases
download_protein_information() {
    local target_dir="$1"
    local base_url="$2"
    local taxon_id="$3"

    log_info "Starting protein information downloads for taxon: ${taxon_id}..."

    local info_dir="${target_dir}/protein_information"
    create_directory "${info_dir}" "protein information" || return 1

    # Define protein information files
    local info_files=(
        "protein.info.${STRING_VERSION}/${taxon_id}.protein.info.${STRING_VERSION}.txt.gz"
        "protein.sequences.${STRING_VERSION}/${taxon_id}.protein.sequences.${STRING_VERSION}.fa.gz"
        "protein.aliases.${STRING_VERSION}/${taxon_id}.protein.aliases.${STRING_VERSION}.txt.gz"
    )

    local descriptions=(
        "protein basic information"
        "protein sequences (FASTA)"
        "protein aliases and cross-references"
    )

    # Download each information file
    for i in "${!info_files[@]}"; do
        local file="${info_files[i]}"
        local desc="${descriptions[i]}"

        download_file \
            "${base_url}/${file}" \
            "${info_dir}" \
            "${desc}" || return 1
    done

    log_info "Protein information download completed successfully"
    return 0
}

#===============================================================================
# EVOLUTIONARY DATA DOWNLOAD FUNCTIONS
#===============================================================================

# Download homology and orthology data
download_evolutionary_data() {
    local target_dir="$1"
    local base_url="$2"
    local taxon_id="$3"

    log_info "Starting evolutionary data downloads for taxon: ${taxon_id}..."

    local evo_dir="${target_dir}/evolutionary_data"
    create_directory "${evo_dir}" "evolutionary data" || return 1

    # Define evolutionary data files
    local evo_files=(
        "protein.homology.${STRING_VERSION}/${taxon_id}.protein.homology.${STRING_VERSION}.txt.gz"
        "protein.orthology.${STRING_VERSION}/${taxon_id}.protein.orthology.${STRING_VERSION}.txt.gz"
    )

    local descriptions=(
        "protein homology relationships"
        "protein orthology relationships"
    )

    # Download each evolutionary data file
    for i in "${!evo_files[@]}"; do
        local file="${evo_files[i]}"
        local desc="${descriptions[i]}"

        download_file \
            "${base_url}/${file}" \
            "${evo_dir}" \
            "${desc}" || return 1
    done

    log_info "Evolutionary data download completed successfully"
    return 0
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

# Validate taxonomy ID format (positive integer)
validate_taxon_id() {
    local taxon_id="$1"

    # Use regex to validate positive integer format
    if [[ ! "${taxon_id}" =~ ^[0-9]+$ ]]; then
        log_error "Invalid taxonomy ID format: '${taxon_id}'"
        log_error "Expected format: positive integer (e.g., 284812 for S. pombe)"
        return 1
    fi

    # Basic validation for reasonable taxonomy ID range
    if (( taxon_id < 1 || taxon_id > 100000000 )); then
        log_error "Invalid taxonomy ID range: ${taxon_id}"
        log_error "Expected range: 1-100000000"
        return 1
    fi

    return 0
}

main() {
    # Parse and validate command-line arguments
    local taxon_id="${1:-${DEFAULT_TAXON_ID}}"
    local download_dir="${2:-./stringdb_datasets_${taxon_id}}"

    # Validate taxonomy ID format
    validate_taxon_id "${taxon_id}" || {
        show_usage
        exit 1
    }

    log_info "Starting STRING database download for taxon: ${taxon_id}"
    log_info "STRING version: ${STRING_VERSION}"
    log_info "Target directory: ${download_dir}"

    # Create main download directory
    create_directory "${download_dir}" "main download" || exit 1

    # Execute download phases in sequence
    log_info "Beginning multi-phase download process..."

    # Phase 1: Download protein interactions
    download_protein_interactions "${download_dir}" "${STRING_BASE_URL}" "${taxon_id}" || exit 1

    # Phase 2: Download protein information
    download_protein_information "${download_dir}" "${STRING_BASE_URL}" "${taxon_id}" || exit 1

    # Phase 3: Download evolutionary data
    download_evolutionary_data "${download_dir}" "${STRING_BASE_URL}" "${taxon_id}" || exit 1

    # Generate comprehensive completion summary
    log_info "All downloads completed successfully!"
    echo ""
    echo "==================== DOWNLOAD SUMMARY ===================="
    echo "STRING Version: ${STRING_VERSION}"
    echo "Target Taxonomy ID: ${taxon_id}"
    echo "Download Location: ${download_dir}"
    echo ""
    echo "Directory Structure:"
    echo "  ${download_dir}/"
    echo "  ├── protein_interactions/           # PPI networks and scores"
    echo "  │   ├── protein.links.*            # Combined interactions"
    echo "  │   └── protein.physical.links.*   # Physical interactions"
    echo "  ├── protein_information/           # Protein metadata"
    echo "  │   ├── protein.info.*            # Basic protein info"
    echo "  │   ├── protein.sequences.*       # Amino acid sequences"
    echo "  │   └── protein.aliases.*        # Alternative identifiers"
    echo "  └── evolutionary_data/             # Homology/orthology"
    echo "      ├── protein.homology.*        # Homologous relationships"
    echo "      └── protein.orthology.*       # Orthologous relationships"
    echo ""
    echo "Total downloaded files: $(find "${download_dir}" -type f | wc -l)"
    echo "Total download size: $(du -sh "${download_dir}" | cut -f1)"
    echo "=========================================================="
}

#===============================================================================
# SCRIPT EXECUTION
#===============================================================================

# Script entry point with comprehensive argument validation
if [[ $# -gt 2 ]]; then
    log_error "Invalid number of arguments (expected 0-2, got $#)"
    show_usage
    exit 1
fi

# Handle help flags
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_usage
    exit 0
fi

# Validate required external tools before starting download
log_info "Checking required tools..."
for tool in wget; do
    if ! command -v "${tool}" &> /dev/null; then
        log_error "Required tool '${tool}' is not installed or not in PATH"
        log_error "Please install ${tool} and ensure it's available in your PATH"
        exit 1
    fi
    log_info "✓ ${tool} is available"
done

# Execute main function with all arguments
main "$@"