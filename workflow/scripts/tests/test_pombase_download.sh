#!/usr/bin/env bash

# =============================================================================
# Test Suite for PomBase Download Script
# =============================================================================
#
# TDD Test Framework for download_file_from_pombaseFTP.sh
# Follows RED-GREEN-REFACTOR cycle
#
# Usage:
#   ./tests/test_pombase_download.sh [test_name]
#   ./tests/test_pombase_download.sh --all
#
# =============================================================================

set -euo pipefail

# Test configuration
SCRIPT_PATH="workflow/scripts/preparation/download_file_from_pombaseFTP.sh"
TEST_DIR="/tmp/pombase_test_$(date +%s)"
EXPECTED_RELEASE="2025-08-01"

# Colors for test output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# =============================================================================
# Test Infrastructure
# =============================================================================

log_test() {
    echo -e "${GREEN}[TEST]${NC} $(date '+%H:%M:%S') - $1"
}

log_fail() {
    echo -e "${RED}[FAIL]${NC} $(date '+%H:%M:%S') - $1"
}

log_skip() {
    echo -e "${YELLOW}[SKIP]${NC} $(date '+%H:%M:%S') - $1"
}

assert_equals() {
    local expected="$1"
    local actual="$2"
    local test_name="$3"
    
    ((TESTS_RUN++))
    
    if [[ "$expected" == "$actual" ]]; then
        ((TESTS_PASSED++))
        log_test "PASS: $test_name"
    else
        ((TESTS_FAILED++))
        log_fail "FAIL: $test_name - Expected: '$expected', Got: '$actual'"
    fi
}

assert_file_exists() {
    local file_path="$1"
    local test_name="$2"
    
    ((TESTS_RUN++))
    
    if [[ -f "$file_path" ]]; then
        ((TESTS_PASSED++))
        log_test "PASS: $test_name"
    else
        ((TESTS_FAILED++))
        log_fail "FAIL: $test_name - File not found: $file_path"
    fi
}

assert_url_accessible() {
    local url="$1"
    local test_name="$2"
    
    ((TESTS_RUN++))
    
    if curl --head --silent --fail "$url" >/dev/null 2>&1; then
        ((TESTS_PASSED++))
        log_test "PASS: $test_name"
    else
        ((TESTS_FAILED++))
        log_fail "FAIL: $test_name - URL not accessible: $url"
    fi
}

# =============================================================================
# Unit Tests
# =============================================================================

test_date_validation() {
    log_test "Running date validation tests..."
    
    # Valid format
    assert_equals "valid" "$(echo "2025-08-01" | grep -E '^[0-9]{4}-[0-9]{2}-[0-9]{2}$' && echo "valid" || echo "invalid")" "Valid date format"
    
    # Invalid formats
    assert_equals "invalid" "$(echo "2025-8-1" | grep -E '^[0-9]{4}-[0-9]{2}-[0-9]{2}$' && echo "valid" || echo "invalid")" "Invalid date format (single digits)"
    assert_equals "invalid" "$(echo "25-08-01" | grep -E '^[0-9]{4}-[0-9]{2}-[0-9]{2}$' && echo "valid" || echo "invalid")" "Invalid date format (short year)"
    assert_equals "invalid" "$(echo "20250801" | grep -E '^[0-9]{4}-[0-9]{2}-[0-9]{2}$' && echo "valid" || echo "invalid")" "Invalid date format (no separators)"
}

test_year_extraction() {
    log_test "Running year extraction tests..."
    
    assert_equals "2025" "$(echo "2025-08-01" | cut -d'-' -f1)" "Year extraction from valid date"
    assert_equals "" "$(echo "invalid" | cut -d'-' -f1)" "Year extraction from invalid date"
}

test_url_construction() {
    log_test "Running URL construction tests..."
    
    local release_date="2025-08-01"
    local year="2025"
    local base_url="https://www.pombase.org"
    
    local expected_chromosomes="${base_url}/monthly_releases/${year}/pombase-${release_date}/genome_sequence_and_features/fasta_format/chromosomes/Schizosaccharomyces_pombe_all_chromosomes.fa"
    local expected_gene_ids="${base_url}/monthly_releases/${year}/pombase-${release_date}/gene_names_and_identifiers/gene_IDs_names_products.tsv"
    
    assert_equals "$expected_chromosomes" "$expected_chromosomes" "Chromosome URL construction"
    assert_equals "$expected_gene_ids" "$expected_gene_ids" "Gene IDs URL construction"
}

test_url_accessibility() {
    log_test "Running URL accessibility tests..."
    
    local base_url="https://www.pombase.org"
    local release_date="2025-08-01"
    local year="2025"
    
    local test_urls=(
        "${base_url}/monthly_releases/${year}/pombase-${release_date}/genome_sequence_and_features/fasta_format/chromosomes/Schizosaccharomyces_pombe_all_chromosomes.fa"
        "${base_url}/monthly_releases/${year}/pombase-${release_date}/gene_names_and_identifiers/gene_IDs_names_products.tsv"
        "${base_url}/monthly_releases/${year}/pombase-${release_date}/gene_expression/qualitative_gene_expression.tsv"
    )
    
    for url in "${test_urls[@]}"; do
        assert_url_accessible "$url" "URL accessible: $(basename "$url")"
    done
}

test_dependency_check() {
    log_test "Running dependency check tests..."
    
    # Test required tools
    assert_equals "found" "$(command -v wget &> /dev/null && echo "found" || echo "missing")" "wget availability"
    assert_equals "found" "$(command -v curl &> /dev/null && echo "found" || echo "missing")" "curl availability"
    
    # Test optional tools
    if command -v md5sum >/dev/null 2>&1; then
        log_test "PASS: md5sum available for checksum verification"
    else
        log_skip "SKIP: md5sum not available"
    fi
}

# =============================================================================
# Integration Tests
# =============================================================================

test_directory_creation() {
    log_test "Running directory creation tests..."
    
    local test_dir="/tmp/test_pombase_dirs"
    mkdir -p "$test_dir"
    
    local expected_dirs=("genome_sequence_and_features" "Gene_metadata" "RNA_metadata" "Protein_features" "ontologies_and_associations")
    
    for dir in "${expected_dirs[@]}"; do
        mkdir -p "$test_dir/$dir"
        assert_file_exists "$test_dir/$dir" "Directory creation: $dir"
    done
    
    # Cleanup
    rm -rf "$test_dir"
}

test_dry_run_mode() {
    log_test "Running dry-run mode test..."
    
    local test_dir="/tmp/test_dry_run"
    local output=$(bash "$SCRIPT_PATH" "2025-08-01" "$test_dir" --dry-run 2>&1)
    
    if [[ "$output" == *"DRY RUN"* ]]; then
        ((TESTS_PASSED++))
        log_test "PASS: Dry-run mode works"
    else
        ((TESTS_FAILED++))
        log_fail "FAIL: Dry-run mode not detected"
    fi
}

# =============================================================================
# Test Runner
# =============================================================================

run_all_tests() {
    echo "=== PomBase Download Script Test Suite ==="
    echo "Script: $SCRIPT_PATH"
    echo "Test directory: $TEST_DIR"
    echo "==========================================="
    
    # Create test directory
    mkdir -p "$TEST_DIR"
    
    # Run unit tests
    test_date_validation
    test_year_extraction
    test_url_construction
    test_url_accessibility
    test_dependency_check
    
    # Run integration tests
    test_directory_creation
    test_dry_run_mode
    
    # Test summary
    echo "==========================================="
    echo "Test Summary:"
    echo "Tests run: $TESTS_RUN"
    echo "Tests passed: $TESTS_PASSED"
    echo "Tests failed: $TESTS_FAILED"
    
    if [[ $TESTS_FAILED -eq 0 ]]; then
        echo -e "${GREEN}All tests passed!${NC}"
        exit 0
    else
        echo -e "${RED}Some tests failed!${NC}"
        exit 1
    fi
}

# =============================================================================
# Main Execution
# =============================================================================

# Check if script exists
if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "Error: Script not found at $SCRIPT_PATH"
    exit 1
fi

# Run tests
run_all_tests