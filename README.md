# DIT-HAP Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/DIT-HAP/DIT_HAP_pipeline)
[![zread](https://img.shields.io/badge/Ask_Zread-_.svg?style=flat&color=00b0aa&labelColor=000000&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB3aWR0aD0iMTYiIGhlaWdodD0iMTYiIHZpZXdCb3g9IjAgMCAxNiAxNiIgZmlsbD0ibm9uZSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KPHBhdGggZD0iTTQuOTYxNTYgMS42MDAxSDIuMjQxNTZDMS44ODgxIDEuNjAwMSAxLjYwMTU2IDEuODg2NjQgMS42MDE1NiAyLjI0MDFWNC45NjAxQzEuNjAxNTYgNS4zMTM1NiAxLjg4ODEgNS42MDAxIDIuMjQxNTYgNS42MDAxSDQuOTYxNTZDNS4zMTUwMiA1LjYwMDEgNS42MDE1NiA1LjMxMzU2IDUuNjAxNTYgNC45NjAxVjIuMjQwMUM1LjYwMTU2IDEuODg2NjQgNS4zMTUwMiAxLjYwMDEgNC45NjE1NiAxLjYwMDFaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik00Ljk2MTU2IDEwLjM5OTlIMi4yNDE1NkMxLjg4ODEgMTAuMzk5OSAxLjYwMTU2IDEwLjY4NjQgMS42MDE1NiAxMS4wMzk5VjEzLjc1OTlDMS42MDE1NiAxNC4xMTM0IDEuODg4MSAxNC4zOTk5IDIuMjQxNTYgMTQuMzk5OUg0Ljk2MTU2QzUuMzE1MDIgMTQuMzk5OSA1LjYwMTU2IDE0LjExMzQgNS42MDE1NiAxMy43NTk5VjExLjAzOTlDNS42MDE1NiAxMC42ODY0IDUuMzE1MDIgMTAuMzk5OSA0Ljk2MTU2IDEwLjM5OTlaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik0xMy43NTg0IDEuNjAwMUgxMS4wMzg0QzEwLjY4NSAxLjYwMDEgMTAuMzk4NCAxLjg4NjY0IDEwLjM5ODQgMi4yNDAxVjQuOTYwMUMxMC4zOTg0IDUuMzEzNTYgMTAuNjg1IDUuNjAwMSAxMS4wMzg0IDUuNjAwMUgxMy43NTg0QzE0LjExMTkgNS42MDAxIDE0LjM5ODQgNS4zMTM1NiAxNC4zOTg0IDQuOTYwMVYyLjI0MDFDMTQuMzk4NCAxLjg4NjY0IDE0LjExMTkgMS42MDAxIDEzLjc1ODQgMS42MDAxWiIgZmlsbD0iI2ZmZiIvPgo8cGF0aCBkPSJNNCAxMkwxMiA0TDQgMTJaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik00IDEyTDEyIDQiIHN0cm9rZT0iI2ZmZiIgc3Ryb2tlLXdpZHRoPSIxLjUiIHN0cm9rZS1saW5lY2FwPSJyb3VuZCIvPgo8L3N2Zz4K&logoColor=ffffff)](https://zread.ai/DIT-HAP/DIT_HAP_pipeline)


**DIT-HAP (Diploid for Insertional Mutagenesis by Transposon and Haploid for Analysis of Phenotype)** is a comprehensive Snakemake workflow for analyzing piggyBac transposon insertion sequencing data. The pipeline processes paired-end sequencing reads to identify transposon insertion sites, performs quality control, and conducts depletion analysis to understand gene essentiality.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Configuration](#configuration)
- [Workflow Architecture](#workflow-architecture)
- [Output Structure](#output-structure)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## Overview

DIT-HAP is designed for genome-wide transposon mutagenesis experiments using the piggyBac transposon system. The workflow:

1. **Processes paired-end sequencing reads** containing piggyBac junction sequences
2. **Separates reads** by junction orientation (PBL/PBR)
3. **Maps reads** to the reference genome with high-stringency filtering
4. **Identifies insertion sites** and annotates them with gene information
5. **Performs depletion analysis** with curve fitting to determine gene essentiality
6. **Generates comprehensive reports** for quality control and analysis

## Features

- **Modular Architecture**: Organized into distinct processing stages for easy maintenance
- **High Stringency Filtering**: Configurable filtering thresholds for mapping quality, alignment, and read pairing
- **Multiple Analysis Modes**: Support for HD/LD density experiments, generation vs. raw analysis, and haploid/diploid organisms
- **Curve Fitting Models**: Logistic, Richards, and sigmoid models for depletion analysis
- **Biological Replicates**: Optional DEseq2 integration for replicate analysis
- **Comprehensive QC**: MultiQC integration and custom quality control reports
- **Automated Data Retrieval**: Automatic download of reference genomes and annotations from PomBase

## Installation

### Prerequisites

- **Python**: ≥ 3.8
- **Snakemake**: ≥ 8.0.0
- **Conda**: Miniconda or Anaconda
- **System Requirements**: Linux/Unix environment with ≥ 8GB RAM recommended

### Clone the Repository

```bash
git clone https://github.com/your-username/DIT-HAP_pipeline.git
cd DIT_HAP_pipeline
```

### Set Up Conda Environments

The workflow uses modular conda environments that are automatically created:

```bash
# Create all environments at once
snakemake --use-conda --conda-create-envs-only

# Or use mamba for faster installation
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only
```

### Manual Environment Setup (Optional)

If you prefer to install dependencies manually:

```bash
# Core environment
conda create -n dit-hap -c conda-forge -c bioconda python=3.10 snakemake=8.0.0 pandas

# Additional environments will be created automatically when needed
```

## Quick Start

### 1. Prepare Sample Sheet

Create a tab-separated sample sheet with the following format:

```tsv
Sample	Timepoint	Condition	read1	read2
sample1	YES0	wildtype	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
sample1	YES3	wildtype	/path/to/sample1_t3_R1.fastq.gz	/path/to/sample1_t3_R2.fastq.gz
```

### 2. Configure the Workflow

Copy and modify an existing configuration file:

```bash
cp config/config_HD_generationPLUS1.yaml config/config_my_experiment.yaml
```

Edit the configuration file to set:
- `project_name`: Your experiment name
- `sample_sheet`: Path to your sample sheet
- `time_points`: Time points for depletion analysis
- Adapter sequences for your piggyBac system

### 3. Run the Workflow

```bash
# Dry run to check the workflow
snakemake -n --use-conda

# Run the full workflow
snakemake --use-conda --cores 16

# Run with specific configuration
snakemake --configfile config/config_my_experiment.yaml --use-conda --cores 16
```

## Usage

### Basic Commands

```bash
# Run entire workflow
snakemake --use-conda --cores [number_of_cores]

# Run with specific configuration
snakemake --configfile config/config_HD_generationPLUS1.yaml --use-conda --cores 16

# Create workflow visualization
snakemake --dag | dot -Tpdf > workflow_diagram.pdf

# Check workflow integrity
snakemake --validate

# Run specific rules
snakemake fastp_preprocessing --use-conda --cores 8
snakemake bwa_mapping --use-conda --cores 8
snakemake multiqc --use-conda --cores 4
```

### Advanced Options

```bash
# Use temporary directory for intermediate files
export TMPDIR=/path/to/temp
snakemake --use-conda --cores 16

# Resume after interruption
snakemake --use-conda --cores 16 --restart-times 2

# Keep intermediate files for debugging
snakemake --use-conda --cores 16 --notemp

# Create detailed log files
snakemake --use-conda --cores 16 --log logs/snakemake_$(date +%Y%m%d_%H%M%S).log
```

## Configuration

### Configuration Files

Configuration files are located in the `config/` directory:

- `config_HD_generationPLUS1.yaml` - High density, generation +1 analysis
- `config_LD_generationPLUS1.yaml` - Low density, generation +1 analysis
- `config_HD_generationRAW.yaml` - High density, raw analysis
- `config_LD_generationRAW.yaml` - Low density, raw analysis
- `config_HD_diploid.yaml` - High density, diploid organism
- `config_LD_haploid.yaml` - Low density, haploid organism

### Key Configuration Parameters

```yaml
# Project settings
project_name: "my_experiment"
sample_sheet: "config/sample_sheet.tsv"

# Adapter sequences for piggyBac
adapter_sequence: "CTGTCTCTTATACACATCT"
PBL_adapter: "CATGCGTCAATTTTACGCAGACTATCTTTCTAGGG"
PBR_adapter: "ACGCATGATTATCTTTAACGTACGTCACAATATGATTATCTTTCTAGGG"

# Read filtering thresholds
aligned_read_filtering:
  read_1_filtering:
    mapq_threshold: 20
    nm_threshold: 3
  read_2_filtering:
    mapq_threshold: 40
    nm_threshold: 15

# Depletion analysis time points
time_points:
  - 0
  - 3.352
  - 6.588
  - 10.104
  - 13.480

# Biological replicates
use_DEseq2_for_biological_replicates: true
```

### Sample Sheet Format

The sample sheet must contain the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| Sample | Sample identifier | `sample1` |
| Timepoint | Time point or condition | `YES0`, `YES3` |
| Condition | Experimental condition | `wildtype`, `treatment` |
| read1 | Path to R1 fastq file | `/path/to/R1.fastq.gz` |
| read2 | Path to R2 fastq file | `/path/to/R2.fastq.gz` |

## Workflow Architecture

The DIT-HAP pipeline consists of four main modules:

### 1. Preparation (`workflow/rules/preparation.smk`)
- Downloads reference genome and annotations from PomBase
- Indexes genome for BWA mapping
- Creates samtools indices

### 2. Preprocessing (`workflow/rules/preprocessing.smk`)
- Quality control with Fastp
- Adapter trimming and junction classification (PBL/PBR separation)
- Read mapping with BWA-MEM2
- Read parsing and filtering
- Insertion site extraction

### 3. Depletion Analysis (`workflow/rules/depletion_analysis.smk`)
- Gene-level aggregation of insertion data
- Hard filtering of low-count insertions
- Curve fitting with multiple models (logistic, Richards, sigmoid)
- Statistical analysis of depletion curves

### 4. Quality Control (`workflow/rules/quality_control.smk`)
- MultiQC report generation
- Insertion density analysis
- Read count distribution analysis
- PBL/PBR correlation analysis
- Insertion orientation analysis

### Data Processing Flow

```
Raw Reads → Fastp QC → PBL/PBR Separation → BWA Mapping → Read Filtering →
Insertion Site Extraction → Gene Annotation → Hard Filtering → Curve Fitting →
Statistical Analysis → Reports
```

## Output Structure

Results are organized in a hierarchical directory structure:

```
DIT_HAP_pipeline/
├── results/{project_name}/
│   ├── 01_fastp_preprocessing/
│   ├── 02_cutadapt_junction_classification/
│   ├── 03_bwa_mapping/
│   ├── 04_parse_bam_to_tsv/
│   ├── 05_aligned_read_filtering/
│   ├── 06_extract_insertion_sites/
│   ├── 07_annotation_concatenation/
│   ├── 08_hard_filtering/
│   ├── 09_insertion_level_depletion_analysis/
│   ├── 10_gene_level_depletion_analysis/
│   ├── 11_insertion_level_curve_fitting/
│   ├── 12_gene_level_curve_fitting/
│   └── 13_final_results/
├── reports/{project_name}/
│   ├── multiqc/
│   ├── mapping_filtering_statistics/
│   ├── PBL_PBR_correlation_analysis/
│   ├── read_count_distribution_analysis/
│   ├── insertion_orientation_analysis/
│   └── insertion_density_analysis/
├── logs/{project_name}/
│   ├── preparation/
│   ├── preprocessing/
│   ├── depletion_analysis/
│   └── quality_control/
└── resources/pombase_data/{release_version}/
    ├── genome_sequence_and_features/
    ├── Gene_metadata/
    └── Protein_features/
```

### Key Output Files

- `results/{project_name}/13_final_results/` - Final depletion analysis results
- `reports/{project_name}/multiqc/quality_control_multiqc_report.html` - Comprehensive QC report
- `reports/{project_name}/mapping_filtering_statistics/mapping_filtering_statistics.tsv` - Mapping statistics
- `results/{project_name}/12_gene_level_curve_fitting/gene_level_fitting_statistics.tsv` - Gene-level curve fitting results

## Examples

### Example 1: High Density Generation +1 Analysis

```bash
# Use the HD generation +1 configuration
snakemake --configfile config/config_HD_generationPLUS1.yaml --use-conda --cores 24
```

### Example 2: Low Density Haploid Analysis

```bash
# Use the LD haploid configuration
snakemake --configfile config/config_LD_haploid.yaml --use-conda --cores 16
```

### Example 3: Running Specific Workflow Stages

```bash
# Run only preprocessing
snakemake bwa_mapping --use-conda --cores 12

# Run only quality control
snakemake multiqc --use-conda --cores 4

# Run only depletion analysis
snakemake gene_level_curve_fitting --use-conda --cores 8
```

### Example 4: Custom Time Point Analysis

For a custom experiment with different time points:

```yaml
# In your config file
time_points:
  - 0
  - 2
  - 4
  - 6
  - 8
  - 12
  - 24
```

## Troubleshooting

### Common Issues

1. **Memory Issues**: Reduce `chunk_size` in configuration or use fewer cores
2. **Conda Environment Errors**: Use `--conda-frontend mamba` for faster resolution
3. **Mapping Failures**: Check adapter sequences and reference genome integrity
4. **Sample Sheet Errors**: Ensure tab-separated format and correct column names

### Debug Mode

```bash
# Enable verbose logging
snakemake --use-conda --cores 8 --printshellcmds --reason

# Keep temporary files for inspection
snakemake --use-conda --cores 8 --notemp

# Print detailed error messages
snakemake --use-conda --cores 8 --show-failed-logs
```

### Log File Locations

- Snakemake logs: `logs/{project_name}/`
- Rule-specific logs: `logs/{project_name}/{stage}/`
- MultiQC report: `reports/{project_name}/multiqc/`

## Citation

If you use DIT-HAP in your research, please cite:

```
DIT-HAP: A comprehensive pipeline for piggyBac transposon insertion sequencing analysis
[Your Name et al., Year]
GitHub Repository: https://github.com/your-username/DIT-HAP_pipeline
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Support

For questions and support:

- Open an issue on GitHub
- Check the [Snakemake documentation](https://snakemake.readthedocs.io/)
- Review the workflow logs for detailed error messages

---

**DIT-HAP Pipeline** - Comprehensive analysis of piggyBac transposon insertion sequencing data