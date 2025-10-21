<!-- OPENSPEC:START -->
# OpenSpec Instructions

These instructions are for AI assistants working in this project.

Always open `@/openspec/AGENTS.md` when the request:
- Mentions planning or proposals (words like proposal, spec, change, plan)
- Introduces new capabilities, breaking changes, architecture shifts, or big performance/security work
- Sounds ambiguous and you need the authoritative spec before coding

Use `@/openspec/AGENTS.md` to learn:
- How to create and apply change proposals
- Spec format and conventions
- Project structure and guidelines

Keep this managed block so 'openspec update' can refresh the instructions.

<!-- OPENSPEC:END -->

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DIT-HAP (Diploid for Insertional Mutagenesis by Transposon and Haploid for Analysis of Phenotype) is a Snakemake workflow for analyzing transposon insertion sequencing data. The workflow processes paired-end sequencing reads to identify piggyBac transposon insertion sites, performs quality control, and conducts depletion analysis to understand gene essentiality.

## Common Commands

### Running the Workflow
```bash
# Run the entire workflow
snakemake --use-conda --cores [number_of_cores]

# Run with specific configuration file
snakemake --configfile config/config_HD_generationPLUS1.yaml --use-conda --cores [number_of_cores]

# Dry run to see what will be executed
snakemake -n --use-conda

# Create workflow graph
snakemake --dag | dot -Tpdf > dag.pdf
```

### Development and Testing
```bash
# Lint Snakefiles
snakemake --lint

# Validate workflow
snakemake --validate

# Run specific rules
snakemake fastp_preprocessing --use-conda --cores [number_of_cores]
snakemake bwa_mapping --use-conda --cores [number_of_cores]
```

### Environment Management
```bash
# Create all conda environments
snakemake --use-conda --conda-create-envs-only

# Update environments
snakemake --use-conda --conda-frontend mamba
```

## Architecture Overview

### Workflow Structure
The pipeline follows a modular Snakemake architecture with distinct processing stages:

1. **Preparation** (`workflow/rules/preparation.smk`): Downloads reference genome and annotations from PomBase
2. **Preprocessing** (`workflow/rules/preprocessing.smk`): Quality control, adapter trimming, junction classification, and read alignment
3. **Depletion Analysis** (`workflow/rules/depletion_analysis.smk`): Curve fitting and statistical analysis of gene essentiality
4. **Quality Control** (`workflow/rules/quality_control.smk`): Generates comprehensive QC reports and analyses

### Configuration System
- **Main config files**: Located in `config/` directory with specific YAML files for different experimental conditions
- **Sample management**: Uses tab-separated sample sheets specifying Sample, Timepoint, Condition, and read file paths
- **Dynamic configuration**: Runtime parameters like adapter sequences, filtering thresholds, and time points are configurable per experiment

### Data Processing Flow
1. **Raw reads** → Fastp preprocessing → Cutadapt junction classification (PBL/PBR separation)
2. **Classification** → BWA mapping → Read parsing and filtering
3. **Filtering** → Insertion site extraction → Annotation concatenation
4. **Analysis** → Hard filtering → Depletion analysis with curve fitting

### Key Components

#### Transposon System
- **piggyBac (PB)** transposon with left (PBL) and right (PBR) junction sequences
- **PBL/PBR adapters** defined in configuration for demultiplexing
- **Orientation-specific processing** for forward and reverse insertions

#### Analysis Modules
- **Insertion-level analysis**: Processes individual transposon insertion sites
- **Gene-level analysis**: Aggregates insertion data to gene-level statistics
- **Curve fitting**: Implements logistic, Richards, and sigmoid models for depletion analysis
- **Replicate handling**: Supports both DEseq2 for biological replicates and no-replicate analysis

#### Visualization and Reporting
- **MultiQC integration**: Consolidates quality control reports
- **Custom plotting**: Uses matplotlib with project-specific styling (`config/DIT_HAP.mplstyle`)
- **Statistical analysis**: Generates PDF reports with correlation analyses and distribution plots

### Python Modules
- **`workflow/src/utils.py`**: Core utility functions for file I/O and data manipulation
- **`workflow/src/plot.py`**: Plotting functions with DIT-HAP styling
- **`workflow/src/enrichment_functions.py`**: Statistical enrichment analysis functions
- **`workflow/scripts/`**: Individual analysis scripts for specific workflow steps

### Conda Environments
Modular conda environments in `workflow/envs/` for each major component:
- **FastQC/Fastp**: Read quality control and preprocessing
- **BWA mapping**: Read alignment to reference genome
- **Statistics**: R and Python packages for depletion analysis
- **MultiQC**: Report generation

## Important Notes

### Configuration Management
- Always specify the correct config file for your experimental condition (HD/LD, generation stage, etc.)
- Sample sheets must match the expected format with Sample, Timepoint, Condition columns
- Time points for depletion analysis must be manually specified in config files

### Data Requirements
- Paired-end sequencing reads with piggyBac junction sequences
- Reference genome and annotation files (automatically downloaded from PomBase)
- Sufficient disk space for intermediate files (consider using `temp()` for temporary outputs)

### Analysis Modes
- **High Density (HD)** vs **Low Density (LD)** insertion experiments
- **Generation vs. Raw** analysis modes
- **Haploid vs. Diploid** organism configurations
- **With/without biological replicates** (DEseq2 vs. no-replicate analysis)

### Output Structure
- `results/{project_name}/`: Main analysis results organized by workflow step
- `reports/{project_name}/`: Quality control reports, visualizations, and statistics
- `logs/{project_name}/`: Log files organized by workflow stage