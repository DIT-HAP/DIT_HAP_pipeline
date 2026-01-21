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

DIT-HAP (Diploid for Insertional Mutagenesis by Transposon and Haploid for Analysis of Phenotype) is a Snakemake workflow for analyzing transposon insertion sequencing data. The workflow processes paired-end sequencing reads to identify piggyBac transposon insertion sites, performs quality control, and conducts depletion analysis to understand gene essentiality. The project uses *Schizosaccharomyces pombe* (fission yeast) as the model organism with data sourced from PomBase.

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

# Update environments (recommended for faster resolution)
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only

# Use mamba for all operations (faster dependency resolution)
snakemake --use-conda --conda-frontend mamba --cores [number_of_cores]
```

### Development and Debugging
```bash
# Enable verbose logging and shell commands
snakemake --use-conda --cores 8 --printshellcmds --reason

# Keep temporary files for debugging
snakemake --use-conda --cores 8 --notemp

# Show detailed error messages and failed logs
snakemake --use-conda --cores 8 --show-failed-logs

# Resume after interruption with restart attempts
snakemake --use-conda --cores 16 --restart-times 2

# Create detailed log files with timestamp
snakemake --use-conda --cores 16 --log logs/snakemake_$(date +%Y%m%d_%H%M%S).log
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
- **`workflow/src/plot.py`**: Plotting functions with DIT-HAP styling (uses custom `config/DIT_HAP.mplstyle`)
- **`workflow/src/enrichment_functions.py`**: Statistical enrichment analysis functions
- **`workflow/src/pombe_feature_functions.py`**: S. pombe-specific feature analysis functions
- **`workflow/src/protein_structure_functions.py`**: Protein structure analysis functions
- **`workflow/scripts/`**: Individual analysis scripts for specific workflow steps

### Custom Plotting Style
The project uses a custom matplotlib style (`config/DIT_HAP.mplstyle`) with:
- **Arial font family** with specific sizing (title: 24pt, axis labels: 20pt, ticks: 18pt)
- **Custom color cycle**: `['dd8369', '6b99df', '98a64e', 'a78bd9', '64af6d', 'd57fbd', '4bb29c', 'e0788f', '4aadce', 'c4954b']`
- **Clean aesthetics**: No top/right spines, subtle grid (alpha: 0.3), publication-ready DPI (300)
- **Legend positioning**: Outside plot area at `loc: 1.05,0.3` by default

### Conda Environments
Modular conda environments in `workflow/envs/` for each major component:
- **FastQC/Fastp**: Read quality control and preprocessing (`fastqc.yml`, `fastp.yml`)
- **Cutadapt**: Adapter trimming and junction classification (`cutadapt.yml`)
- **BWA mapping**: Read alignment to reference genome (`bwa_mapping.yml`)
- **Statistics**: R and Python packages for depletion analysis (`statistics_and_figure_plotting.yml`, `pydeseq2.yml`)
- **Bioinformatics**: Biopython, pybedtools, pysam for biological data processing (`biopython.yml`, `pybedtools.yml`, `pysam.yml`)
- **MultiQC**: Report generation (`multiqc.yml`)

## Python Development Guidelines

When writing Python scripts for this workflow, follow these standards based on the project's Cursor rules:

### Code Structure and Standards
- **Use `pathlib.Path`** instead of `os.path` for all file operations
- **Use `logging` module** instead of print statements for output (except final results)
- **Use `argparse`** for command-line argument parsing with clear help messages
- **Prefer f-strings** for string formatting
- **Use context managers** (`with` statements) for file operations
- **Leverage scientific library functions** over manual loops: use numpy, pandas, scipy, sklearn built-in functions
- **Prefer vectorized operations** (numpy arrays, pandas operations) over iterative approaches
- **Use type hints** for all function parameters and return values

### Function Design and Documentation
- **Comprehensive module docstring** describing script purpose, input/output, and usage
- **Clear function organization** with descriptive names and single responsibilities
- **Function docstrings** with Args, Returns, and brief description
- **Meaningful variable names** that reflect biological/computational concepts
- **Type hints**: `from typing import List, Dict, Tuple, Optional, Union`

### Statistics and Logging
- **Set up logging** with appropriate levels (INFO, DEBUG, WARNING, ERROR)
- **Track key metrics**: processing time, record counts, success/failure rates, biologically relevant metrics
- **Store statistics** in structured format (dict or dataclass)
- **Present final statistics** as formatted table using tabulate library
- **Log processing steps** and intermediate results

### Example Template Structure
```python
#!/usr/bin/env python3
"""
Brief description of the bioinformatics script.

This script processes [input type] and generates [output type].
Typical usage: python script.py -i input.file -o output.file

Author: [Your name]
Date: [Date]
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional
import argparse

def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""

def process_data(input_file: Path) -> Dict[str, int]:
    """Main data processing function using vectorized operations."""

def main() -> None:
    """Main execution function with proper error handling."""

if __name__ == "__main__":
    main()
```

### Plotting Guidelines
Follow the project's custom plotting style (see Custom Plotting Style section above). Key principles:
- **Use the custom color palette** from `config/DIT_HAP.mplstyle`
- **Remove top and right spines** for cleaner appearance
- **Set 300 DPI** for publication quality
- **Include statistical annotations** when relevant (correlation, p-values, sample size)
- **Position legends outside plot area** to avoid data obscuration
- **Use consistent styling** across related figures

## Important Notes

### Configuration Management
- Always specify the correct config file for your experimental condition (HD/LD, generation stage, etc.)
- Sample sheets must match the expected format with Sample, Timepoint, Condition columns
- Time points for depletion analysis must be manually specified in config files
- **Default working directory**: `/data/c/yangyusheng_optimized/DIT_HAP_pipeline`
- **Active configuration**: `config/config_HD_generationRAW.yaml` (currently set in Snakefile:15)

### Available Configuration Files
- `config/config_HD_generationRAW.yaml` - High density, raw analysis (active)
- `config/config_HD_generationPLUS1.yaml` - High density, generation +1 analysis
- `config/config_LD_generationRAW.yaml` - Low density, raw analysis
- `config/config_LD_generationPLUS1.yaml` - Low density, generation +1 analysis
- `config/config_HD_diploid.yaml` - High density, diploid organism
- `config/config_LD_haploid.yaml` - Low density, haploid organism
- `config/config_spikein.yaml` - Spike-in analysis
- `config/config_1328_spore2YES6.yaml` - Specific experiment configuration

### Data Requirements
- Paired-end sequencing reads with piggyBac junction sequences
- Reference genome and annotation files (automatically downloaded from PomBase)
- Sufficient disk space for intermediate files (consider using `temp()` for temporary outputs)

### Analysis Modes
- **High Density (HD)** vs **Low Density (LD)** insertion experiments
- **Generation vs. Raw** analysis modes
- **Haploid vs. Diploid** organism configurations
- **With/without biological replicates** (DEseq2 vs. no-replicate analysis)

### Snakemake Workflow Organization
The main `Snakefile` uses a modular structure with wildcard constraints based on sample sheet:
- **Wildcard constraints**: Automatically generated from sample sheet for sample, timepoint, condition
- **Sample sheet dict**: Nested dictionary structure `sample[timepoint][condition][fq1/fq2]` for file paths
- **Rule organization**: Four main modules included from `workflow/rules/`
- **Target rules**: Currently commented out specific outputs for flexible execution

### Output Structure
- `results/{project_name}/`: Main analysis results organized by workflow step (01-13 sequential numbering)
- `reports/{project_name}/`: Quality control reports, visualizations, and statistics
- `logs/{project_name}/`: Log files organized by workflow stage
- `resources/pombase_data/{release_version}/`: Automatically downloaded reference data