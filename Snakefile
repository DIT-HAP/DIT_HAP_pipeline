"""
DIT-HAP pipeline Snakefile for high-throughput insertion sequencing analysis.

This workflow processes high-throughput insertion sequencing data to identify
fitness effects of gene disruptions under various conditions and timepoints.

Pipeline Overview:
1. Data preparation and preprocessing
2. Quality control and filtering
3. Insertion site identification
4. Depletion analysis and curve fitting
5. Gene-level statistical analysis

Usage:
    snakemake --cores 32 --use-conda
    snakemake --cores 32 --use-conda --configfile config/custom_config.yaml

Configuration:
    All parameters are defined in config/config_generationPlus1.yaml
    See config/README.md for detailed parameter descriptions.

Outputs:
    results/{project_name}/ - All analysis results
    reports/{project_name}/ - Generated reports and visualizations
    logs/{project_name}/ - Execution logs and metrics
"""

from snakemake.utils import min_version
from pathlib import Path
import pandas as pd
import logging
from typing import Dict, List, Any

min_version("8.0")

# Configure logging for workflow monitoring
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Load configuration with validation
# -----------------------------------------------------
try:
    configfile: "config/config_generationPlus1.yaml"
    workdir: "/data/c/yangyusheng_optimized/DIT_HAP_pipeline"
except Exception as e:
    logger.error(f"Failed to load configuration: {e}")
    raise SystemExit(1)

# Event handlers for workflow lifecycle
# -----------------------------------------------------
onstart:
    """Log workflow initialization with configuration details."""
    logger.info("=" * 60)
    logger.info("DIT-HAP Pipeline Started")
    logger.info(f"Project: {config.get('project_name', 'Unknown')}")
    logger.info(f"Working directory: {Path.cwd()}")
    logger.info(f"Config file: config/config_generationPlus1.yaml")
    logger.info("=" * 60)

onsuccess:
    """Log successful workflow completion."""
    logger.info("=" * 60)
    logger.info("DIT-HAP Pipeline completed successfully!")
    logger.info("All analysis outputs are available in results/ directory")
    logger.info("=" * 60)

onerror:
    """Log detailed error information for debugging."""
    logger.error("=" * 60)
    logger.error("DIT-HAP Pipeline FAILED!")
    logger.error("Check logs for detailed error information")
    logger.error("=" * 60)

# Project configuration and sample handling
# -----------------------------------------------------
project_name = config["project_name"]

def load_sample_sheet(config_path: str) -> pd.DataFrame:
    """
    Load and validate sample sheet configuration.
    
    Parameters
    ----------
    config_path : str
        Path to sample sheet TSV file
        
    Returns
    -------
    pd.DataFrame
        Validated sample sheet with required columns
        
    Raises
    ------
    ValueError
        If required columns are missing or data is malformed
    FileNotFoundError
        If sample sheet file does not exist
    """
    required_columns = {"Sample", "Timepoint", "Condition", "read1", "read2"}
    
    try:
        sample_sheet = pd.read_csv(config_path, sep="\t")
        
        # Validate required columns
        missing_cols = required_columns - set(sample_sheet.columns)
        if missing_cols:
            raise ValueError(
                f"Sample sheet missing required columns: {missing_cols}"
            )
        
        # Validate file paths exist
        for _, row in sample_sheet.iterrows():
            for read_col in ["read1", "read2"]:
                file_path = Path(row[read_col])
                if not file_path.exists():
                    logger.warning(f"File not found: {file_path}")
        
        logger.info(f"Loaded sample sheet with {len(sample_sheet)} samples")
        return sample_sheet
        
    except FileNotFoundError:
        logger.error(f"Sample sheet not found: {config_path}")
        raise
    except Exception as e:
        logger.error(f"Error loading sample sheet: {e}")
        raise

# Load and validate sample data
sample_sheet = load_sample_sheet(config["sample_sheet"])
samples = sample_sheet["Sample"].unique().tolist()
timepoints = sample_sheet["Timepoint"].unique().tolist()
conditions = sample_sheet["Condition"].unique().tolist()

# Create validated sample dictionary structure
def create_sample_dictionary(
    sample_sheet: pd.DataFrame,
    samples: List[str],
    timepoints: List[str],
    conditions: List[str]
) -> Dict[str, Dict[str, Dict[str, Dict[str, Path]]]]:
    """
    Create nested dictionary structure for sample file organization.
    
    Parameters
    ----------
    sample_sheet : pd.DataFrame
        Validated sample sheet data
    samples : List[str]
        List of unique sample identifiers
    timepoints : List[str]
        List of unique timepoint identifiers
    conditions : List[str]
        List of unique condition identifiers
        
    Returns
    -------
    Dict
        Nested dictionary mapping sample -> timepoint -> condition -> read files
    """
    sample_dict = {
        sample: {
            timepoint: {
                condition: {"fq1": None, "fq2": None} 
                for condition in conditions
            } 
            for timepoint in timepoints
        } 
        for sample in samples
    }
    
    # Populate with actual file paths
    for _, row in sample_sheet.iterrows():
        sample = row["Sample"]
        timepoint = row["Timepoint"]
        condition = row["Condition"]
        
        if sample in sample_dict and timepoint in sample_dict[sample] and condition in sample_dict[sample][timepoint]:
            sample_dict[sample][timepoint][condition]["fq1"] = Path(row["read1"])
            sample_dict[sample][timepoint][condition]["fq2"] = Path(row["read2"])
        else:
            logger.warning(f"Skipping invalid sample entry: {sample}-{timepoint}-{condition}")
    
    return sample_dict

sample_sheet_dict = create_sample_dictionary(sample_sheet, samples, timepoints, conditions)

# Wildcard constraints for sample identification
# -----------------------------------------------------
wildcard_constraints:
    sample = "|".join(map(str, samples)),
    timepoint = "|".join(map(str, timepoints)),
    condition = "|".join(map(str, conditions))

# Target rule defining final workflow outputs
# -----------------------------------------------------
rule all:
    """
    Main target rule defining all final outputs of the DIT-HAP pipeline.
    
    This rule ensures all analysis steps are completed and final results
    are generated in the appropriate directory structure.
    """
    input:
        # Insertion-level depletion analysis results
        insertion_curve_fitting = f"results/{project_name}/16_insertion_level_curve_fitting/insertions_LFC_fitted.tsv",
        
        # Gene-level statistical analysis results
        gene_level_statistics = f"results/{project_name}/18_gene_level_curve_fitting/Gene_level_statistics_fitted.tsv",
        
        # Quality control reports (when QC rules are enabled)
        # directory(f"reports/{project_name}/multiqc_report.html"),
        
        # Optional outputs (commented out for current pipeline)
        # non_coding_genes = f"results/{project_name}/19_insertion_in_non_coding_genes/Non_coding_genes_Gene_level_statistics_fitted.tsv",
        # annotations = f"results/{project_name}/12_concatenated/annotations.tsv",
        # gene_coverage = directory(f"reports/{project_name}/gene_coverage_analysis"),
    output:
        # Completion marker file
        workflow_complete = touch(f"logs/{project_name}/workflow_complete.log"),
    run:
        """Log completion of all target outputs."""
        logger.info("All target outputs successfully generated")
        logger.info(f"Insertion curve fitting: {input.insertion_curve_fitting}")
        logger.info(f"Gene level statistics: {input.gene_level_statistics}")

# Workflow rule modules
# -----------------------------------------------------
# Each module handles a specific stage of the analysis pipeline
include: "workflow/rules/preparation.smk"      # Reference genome preparation and indexing
include: "workflow/rules/preprocessing.smk"    # Read preprocessing and quality control
include: "workflow/rules/depletion_analysis.smk" # Main depletion and curve fitting analysis
# include: "workflow/rules/quality_control.smk"  # Optional QC and reporting
# include: "workflow/rules/misc.smk"            # Miscellaneous utility rules

# Validation helpers
# -----------------------------------------------------
def validate_workflow_config() -> None:
    """
    Validate critical configuration parameters before execution.
    
    Raises
    ------
    ValueError
        If critical configuration parameters are missing or invalid
    """
    required_keys = [
        "project_name", "sample_sheet", "Pombase_release_version",
        "adapter_sequence", "r1_mapq_threshold", "hard_filtering_cutoff"
    ]
    
    missing_keys = [k for k in required_keys if k not in config]
    if missing_keys:
        raise ValueError(f"Missing required configuration keys: {missing_keys}")
    
    # Validate numeric parameters
    numeric_params = [
        "r1_mapq_threshold", "r1_ncigar_value", "r1_nm_threshold",
        "r2_mapq_threshold", "hard_filtering_cutoff"
    ]
    
    for param in numeric_params:
        if param in config and not isinstance(config[param], (int, float)):
            raise ValueError(f"Parameter {param} must be numeric")
    
    logger.info("Configuration validation completed successfully")

# Run configuration validation
validate_workflow_config()


