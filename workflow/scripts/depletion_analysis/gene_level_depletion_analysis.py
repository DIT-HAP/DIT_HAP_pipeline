#!/usr/bin/env python3
"""
Gene-Level Depletion Analysis for Transposon Insertion Sequencing

This script performs gene-level aggregation of transposon insertion depletion data,
combining multiple insertion sites within genes using weighted averages and statistical
meta-analysis methods. It processes log2 fold changes (LFC) and adjusted p-values
from insertion-level analysis to generate gene-level fitness scores.

The analysis pipeline includes:
1. Loading insertion-level LFC and p-value data
2. Filtering for in-gene insertions based on annotations
3. Calculating statistical weights from p-values
4. Aggregating insertions per gene using weighted averages
5. Combining p-values using Fisher's method or alternatives
6. Generating gene-level fitness profiles across time points

Statistical Methods:
- Fisher's method: Combines p-values assuming independence
- Stouffer's method: Z-score based combination with equal weights
- Brown's method: Fisher's method adjusted for correlation structure

Typical usage:
    python gene_level_depletion_analysis.py -l lfc_results.csv -a annotations.csv -o gene_results.csv

Input: CSV files with insertion-level LFC/p-values and genomic annotations
Output: CSV file with gene-level fitness scores and combined p-values

Author: Bioinformatics Pipeline
Date: 2024
"""

import logging
import numpy as np
import pandas as pd
import sys
import multiprocessing as mp
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import argparse
from functools import partial
from collections import defaultdict
import time
import warnings

import scipy.stats as stats
from scipy.stats import permutation_test, ConstantInputWarning
from sklearn.utils import resample

from loguru import logger
from pydantic import BaseModel, field_validator, Field

# ======================== Configuration & Models ========================

class GeneLevelDepletionAnalysisConfig(BaseModel):

    lfc_path: Path = Field(..., description="Path to the LFC file")
    weights_path: Path = Field(..., description="Path to the weights file")
    annotations_path: Path = Field(..., description="Path to the annotations file")

    # ======================== Validation ========================
    @field_validator("lfc_path", "weights_path", "annotations_path")
    def validate_input_exists(cls, v, field):
        """Validate that input files exist."""
        if not v.exists():
            raise ValueError(f"Input file {v} does not exist")
        return v

    class Config:
        frozen = True

# ======================== Logging Setup ========================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for gene-level depletion analysis."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )

# ======================== Core Functions ========================

@logger.catch
def load_data(lfc_path: Path, weights_path: Path, annotations_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Index]:
    """
    Load and validate insertion-level LFC data and genomic annotations.
    
    Args:
        lfc_path: Path to CSV file with LFC and p-value results
        annotations_path: Path to CSV file with insertion annotations
        
    Returns:
        Tuple of (results_df, annotations_df, in_gene_insertions_index)
        
    Raises:
        FileNotFoundError: If input files don't exist
        ValueError: If required columns are missing
    """
    logger.info(f"Loading LFC data from {lfc_path}")
    logger.info(f"Loading weights data from {weights_path}")
    logger.info(f"Loading annotations from {annotations_path}")
    
    if not lfc_path.exists():
        raise FileNotFoundError(f"LFC file not found: {lfc_path}")
    if not weights_path.exists():
        raise FileNotFoundError(f"Weights file not found: {weights_path}")
    if not annotations_path.exists():
        raise FileNotFoundError(f"Annotations file not found: {annotations_path}")
    
    try:
        # Load results with multi-level columns and index
        lfc_df = pd.read_csv(lfc_path, index_col=[0, 1, 2, 3], sep="\t")
        weights_df = pd.read_csv(weights_path, index_col=[0, 1, 2, 3], sep="\t")
        insertion_annotations = pd.read_csv(annotations_path, index_col=[0, 1, 2, 3], sep="\t")
        
        # Validate required columns exist
        
        # Filter for in-gene insertions
        in_gene_insertions = insertion_annotations.query(
            "(Type != 'Intergenic region') and (Distance_to_stop_codon > 4)"
        ).index
        in_gene_insertion_count = lfc_df[lfc_df.index.isin(in_gene_insertions)].shape[0]
        
        logging.info(f"Loaded {lfc_df.shape[0]} total insertions")
        logging.info(f"Found {in_gene_insertion_count} in-gene insertions")
        logging.info(f"Time points in data: {lfc_df.columns.unique().tolist()}")
        
        return lfc_df, weights_df, insertion_annotations, in_gene_insertions
        
    except pd.errors.EmptyDataError:
        raise ValueError("Input files are empty")
    except Exception as e:
        raise ValueError(f"Error loading data: {e}")

@logger.catch
def process_lfc_data(lfc_df: pd.DataFrame, weights_df: pd.DataFrame, in_gene_insertions: pd.Index) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract and filter LFC and p-value data for in-gene insertions.
    
    Args:
        results: Multi-level DataFrame with LFC and p-value data
        in_gene_insertions: Index of insertions located within genes
        
    Returns:
        Tuple of (in_gene_LFCs, in_gene_pvalues) DataFrames
    """
    logger.info("Processing LFC and p-value data")
    
    # Extract LFC and p-value data
    LFCs = lfc_df.copy().rename_axis("Timepoint", axis=1)
    weights = weights_df.copy().rename_axis("Timepoint", axis=1)
    
    # Filter for in-gene insertions
    in_gene_LFCs = LFCs[LFCs.index.isin(in_gene_insertions)].copy()
    in_gene_weights = weights[weights.index.isin(in_gene_insertions)].copy()
    
    # Data quality checks
    lfc_nan_count = in_gene_LFCs.isna().any(axis=1).sum()
    weight_nan_count = in_gene_weights.isna().any(axis=1).sum()
    
    if lfc_nan_count > 0:
        logger.warning(f"Found {lfc_nan_count} rows with NaN values in LFC data")
    if weight_nan_count > 0:
        logger.warning(f"Found {weight_nan_count} rows with NaN values in weights data")
        logger.warning(f"Fill NaN values in weights data with 1 since they are mainly due to the large cooks distance")
        logger.warning(f"Note on p-values set to NA: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#more-information-on-results-columns")
        in_gene_weights = in_gene_weights.fillna(1)
    
    logging.info(f"Processed {in_gene_LFCs.shape[0]} in-gene LFC measurements")
    
    return in_gene_LFCs, in_gene_weights


@logger.catch
def calculate_statistical_weights(GMs: pd.DataFrame, min_weight: float = 1e-10) -> pd.DataFrame:
    """
    Calculate statistical weights from adjusted p-values using -log10 transformation.
    
    Args:
        GMs: DataFrame with 'Padj' column
        min_pvalue: Minimum p-value to prevent infinite weights
        
    Returns:
        DataFrame with added 'Weights' column
    """
    logger.info("Calculating statistical weights from p-values")
    
    # Clip p-values to prevent infinite weights
    logger.info(f"Clipping weights to prevent infinite weights, min_weight={min_weight}, max_weight={1-min_weight}")
    clipped_weights = GMs["Weights"].clip(lower=min_weight, upper=1-min_weight)

    # Calculate -log10 weights
    GMs["Weights"] = -np.log10(clipped_weights)
    
    # Log weight statistics
    weight_stats = GMs["Weights"].describe()
    logger.info(f"Weight statistics: mean={weight_stats['mean']:.2f}, "
                f"median={weight_stats['50%']:.2f}, max={weight_stats['max']:.2f}")
    
    return GMs


@logger.catch
def merge_and_normalize_weights(GMs: pd.DataFrame, insertion_annotations: pd.DataFrame) -> pd.DataFrame:
    """
    Merge LFC data with annotations and normalize weights within genes.
    
    Args:
        GMs: DataFrame with LFC, p-values, and weights
        insertion_annotations: DataFrame with genomic annotations
        
    Returns:
        Annotated DataFrame with normalized weights per gene per timepoint
    """
    logger.info("Merging data with annotations and normalizing weights")
    
    # Merge with annotations
    GMs_annotated = pd.merge(
        GMs, insertion_annotations, 
        how="left", left_index=True, right_index=True
    )
    
    # Remove rows with missing weights
    initial_count = GMs_annotated.shape[0]
    GMs_annotated = GMs_annotated[GMs_annotated["Weights"].notna()].copy()
    final_count = GMs_annotated.shape[0]
    
    if initial_count != final_count:
        logger.warning(f"GMs_annotated: Removed {initial_count - final_count} rows with missing weights")
    
    # Normalize weights within each gene and timepoint
    GMs_annotated["Normalized_weights"] = GMs_annotated.groupby(
        ["Systematic ID", "Timepoint"]
    )["Weights"].transform(lambda x: x / x.sum())
    
    # Validate normalization
    weight_sums = GMs_annotated.groupby(["Systematic ID", "Timepoint"])["Normalized_weights"].sum()
    if not np.allclose(weight_sums, 1.0, rtol=1e-10):
        logger.warning("Weight normalization may have numerical issues")
    
    logger.info(f"GMs_annotated: Successfully annotated {GMs_annotated.shape[0]} insertions")
    
    return GMs_annotated


@logger.catch
def gene_level_lfc_and_weights(gene_df: pd.DataFrame) -> Tuple[List[float], List[float]]:
    """
    Calculate gene-level LFC and combined p-values across timepoints.
    
    Args:
        gene_df: DataFrame with insertion data for a single gene
        method: P-value combination method ('fisher', 'stouffer', 'brown')
        n_permutations: Number of permutations for permutation test
        
    Returns:
        Tuple of (gene_level_lfcs, gene_level_pvalues) across timepoints
    """
    # Remove duplicate transcript entries, keeping only .1 transcripts
    gene_df_clean = gene_df.reset_index().set_index(["Chr", "Coordinate", "Strand", "Target", "Timepoint"])
    
    gene_level_results = pd.DataFrame()
    
    for tp, tp_df in gene_df_clean.groupby("Timepoint"):
        lfcs = tp_df["LFC"].values
        weights = tp_df["Weights"].values
        normalized_weights = tp_df["Normalized_weights"].values
        
        # Calculate weighted average LFC
        observed_gene_level_lfc = np.average(lfcs, weights=normalized_weights)

        
        gene_level_results.loc[tp, "LFC"] = observed_gene_level_lfc
    
    # Sort by timepoint and return as lists
    sorted_results = gene_level_results.sort_index()
    return sorted_results.round(3)


@logger.catch
def calculate_gene_level_statistics(GMs_annotated: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate gene-level statistics for all genes.
    
    Args:
        GMs_annotated: Annotated DataFrame with insertion-level data
        method: P-value combination method
        
    Returns:
        DataFrame with gene-level LFC and p-values across timepoints
    """
    logger.info("Calculating gene-level statistics")
    
    Gene_level_statistics = pd.DataFrame()
    unique_genes = GMs_annotated["Systematic ID"].unique()
    total_genes = len(unique_genes)
    
    start_time = time.time()
    
    for n, (sysID, gene, viability, ess) in enumerate(
        GMs_annotated.groupby(["Systematic ID", "Name", "FYPOviability", "DeletionLibrary_essentiality"]).groups.keys(), 1
    ):
        gene_df = GMs_annotated.query(
            f"`Systematic ID` == '{sysID}' and Name == '{gene}' and FYPOviability == '{viability}' and DeletionLibrary_essentiality == '{ess}'"
        )
        
        try:
            lfc_and_weights = gene_level_lfc_and_weights(gene_df)
            
            Gene_level_statistics.loc[sysID, "Name"] = gene
            Gene_level_statistics.loc[sysID, "FYPOviability"] = viability
            Gene_level_statistics.loc[sysID, "DeletionLibrary_essentiality"] = ess
            
            # Store LFC and p-values for each timepoint
            timepoints = lfc_and_weights.index.tolist()
            for i, tp in enumerate(timepoints):
                if i < len(lfc_and_weights["LFC"]):
                    Gene_level_statistics.loc[sysID, tp] = lfc_and_weights.loc[tp, "LFC"]
                else:
                    Gene_level_statistics.loc[sysID, tp] = np.nan
            
        except Exception as e:
            logger.error(f"Error processing gene {sysID} ({gene}): {e}")
            continue
        
        # Progress reporting
        if n % 100 == 0 or n == total_genes:
            elapsed = time.time() - start_time
            rate = n / elapsed
            eta = (total_genes - n) / rate if rate > 0 else 0
            logger.info(f"Processed {n}/{total_genes} genes "
                        f"({n/total_genes*100:.1f}%) - ETA: {eta:.0f}s")
    
    # Remove genes with all NaN values
    complete_genes = Gene_level_statistics.dropna(how='all', axis=0)
    logger.info(f"Completed analysis for {len(complete_genes)} genes")
    
    return complete_genes


@logger.catch
def generate_summary_statistics(Gene_level_statistics: pd.DataFrame) -> Dict[str, Union[int, float]]:
    """
    Generate summary statistics for gene-level analysis results.
    
    Args:
        Gene_level_statistics: DataFrame with gene-level results
        
    Returns:
        Dictionary with summary statistics
    """
    
    stats = {
        'Total genes analyzed': len(Gene_level_statistics),
        'FYPOviability: Essential genes': len(Gene_level_statistics[Gene_level_statistics['FYPOviability'] == 'inviable']),
        'FYPOviability: Non-essential genes': len(Gene_level_statistics[Gene_level_statistics['FYPOviability'] == 'viable']),
        'DeletionLibrary_essentiality: Essential genes': len(Gene_level_statistics[Gene_level_statistics['DeletionLibrary_essentiality'] == 'E']),
        'DeletionLibrary_essentiality: Non-essential genes': len(Gene_level_statistics[Gene_level_statistics['DeletionLibrary_essentiality'] == 'V'])
    }
    
    return stats

@logger.catch
def display_summary_table(stats: Dict[str, Union[int, float]]) -> None:
    """Display summary statistics in formatted table."""
    logger.info("\n" + "="*60)
    logger.info("GENE-LEVEL DEPLETION ANALYSIS SUMMARY")
    logger.info("="*60)
    
    for key, value in stats.items():
        if isinstance(value, float):
            logger.info(f"{key:<40}: {value:.3f}")
        else:
            logger.info(f"{key:<40}: {value}")
    
    logger.info("="*60)


@logger.catch
def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='''
        Gene-level depletion analysis for transposon insertion sequencing.
        
        This script aggregates insertion-level LFC and p-value data to generate
        gene-level fitness scores using weighted averages and statistical
        meta-analysis methods.
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-l', '--lfc_path', type=Path, required=True,
                       help='Path to CSV file with LFC and p-value results')
    parser.add_argument('-a', '--annotations_path', type=Path, required=True,
                       help='Path to CSV file with insertion annotations')
    parser.add_argument('-w', '--weights_path', type=Path, required=True,
                       help='Path to CSV file with weights')
    parser.add_argument('-o', '--output_path', type=Path, required=True,
                       help='Path for output CSV file with gene-level results')
    
    return parser.parse_args()


def main() -> None:
    """Main execution function."""
    start_time = time.time()
    
    # Parse arguments and setup
    args = parse_arguments()
    setup_logging()
    
    logger.info("Starting gene-level depletion analysis")
    logger.info(f"LFC file: {args.lfc_path}")
    logger.info(f"Weights file: {args.weights_path}")
    logger.info(f"Annotations file: {args.annotations_path}")
    
    try:
        # Create output directory
        args.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Load and process data
        results, weights, insertion_annotations, in_gene_insertions = load_data(
            args.lfc_path, args.weights_path, args.annotations_path
        )
        
        in_gene_LFCs, in_gene_weights = process_lfc_data(results, weights, in_gene_insertions)
        
        # Merge LFC and p-value data
        GMs = pd.merge(
            in_gene_LFCs.stack().to_frame("LFC"), 
            in_gene_weights.stack().to_frame("Weights"), 
            how="left", left_index=True, right_index=True
        )
        
        # Calculate weights and merge with annotations
        GMs = calculate_statistical_weights(GMs)
        GMs_annotated = merge_and_normalize_weights(GMs, insertion_annotations)
        
        # Calculate gene-level statistics
        Gene_level_statistics = calculate_gene_level_statistics(GMs_annotated)
        
        # Generate and display summary statistics
        stats = generate_summary_statistics(Gene_level_statistics)
        display_summary_table(stats)
        
        # Save results
        Gene_level_statistics = Gene_level_statistics.rename_axis("Systematic ID")
        Gene_level_statistics.to_csv(args.output_path.parent/"LFC.tsv", index=True, sep="\t")
        Gene_level_statistics.to_csv(args.output_path, index=True, sep="\t")
        
        # Final summary
        elapsed_time = time.time() - start_time
        logger.info(f"\nAnalysis completed successfully in {elapsed_time:.1f} seconds")
        logger.info(f"Results saved to: {args.output_path}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
