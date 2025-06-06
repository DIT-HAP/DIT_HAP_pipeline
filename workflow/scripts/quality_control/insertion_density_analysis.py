"""
Insertion Density Analysis for Transposon Insertion Sequencing

This script analyzes the density and distribution patterns of transposon insertions
within genes to assess mutagenesis quality and coverage. It calculates comprehensive
statistics including insertion density, gap analysis, read distribution inequality,
and strand preferences.

The analysis includes:
1. Loading insertion LFC data and genomic annotations
2. Filtering for in-gene insertions using established criteria
3. Calculating insertion and site density metrics
4. Analyzing gap regions between insertions
5. Measuring read distribution inequality (concentration vs. dispersion)
6. Calculating strand preference and pairing statistics

Key Metrics:
- Insertion density: Number of insertions per gene length
- Site density: Number of unique sites (Forward/Reverse counted as one)
- Gap analysis: Statistics on regions without insertions
- Gini coefficient: Measure of read distribution inequality
- Strand bias: Preference for Forward vs. Reverse insertions

Typical usage:
    python insertion_density_analysis.py -l LFC.csv -a annotations.tsv -o density_stats.csv

Input: CSV files with insertion LFC data and genomic annotations
Output: CSV file with comprehensive insertion density statistics per gene
"""

import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import argparse
from collections import defaultdict
import time
import warnings

# Suppress pandas warnings for cleaner output
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

# Configure matplotlib for publication quality
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 10,
    'axes.linewidth': 1.2,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.8
})

# Color palette following cursor rules
COLOR_PALETTE = {
    'primary_purple': '#962955',
    'primary_green': '#7fb775', 
    'primary_blue': '#6479cc',
    'primary_gold': '#ad933c',
    'neutral_gray': '#89afba'
}


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def load_insertion_data(insertion_data_path: Path, initial_timepoint: str) -> pd.DataFrame:
    """
    Load insertion data and extract initial timepoint read counts.
    
    Args:
        insertion_data_path: Path to CSV file with the insertion data
        initial_timepoint: Column name for initial timepoint
        
    Returns:
        DataFrame with insertion data indexed by genomic coordinates
        
    Raises:
        FileNotFoundError: If insertion data file doesn't exist
        ValueError: If required columns are missing
    """
    logging.info(f"Loading insertion data from {insertion_data_path}")
    
    if not insertion_data_path.exists():
        raise FileNotFoundError(f"Insertion data file not found: {insertion_data_path}")
    
    try:
        # Load with multi-level index
        insertion_data = pd.read_csv(insertion_data_path, index_col=[0, 1, 2, 3], header=[0, 1])
        
        # Validate that initial timepoint exists
        if initial_timepoint not in insertion_data.columns:
            available_cols = insertion_data.columns.tolist()
            raise ValueError(
                f"Initial timepoint '{initial_timepoint}' not found. "
                f"Available columns: {available_cols}"
            )
        
        # Extract initial timepoint data (used as read counts proxy)
        insertion_data = insertion_data[initial_timepoint]
        
        logging.info(f"Loaded {len(insertion_data)} insertions with read count data")
        
        return insertion_data
        
    except Exception as e:
        raise ValueError(f"Error loading insertion data: {e}")


def load_annotation_data(annotations_path: Path) -> pd.DataFrame:
    """
    Load genomic annotations for insertions.
    
    Args:
        annotations_path: Path to TSV file with annotations
        
    Returns:
        DataFrame with annotations indexed by genomic coordinates
        
    Raises:
        FileNotFoundError: If annotations file doesn't exist
        ValueError: If required columns are missing
    """
    logging.info(f"Loading annotation data from {annotations_path}")
    
    if not annotations_path.exists():
        raise FileNotFoundError(f"Annotations file not found: {annotations_path}")
    
    try:
        # Load with multi-level index and tab separator
        annotations = pd.read_csv(annotations_path, index_col=[0, 1, 2, 3], sep="\t")
        
        # Validate required columns exist
        required_cols = ['Type', 'Distance_to_stop_codon', 'Systematic ID', 'Name', 'Chr_Interval', 'Strand_Interval', 'ParentalRegion_start', 'ParentalRegion_end', 'ParentalRegion_length', 'Insertion_direction']
        missing_cols = [col for col in required_cols if col not in annotations.columns]
        if missing_cols:
            raise ValueError(f"Missing required annotation columns: {missing_cols}")
        
        logging.info(f"Loaded annotations for {len(annotations)} insertions")
        logging.info(f"Annotation columns: {annotations.columns.tolist()}")
        
        return annotations
        
    except Exception as e:
        raise ValueError(f"Error loading annotation data: {e}")


def filter_in_gene_insertions(insertion_data: pd.DataFrame, 
                             annotations: pd.DataFrame) -> pd.DataFrame:
    """
    Filter insertions to include only those within genes using established criteria.
    
    Args:
        insertion_data: DataFrame with insertion read counts
        annotations: DataFrame with genomic annotations
        
    Returns:
        DataFrame with in-gene insertions and their annotations
    """
    logging.info("Filtering for in-gene insertions")
    
    # Merge insertion data with annotations
    merged_data = pd.merge(
        insertion_data, annotations,
        left_index=True, right_index=True,
        how='inner'
    )
    
    # Apply in-gene filtering criteria (same as gene_level_depletion_analysis.py)
    in_gene_mask = (
        (merged_data['Type'] != 'Intergenic region') & 
        (merged_data['Distance_to_stop_codon'] > 4)
    )
    
    in_gene_insertions = merged_data[in_gene_mask].copy()
    
    logging.info(f"Found {len(in_gene_insertions)} in-gene insertions")
    logging.info(f"Filtered out {len(merged_data) - len(in_gene_insertions)} intergenic/near-stop insertions")
    
    return in_gene_insertions

def calculate_insertion_statistics(gene_insertions: pd.DataFrame) -> Dict[str, Union[int, float]]:
    """
    Calculate basic insertion and site statistics for a gene.
    
    Args:
        gene_insertions: DataFrame with insertions for a single gene
        
    Returns:
        Dictionary with insertion statistics
    """
    # Basic counts
    total_insertions = len(gene_insertions)
    
    # Count unique sites (same coordinate regardless of strand)
    coordinates = gene_insertions.index.get_level_values(1)
    unique_sites = len(coordinates.unique())
    
    # Calculate normalized densities (per 1000 bp)
    gene_length = gene_insertions["ParentalRegion_length"].iloc[0]
    insertion_density = (total_insertions / gene_length) * 1000 if gene_length > 0 else 0
    site_density = (unique_sites / gene_length) * 1000 if gene_length > 0 else 0
    
    return {
        'total_insertions': total_insertions,
        'unique_sites': unique_sites,
        'gene_length': gene_length,
        'insertion_density_per_kb': round(insertion_density, 3),
        'site_density_per_kb': round(site_density, 3)
    }


def calculate_gap_statistics(gene_insertions: pd.DataFrame) -> Dict[str, Union[int, float, str]]:
    """
    Calculate statistics about gaps between insertions within a gene.
    
    Args:
        gene_insertions: DataFrame with insertions for a single gene
        
    Returns:
        Dictionary with gap statistics
    """
    coordinates = sorted(gene_insertions.index.get_level_values(1).unique())
    start_coordinate = gene_insertions["ParentalRegion_start"].iloc[0]
    end_coordinate = gene_insertions["ParentalRegion_end"].iloc[0]
    gene_length = gene_insertions["ParentalRegion_length"].iloc[0]

    coordinates_with_start_and_end = sorted(set([start_coordinate] + coordinates + [end_coordinate]))
    
    # Calculate gaps between consecutive insertions
    gaps = [coordinates_with_start_and_end[i+1] - coordinates_with_start_and_end[i] - 1 for i in range(len(coordinates_with_start_and_end)-1)]
    gaps = [gap for gap in gaps if gap > 0]  # Only count actual gaps

    normalized_gaps = [round(gap / gene_length, 3) for gap in gaps]
    gini_coefficient_of_location = calculate_gini_coefficient(normalized_gaps)
    
    if not gaps:
        return {
            'num_gaps': 0,
            'largest_gap': 0,
            'largest_gap_fraction': 0,
            'smallest_gap': 0,
            'smallest_gap_fraction': 0,
            'mean_gap_length': 0,
            'mean_gap_length_fraction': 0,
            'median_gap_length': 0,
            'median_gap_length_fraction': 0,
            'gap_length_sd': 0,
            'gap_length_sd_fraction': 0,
            'all_gap_lengths': "",
            'all_gap_lengths_fraction': "",
            'gini_coefficient_of_location': np.nan
        }
    
    return {
        'num_gaps': len(gaps),
        'largest_gap': max(gaps),
        'largest_gap_fraction': max(normalized_gaps),
        'smallest_gap': min(gaps),
        'smallest_gap_fraction': min(normalized_gaps),
        'mean_gap_length': round(np.mean(gaps), 2),
        'mean_gap_length_fraction': round(np.mean(normalized_gaps), 2),
        'median_gap_length': round(np.median(gaps), 2),
        'median_gap_length_fraction': round(np.median(normalized_gaps), 2),
        'gap_length_sd': round(np.std(gaps), 2),
        'gap_length_sd_fraction': round(np.std(normalized_gaps), 2),
        'all_gap_lengths': ",".join(map(str, gaps)),
        'all_gap_lengths_fraction': ",".join(map(str, normalized_gaps)),
        'gini_coefficient_of_location': round(gini_coefficient_of_location, 3)
    }


def calculate_gini_coefficient(values: np.ndarray) -> float:
    """
    Calculate Gini coefficient to measure inequality in read distribution.
    
    The Gini coefficient ranges from 0 (perfect equality) to 1 (maximum inequality).
    Higher values indicate more concentrated/uneven distribution of reads.
    
    Args:
        values: Array of read counts
        
    Returns:
        Gini coefficient (0-1)
    """
    if len(values) == 0:
        return 0.0
    
    # Sort values
    sorted_values = np.sort(values)
    n = len(sorted_values)
    
    # Calculate Gini coefficient
    cumsum = np.cumsum(sorted_values)
    gini = (2 * np.sum((np.arange(1, n+1) * sorted_values))) / (n * cumsum[-1]) - (n + 1) / n
    
    return max(0.0, min(1.0, gini))  # Ensure result is between 0 and 1


def calculate_read_statistics(gene_insertions: pd.DataFrame) -> Dict[str, Union[int, float]]:
    """
    Calculate read distribution statistics for insertions within a gene.
    
    Args:
        gene_insertions: DataFrame with insertions for a single gene
        
    Returns:
        Dictionary with read statistics
    """
    read_counts = gene_insertions['baseMean'].values
    
    if len(read_counts) == 0:
        return {
            'total_reads': 0,
            'mean_reads_per_insertion': 0,
            'median_reads_per_insertion': 0,
            'read_count_sd': 0,
            'gini_coefficient_of_depth': 0
        }
    
    # Calculate basic statistics
    total_reads = read_counts.sum()
    mean_reads = np.mean(read_counts)
    median_reads = np.median(read_counts)
    read_sd = np.std(read_counts)
    
    # Calculate Gini coefficient for read distribution inequality
    gini_coeff_of_depth = calculate_gini_coefficient(read_counts)
    
    return {
        'total_reads': int(total_reads),
        'mean_reads_per_insertion': round(mean_reads, 2),
        'median_reads_per_insertion': round(median_reads, 2),
        'read_count_sd': round(read_sd, 2),
        'gini_coefficient_of_depth': round(gini_coeff_of_depth, 3)
    }


def calculate_strand_statistics(gene_insertions: pd.DataFrame) -> Dict[str, Union[int, float]]:
    """
    Calculate strand preference and pairing statistics.
    
    Args:
        gene_insertions: DataFrame with insertions for a single gene
        
    Returns:
        Dictionary with strand statistics
    """
    strands = gene_insertions["Insertion_direction"].values
    coordinates = gene_insertions.index.get_level_values(1)
    
    # Count forward and reverse insertions
    forward_count = (strands == 'Forward').sum()
    reverse_count = (strands == 'Reverse').sum()
    total_insertions = len(strands)
    total_sites = len(coordinates.unique())
    
    # Calculate strand preference
    forward_preference = forward_count / total_insertions if total_insertions > 0 else 0
    reverse_preference = reverse_count / total_insertions if total_insertions > 0 else 0
    
    # Calculate strand bias (absolute difference from 50:50)
    strand_bias = abs(forward_preference - 0.5)
    
    # Count paired insertions (same coordinate, different strands)
    coord_strand_pairs = list(zip(coordinates, strands))
    coord_counts = defaultdict(lambda: {'forward': 0, 'reverse': 0})
    
    for coord, strand in coord_strand_pairs:
        if strand == 'Forward':
            coord_counts[coord]['forward'] += 1
        else:
            coord_counts[coord]['reverse'] += 1
    
    # Count sites with both forward and reverse insertions
    paired_sites = sum(1 for counts in coord_counts.values() 
                      if counts['forward'] > 0 and counts['reverse'] > 0)
    
    paired_sites_fraction = paired_sites / total_sites if total_sites > 0 else 0
    
    return {
        'forward_insertions': forward_count,
        'reverse_insertions': reverse_count,
        'forward_preference': round(forward_preference, 3),
        'reverse_preference': round(reverse_preference, 3),
        'strand_bias': round(strand_bias, 3),
        'paired_sites': paired_sites,
        'paired_sites_fraction': round(paired_sites_fraction, 3)
    }


def analyze_gene_insertions(gene_id: str, gene_insertions: pd.DataFrame) -> Dict[str, Union[str, int, float]]:
    """
    Perform comprehensive analysis of insertions within a single gene.
    
    Args:
        gene_id: Systematic gene ID
        gene_insertions: DataFrame with insertions for this gene
        
    Returns:
        Dictionary with all calculated statistics for the gene
    """
    # Calculate all statistics
    insertion_stats = calculate_insertion_statistics(gene_insertions)
    gap_stats = calculate_gap_statistics(gene_insertions)
    read_stats = calculate_read_statistics(gene_insertions)
    strand_stats = calculate_strand_statistics(gene_insertions)
    
    # Combine all statistics
    gene_analysis = {
        'Systematic ID': gene_id,
        'Name': gene_insertions['Name'].iloc[0],
        'Chr': gene_insertions['Chr_Interval'].iloc[0],
        'Start': gene_insertions['ParentalRegion_start'].iloc[0],
        'End': gene_insertions['ParentalRegion_end'].iloc[0],
        'Length': gene_insertions['ParentalRegion_length'].iloc[0],
        'Strand': gene_insertions['Strand_Interval'].iloc[0],
    }
    
    gene_analysis.update(insertion_stats)
    gene_analysis.update(gap_stats)
    gene_analysis.update(read_stats)
    gene_analysis.update(strand_stats)
    
    return gene_analysis


def generate_summary_statistics(results_df: pd.DataFrame) -> Dict[str, Union[int, float]]:
    """
    Generate summary statistics across all analyzed genes.
    
    Args:
        results_df: DataFrame with gene-level analysis results
        
    Returns:
        Dictionary with summary statistics
    """
    stats = {
        'total_genes_analyzed': len(results_df),
        'total_insertions_analyzed': results_df['total_insertions'].sum(),
        'total_unique_sites': results_df['unique_sites'].sum(),
        'mean_insertions_per_gene': results_df['total_insertions'].mean(),
        'median_insertions_per_gene': results_df['total_insertions'].median(),
        'mean_insertion_density_per_kb': results_df['insertion_density_per_kb'].mean(),
        'median_insertion_density_per_kb': results_df['insertion_density_per_kb'].median(),
        'mean_gini_coefficient_of_location': results_df['gini_coefficient_of_location'].mean(),
        'genes_with_high_inequality_of_location': len(results_df[results_df['gini_coefficient_of_location'] > 0.5]),
        'mean_gini_coefficient_of_depth': results_df['gini_coefficient_of_depth'].mean(),
        'genes_with_high_inequality_of_depth': len(results_df[results_df['gini_coefficient_of_depth'] > 0.5]),
        'mean_strand_bias': results_df['strand_bias'].mean(),
        'genes_with_strong_strand_bias': len(results_df[results_df['strand_bias'] > 0.2]),
        'mean_paired_sites_fraction': results_df['paired_sites_fraction'].mean(),
        'genes_with_high_paired_sites_fraction': len(results_df[results_df['paired_sites_fraction'] > 0.5])
    }
    
    return stats


def plot_numeric_distributions_to_pdf(results_df: pd.DataFrame, output_path: Path) -> None:
    """
    Generate histograms for all numeric columns and save to a multi-page PDF.
    
    Args:
        results_df: DataFrame with gene-level analysis results
        output_path: Path for the output PDF file
    """
    logging.info("Generating histograms for numeric columns in PDF format")
    
    # Identify numeric columns (excluding string columns like gene names)
    numeric_columns = results_df.select_dtypes(include=[np.number]).columns.tolist()
    
    # Remove columns that are identifiers or coordinates
    exclude_patterns = ['Start', 'End', 'Length', 'Chr']
    numeric_columns = [col for col in numeric_columns 
                      if not any(pattern in col for pattern in exclude_patterns)]
    
    if not numeric_columns:
        logging.warning("No numeric columns found for plotting")
        return
    
    # Group columns by category for better organization
    column_groups = {
        'Density Metrics': [col for col in numeric_columns if 'density' in col.lower()],
        'Gap Statistics': [col for col in numeric_columns if 'gap' in col.lower()],
        'Read Statistics': [col for col in numeric_columns if any(term in col.lower() 
                           for term in ['read', 'gini_coefficient_of_depth'])],
        'Strand Statistics': [col for col in numeric_columns if any(term in col.lower() 
                             for term in ['forward', 'reverse', 'strand', 'paired'])],
        'Location Statistics': [col for col in numeric_columns if 'gini_coefficient_of_location' in col.lower()],
        'Count Statistics': [col for col in numeric_columns if any(term in col.lower() 
                            for term in ['total_insertions', 'unique_sites', 'num_gaps'])]
    }
    
    # Add any remaining columns to a general category
    all_categorized = set()
    for group_cols in column_groups.values():
        all_categorized.update(group_cols)
    
    remaining_cols = [col for col in numeric_columns if col not in all_categorized]
    if remaining_cols:
        column_groups['Other Statistics'] = remaining_cols
    
    # Create PDF file
    pdf_path = output_path.parent / f"{output_path.stem}_histograms.pdf"
    
    with PdfPages(pdf_path) as pdf:
        # Create title page
        create_title_page(pdf, results_df)
        
        # Create summary plot with key metrics first
        create_summary_histogram_plot_pdf(results_df, pdf)
        
        # Generate plots for each group
        for group_name, group_columns in column_groups.items():
            if not group_columns:
                continue
                
            logging.info(f"Plotting {group_name}: {len(group_columns)} columns")
            
            # Calculate subplot layout
            n_cols = min(3, len(group_columns))
            n_rows = int(np.ceil(len(group_columns) / n_cols))
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
            if n_rows == 1 and n_cols == 1:
                axes = [axes]
            elif n_rows == 1 or n_cols == 1:
                axes = axes.flatten()
            else:
                axes = axes.flatten()
            
            for idx, column in enumerate(group_columns):
                ax = axes[idx] if len(group_columns) > 1 else axes[0]
                
                # Get data and remove NaN values
                data = results_df[column].dropna()
                
                if len(data) == 0:
                    ax.text(0.5, 0.5, 'No Data', transform=ax.transAxes, 
                           ha='center', va='center', fontsize=12)
                    ax.set_title(column, fontsize=11, fontweight='semibold')
                    continue
                
                # Choose color based on column type
                if 'density' in column.lower():
                    color = COLOR_PALETTE['primary_blue']
                elif 'gap' in column.lower():
                    color = COLOR_PALETTE['primary_green']
                elif any(term in column.lower() for term in ['read', 'gini']):
                    color = COLOR_PALETTE['primary_purple']
                else:
                    color = COLOR_PALETTE['primary_gold']
                
                # Check if this is a read depth related column that should be log transformed
                is_read_depth = any(term in column.lower() for term in 
                                  ['read', 'basemean', 'count', 'depth'])
                
                # Transform data if it's read depth related
                if is_read_depth and data.min() > 0:
                    # Log transform the data values (add small constant to handle zeros)
                    plot_data = np.log10(data + 1)
                    xlabel = 'log10(Value + 1)'
                    
                    # Calculate statistics on original data for annotation
                    mean_val = data.mean()
                    median_val = data.median()
                    
                    # Calculate mean and median of transformed data for reference lines
                    plot_mean = plot_data.mean()
                    plot_median = plot_data.median()
                else:
                    # Use original data for other metrics
                    plot_data = data
                    xlabel = 'Value'
                    mean_val = data.mean()
                    median_val = data.median()
                    plot_mean = mean_val
                    plot_median = median_val
                
                # Create histogram
                n_bins = min(30, max(10, int(np.sqrt(len(plot_data)))))
                ax.hist(plot_data, bins=n_bins, color=color, alpha=0.7, 
                       edgecolor='white', linewidth=0.5)
                
                # Add statistics text (always show original data statistics)
                stats_text = f'Mean: {mean_val:.3f}\nMedian: {median_val:.3f}\nN: {len(data)}'
                
                ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                       fontsize=8, verticalalignment='top',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
                
                # Formatting
                ax.set_title(column.replace('_', ' ').title(), fontsize=11, fontweight='semibold')
                ax.set_xlabel(xlabel, fontsize=10)
                ax.set_ylabel('Frequency', fontsize=10)
                ax.grid(True, alpha=0.3)
                
                # Add vertical lines for mean and median (on transformed scale if applicable)
                ax.axvline(plot_mean, color='red', linestyle='--', alpha=0.7, linewidth=1)
                ax.axvline(plot_median, color='orange', linestyle=':', alpha=0.7, linewidth=1)
            
            # Hide unused subplots
            for idx in range(len(group_columns), len(axes)):
                axes[idx].set_visible(False)
            
            # Add group title
            fig.suptitle(f'{group_name} - Distribution Analysis', 
                        fontsize=16, fontweight='bold', y=0.98)
            plt.tight_layout()
            
            # Save page to PDF
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    
    logging.info(f"Multi-page histogram PDF saved to {pdf_path}")


def create_title_page(pdf: PdfPages, results_df: pd.DataFrame) -> None:
    """
    Create a title page for the PDF with analysis summary.
    
    Args:
        pdf: PdfPages object to save the page to
        results_df: DataFrame with analysis results for summary stats
    """
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    
    # Title
    ax.text(0.5, 0.9, 'Insertion Density Analysis', 
           transform=ax.transAxes, fontsize=24, fontweight='bold', 
           ha='center', va='center')
    
    ax.text(0.5, 0.85, 'Histogram Distribution Report', 
           transform=ax.transAxes, fontsize=16, 
           ha='center', va='center', style='italic')
    
    # Analysis summary
    summary_text = f"""
Analysis Summary:
• Total genes analyzed: {len(results_df):,}
• Numeric metrics calculated: {len(results_df.select_dtypes(include=[np.number]).columns)}
• Analysis includes: insertion density, gap statistics, read distribution, and strand preferences

Report Contents:
1. Key Metrics Summary (6 most important measures)
2. Density Metrics (insertion and site density per kb)
3. Gap Statistics (gap lengths, counts, and distributions)
4. Read Statistics (read counts and depth inequality)
5. Strand Statistics (forward/reverse preferences and pairing)
6. Location Statistics (spatial inequality measures)
7. Count Statistics (total insertions, unique sites, gap counts)

Statistical Annotations:
• Red dashed line: Mean value
• Orange dotted line: Median value
• Text box: Mean, Median, and sample size (N)
• Bin count optimized using square root rule
• Read depth metrics: log10 transformed for better visualization
• Statistics shown: Original data values (not transformed)
"""
    
    ax.text(0.1, 0.7, summary_text, transform=ax.transAxes, fontsize=12, 
           va='top', ha='left', linespacing=1.5)
    
    # Footer
    from datetime import datetime
    ax.text(0.5, 0.1, f'Generated on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', 
           transform=ax.transAxes, fontsize=10, ha='center', va='center',
           style='italic', color='gray')
    
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def create_summary_histogram_plot_pdf(results_df: pd.DataFrame, pdf: PdfPages) -> None:
    """
    Create a summary plot with the most important metrics for PDF.
    
    Args:
        results_df: DataFrame with gene-level analysis results
        pdf: PdfPages object to save the plot to
    """
    # Select key metrics for summary plot
    key_metrics = [
        'insertion_density_per_kb',
        'site_density_per_kb', 
        'gini_coefficient_of_depth',
        'gini_coefficient_of_location',
        'strand_bias',
        'paired_sites_fraction'
    ]
    
    # Filter to only include columns that exist in the data
    available_metrics = [col for col in key_metrics if col in results_df.columns]
    
    if not available_metrics:
        logging.warning("No key metrics available for summary plot")
        return
    
    # Create summary plot
    n_cols = 3
    n_rows = int(np.ceil(len(available_metrics) / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
    if n_rows == 1:
        axes = axes.flatten() if n_cols > 1 else [axes]
    else:
        axes = axes.flatten()
    
    colors = [COLOR_PALETTE['primary_blue'], COLOR_PALETTE['primary_green'], 
              COLOR_PALETTE['primary_purple'], COLOR_PALETTE['primary_gold']] * 2
    
    for idx, metric in enumerate(available_metrics):
        ax = axes[idx]
        data = results_df[metric].dropna()
        
        if len(data) == 0:
            ax.text(0.5, 0.5, 'No Data', transform=ax.transAxes, 
                   ha='center', va='center', fontsize=12)
            ax.set_title(metric.replace('_', ' ').title(), fontsize=12, fontweight='bold')
            continue
        
        # Check if this is a read depth related metric that should be log transformed
        is_read_depth = any(term in metric.lower() for term in 
                          ['read', 'basemean', 'count', 'depth', 'gini_coefficient_of_depth'])
        
        # Transform data if it's read depth related
        if is_read_depth and data.min() > 0:
            # Log transform the data values (add small constant to handle zeros)
            plot_data = np.log10(data + 1)
            xlabel = 'log10(Value + 1)'
            
            # Calculate statistics on original data for annotation
            mean_val = data.mean()
            median_val = data.median()
            std_val = data.std()
            
            # Calculate mean and median of transformed data for reference lines
            plot_mean = plot_data.mean()
            plot_median = plot_data.median()
        else:
            # Use original data for other metrics
            plot_data = data
            xlabel = 'Value'
            mean_val = data.mean()
            median_val = data.median()
            std_val = data.std()
            plot_mean = mean_val
            plot_median = median_val
        
        # Create histogram
        n_bins = min(25, max(10, int(np.sqrt(len(plot_data)))))
        ax.hist(plot_data, bins=n_bins, color=colors[idx], alpha=0.7, 
               edgecolor='white', linewidth=0.5)
        
        # Add statistics (always show original data statistics)
        stats_text = f'Mean: {mean_val:.3f}\nMedian: {median_val:.3f}\nStd: {std_val:.3f}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9))
        
        # Formatting
        title = metric.replace('_', ' ').replace('per kb', '(per kb)').title()
        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel('Number of Genes', fontsize=11)
        ax.grid(True, alpha=0.3)
        
        # Add reference lines (on transformed scale if applicable)
        ax.axvline(plot_mean, color='red', linestyle='--', alpha=0.8, linewidth=1.5, label='Mean')
        ax.axvline(plot_median, color='orange', linestyle=':', alpha=0.8, linewidth=1.5, label='Median')
        
        # Add legend for first plot
        if idx == 0:
            ax.legend(fontsize=9, frameon=False)
    
    # Hide unused subplots
    for idx in range(len(available_metrics), len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle('Key Insertion Density Metrics Distribution', 
                fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    # Save to PDF
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def display_summary_table(stats: Dict[str, Union[int, float]]) -> None:
    """Display summary statistics in formatted table."""
    logging.info("\n" + "="*60)
    logging.info("INSERTION DENSITY ANALYSIS SUMMARY")
    logging.info("="*60)
    
    for key, value in stats.items():
        if isinstance(value, float):
            logging.info(f"{key:<35}: {value:.3f}")
        else:
            logging.info(f"{key:<35}: {value}")
    
    logging.info("="*60)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='''
        Insertion density analysis for transposon insertion sequencing.
        
        This script analyzes the density and distribution patterns of insertions
        within genes, calculating comprehensive statistics including gap analysis,
        read distribution inequality, and strand preferences.
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--insertion_data_path', type=Path, required=True,
                       help='Path to CSV file with insertion data')
    parser.add_argument('-a', '--annotations_path', type=Path, required=True,
                       help='Path to TSV file with insertion annotations')
    parser.add_argument('-o', '--output_path', type=Path, required=True,
                       help='Path for output CSV file with density statistics')
    parser.add_argument('-t', '--initial_timepoint', type=str, required=True,
                       help='Initial timepoint column name for read counts')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose logging')
    
    return parser.parse_args()


def main() -> None:
    """Main execution function."""
    start_time = time.time()
    
    # Parse arguments and setup
    args = parse_arguments()
    setup_logging(args.verbose)
    
    logging.info("Starting insertion density analysis")
    logging.info(f"Insertion data file: {args.insertion_data_path}")
    logging.info(f"Annotations file: {args.annotations_path}")
    logging.info(f"Initial timepoint: {args.initial_timepoint}")
    
    try:
        # Create output directory
        args.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Load data
        insertion_data = load_insertion_data(args.insertion_data_path, args.initial_timepoint)
        annotations = load_annotation_data(args.annotations_path)
        
        # Filter for in-gene insertions
        in_gene_insertions = filter_in_gene_insertions(insertion_data, annotations)

        # valid genes
        valid_genes = in_gene_insertions["Systematic ID"].unique().tolist()
        logging.info(f"Analyzing {len(valid_genes)} genes")
        
        # Analyze each gene
        gene_results = []
        
        for gene_id in valid_genes:
            gene_insertions = in_gene_insertions[
                in_gene_insertions['Systematic ID'] == gene_id
            ]
            
            if len(gene_insertions) > 0:
                gene_analysis = analyze_gene_insertions(gene_id, gene_insertions)
                gene_results.append(gene_analysis)
        
        # Create results DataFrame
        results_df = pd.DataFrame(gene_results)
        results_df = results_df.set_index('Systematic ID')
        
        # Generate and display summary statistics
        stats = generate_summary_statistics(results_df)
        display_summary_table(stats)
        
        # Generate histogram plots in PDF format
        plot_numeric_distributions_to_pdf(results_df, args.output_path)
        
        # Save results
        results_df.to_csv(args.output_path, index=True)
        
        # Final summary
        elapsed_time = time.time() - start_time
        logging.info(f"\nAnalysis completed successfully in {elapsed_time:.1f} seconds")
        logging.info(f"Results saved to: {args.output_path}")
        logging.info(f"Histogram PDF saved to: {args.output_path.parent / f'{args.output_path.stem}_histograms.pdf'}")
        logging.info(f"Analyzed {len(results_df)} genes with insertion data")
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
