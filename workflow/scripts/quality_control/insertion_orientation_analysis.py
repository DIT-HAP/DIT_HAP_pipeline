"""
Insertion Orientation Analysis Script

This script analyzes strand orientation (+/-) pairs from multiple TSV files with multi-level
indexing. For each file, it creates a single figure with subplots arranged in 1 row × n columns
(where n = number of numeric columns). Each subplot shows all +/- strand pairs for that column
using log-scale scatter plots. Results are saved as a multi-page PDF report.

Typical usage:
  python insertion_orientation_analysis.py -i file1.tsv file2.tsv file3.tsv -o orientation_analysis.pdf
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any, Union
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as mticker
from tabulate import tabulate
import warnings

# Global logger instance
logger = logging.getLogger(__name__)

def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration.

    Args:
        verbose (bool): If True, sets logging level to DEBUG, otherwise INFO.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    # Suppress specific matplotlib warnings
    warnings.filterwarnings('ignore', category=UserWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning)

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: An object holding the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Analyze insertion orientation (+/-) strand pairs from multiple TSV files.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    
    io_group = parser.add_argument_group('Input/Output Arguments')
    io_group.add_argument(
        '-i', '--input',
        nargs='+',
        type=Path,
        required=True,
        help='One or more input TSV files with multi-level indexing.'
    )
    io_group.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output PDF file path for the plots.'
    )
    
    plot_group = parser.add_argument_group('Plotting Parameters')
    plot_group.add_argument(
        '--subplot_width',
        type=float,
        default=4.0,
        help='Width of each subplot in inches (default: %(default)s).'
    )
    plot_group.add_argument(
        '--subplot_height',
        type=float,
        default=4.0,
        help='Height of each subplot in inches (default: %(default)s).'
    )
    plot_group.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Resolution for plots (default: %(default)s).'
    )
    plot_group.add_argument(
        '--point_size',
        type=int,
        default=30,
        help='Size of scatter plot points (default: %(default)s).'
    )
    
    control_group = parser.add_argument_group('Logging and Control')
    control_group.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose (DEBUG level) logging.'
    )
    
    parser.epilog = f"Example: python {parser.prog} -i file1.tsv file2.tsv -o strand_comparison.pdf"
    
    return parser.parse_args()

def read_tsv_with_multiindex(file_path: Path) -> pd.DataFrame:
    """Read TSV file with multi-level indexing [0,1,2,3].

    Args:
        file_path (Path): Path to the input TSV file.

    Returns:
        pd.DataFrame: DataFrame with multi-level index.

    Raises:
        FileNotFoundError: If the input file doesn't exist.
        pd.errors.EmptyDataError: If the file is empty.
        pd.errors.ParserError: If the file cannot be parsed.
    """
    logger.info(f"Reading TSV file: {file_path}")
    
    try:
        df = pd.read_csv(file_path, sep='\t', index_col=[0, 1, 2, 3], engine='python')
        logger.info(f"Successfully read {len(df)} rows with {len(df.columns)} columns")
        logger.debug(f"Index levels: {df.index.names}")
        logger.debug(f"Columns: {list(df.columns)}")
        return df
    except Exception as e:
        logger.error(f"Failed to read {file_path}: {e}")
        raise

def extract_strand_pairs(df: pd.DataFrame) -> List[Tuple[pd.DataFrame, pd.DataFrame, str]]:
    """Extract +/- strand pairs from the DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame with multi-level index where level 2 is strand.

    Returns:
        List[Tuple[pd.DataFrame, pd.DataFrame, str]]: List of tuples containing 
            (positive_strand_data, negative_strand_data, pair_identifier).
    """
    logger.debug("Extracting +/- strand pairs...")
    
    # Get the third level (index 2) which should be strand
    strand_level = 2
    strand_pairs = []
    
    # Group by all levels except strand to find matching pairs
    grouping_levels = [0, 1, 3]  # All levels except strand (level 2)
    
    for group_key, group_data in df.groupby(level=grouping_levels):
        # Check if we have both + and - strands for this group
        strands_available = group_data.index.get_level_values(strand_level).unique()
        
        if '+' in strands_available and '-' in strands_available:
            # Extract data for each strand
            pos_data = group_data.xs('+', level=strand_level)
            neg_data = group_data.xs('-', level=strand_level)
            
            # Create identifier for this pair
            pair_id = f"{group_key}" if isinstance(group_key, str) else "_".join(map(str, group_key))
            
            strand_pairs.append((pos_data, neg_data, pair_id))
            logger.debug(f"Found strand pair: {pair_id}")
    
    logger.debug(f"Found {len(strand_pairs)} +/- strand pairs")
    return strand_pairs

def create_file_comparison_figure(
    df: pd.DataFrame,
    filename: str,
    subplot_width: float,
    subplot_height: float,
    dpi: int,
    point_size: int
) -> Tuple[plt.Figure, Dict[str, Any]]:
    """Create a figure with subplots comparing +/- strand values for all columns.

    Args:
        df (pd.DataFrame): Input DataFrame with multi-level index.
        filename (str): Name of the source file.
        subplot_width (float): Width of each subplot.
        subplot_height (float): Height of each subplot.
        dpi (int): Figure resolution.
        point_size (int): Size of scatter plot points.

    Returns:
        Tuple[plt.Figure, Dict[str, Any]]: The generated figure and statistics.
    """
    # Plotting style parameters from "python-plotting" rule
    plt.rcParams['font.family'] = 'sans-serif'
    primary_color = '#6479cc'  # Medium blue from the palette
    axis_linewidth = 1.2
    grid_alpha = 0.3
    
    # Typography
    title_fontsize = 14
    axis_label_fontsize = 12
    tick_label_fontsize = 10
    annotation_fontsize = 9
    
    # Extract strand pairs
    strand_pairs = extract_strand_pairs(df)
    
    if not strand_pairs:
        logger.warning(f"No +/- strand pairs found in {filename}")
        # Create empty figure
        fig, ax = plt.subplots(figsize=(subplot_width, subplot_height), dpi=dpi)
        ax.text(0.5, 0.5, f'No +/- strand pairs found\nin {filename}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=tick_label_fontsize)
        ax.set_title(filename, fontsize=title_fontsize, fontweight='bold')
        return fig, {"filename": filename, "pairs_found": 0, "columns_analyzed": 0, "correlations": []}
    
    # Get all numeric columns
    numeric_columns = df.select_dtypes(include=[np.number]).columns.tolist()
    
    if not numeric_columns:
        logger.warning(f"No numeric columns found in {filename}")
        # Create empty figure
        fig, ax = plt.subplots(figsize=(subplot_width, subplot_height), dpi=dpi)
        ax.text(0.5, 0.5, f'No numeric columns found\nin {filename}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=tick_label_fontsize)
        ax.set_title(filename, fontsize=title_fontsize, fontweight='bold')
        return fig, {"filename": filename, "pairs_found": len(strand_pairs), "columns_analyzed": 0, "correlations": []}
    
    logger.info(f"Creating figure for {filename} with {len(numeric_columns)} columns: {numeric_columns}")
    
    # Create figure with subplots: 1 row × n columns
    n_cols = len(numeric_columns)
    fig_width = n_cols * subplot_width
    fig_height = subplot_height
    
    fig, axes = plt.subplots(1, n_cols, figsize=(fig_width, fig_height), dpi=dpi)
    if n_cols == 1:
        axes = [axes]  # Ensure axes is always a list
    
    # Set overall title
    fig.suptitle(f'Strand Orientation Analysis: {filename}', 
                 fontsize=title_fontsize + 2, fontweight='bold', y=0.95)
    
    stats = {
        "filename": filename,
        "pairs_found": len(strand_pairs),
        "columns_analyzed": len(numeric_columns),
        "correlations": []
    }
    
    # Process each column
    for col_idx, column in enumerate(numeric_columns):
        ax = axes[col_idx]
        
        # Remove top and right spines
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines['left'].set_linewidth(axis_linewidth)
        ax.spines['bottom'].set_linewidth(axis_linewidth)
        
        # Collect all data points for this column across all strand pairs
        all_pos_values = []
        all_neg_values = []
        pair_correlations = []
        
        for pos_data, neg_data, pair_id in strand_pairs:
            # Get values for the specific column
            pos_values = pos_data.get(column, pd.Series(dtype=float))
            neg_values = neg_data.get(column, pd.Series(dtype=float))
            
            # Align the data by index and remove NaN values
            aligned_data = pd.DataFrame({'pos': pos_values, 'neg': neg_values}).dropna()
            
            if not aligned_data.empty:
                pos_clean = aligned_data['pos']
                neg_clean = aligned_data['neg']
                
                # Filter positive values for log scale
                positive_mask = (pos_clean > 0) & (neg_clean > 0)
                pos_log = pos_clean[positive_mask]
                neg_log = neg_clean[positive_mask]
                
                if len(pos_log) > 0:
                    all_pos_values.extend(pos_log.tolist())
                    all_neg_values.extend(neg_log.tolist())
                    
                    # Calculate correlation for this pair
                    if len(pos_log) > 1:
                        correlation = np.corrcoef(pos_log, neg_log)[0, 1]
                        if not np.isnan(correlation):
                            pair_correlations.append(correlation)
                            stats["correlations"].append({
                                "filename": filename,
                                "column": column,
                                "pair_id": pair_id,
                                "correlation": correlation,
                                "n_points": len(pos_log)
                            })
        
        # Plot all data points for this column
        if all_pos_values and all_neg_values:
            pos_array = np.array(all_pos_values)
            neg_array = np.array(all_neg_values)
            
            # Create scatter plot
            ax.scatter(pos_array, neg_array, s=point_size, color=primary_color, alpha=0.6,
                      edgecolors='white', linewidth=0.5, rasterized=True)
            
            # Set log scale for both axes
            ax.set_xscale('log')
            ax.set_yscale('log')
            
            # Add diagonal line (y=x)
            min_val = min(pos_array.min(), neg_array.min())
            max_val = max(pos_array.max(), neg_array.max())
            ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.8, linewidth=2.0)
            
            # Calculate overall statistics for this column
            if pair_correlations:
                mean_correlation = np.mean(pair_correlations)
                n_pairs = len(pair_correlations)
                n_points = len(all_pos_values)
                
                # Add statistics text box
                stats_text = f"Pairs: {n_pairs}\nPoints: {n_points:,}\nMean r: {mean_correlation:.3f}"
                ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=annotation_fontsize,
                       verticalalignment='top', bbox=dict(boxstyle='round,pad=0.3', 
                       facecolor='white', alpha=0.8, edgecolor='gray'))
        else:
            # No data to plot
            ax.text(0.5, 0.5, 'No positive data\nfor log scale', ha='center', va='center',
                   transform=ax.transAxes, fontsize=tick_label_fontsize)
        
        # Customize subplot
        ax.set_xlabel('Positive Strand (+)', fontsize=axis_label_fontsize, fontweight='semibold')
        if col_idx == 0:  # Only label y-axis on leftmost subplot
            ax.set_ylabel('Negative Strand (-)', fontsize=axis_label_fontsize, fontweight='semibold')
        ax.set_title(column, fontsize=axis_label_fontsize, fontweight='semibold')
        ax.grid(True, linestyle=':', alpha=grid_alpha, linewidth=0.8)
        ax.tick_params(labelsize=tick_label_fontsize)
    
    plt.tight_layout()
    return fig, stats

def analyze_multiple_files(
    input_files: List[Path],
    output_path: Path,
    subplot_width: float,
    subplot_height: float,
    dpi: int,
    point_size: int
) -> Dict[str, Any]:
    """Analyze strand orientations across multiple files and generate PDF report.

    Args:
        input_files (List[Path]): List of input TSV files.
        output_path (Path): Path for output PDF.
        subplot_width (float): Width of each subplot.
        subplot_height (float): Height of each subplot.
        dpi (int): Figure resolution.
        point_size (int): Scatter plot point size.

    Returns:
        Dict[str, Any]: Aggregated analysis statistics.
    """
    logger.info("Starting multi-file strand orientation analysis...")
    
    # Sort files by name for consistent processing order
    sorted_files = sorted(input_files, key=lambda p: p.name)
    
    all_stats = []
    total_figures = 0
    
    # Generate plots and save to PDF
    with PdfPages(output_path) as pdf:
        for file_path in sorted_files:
            filename = file_path.name
            logger.info(f"--- Processing file: {filename} ---")
            
            try:
                # Read the file
                df = read_tsv_with_multiindex(file_path)
                
                # Create figure for this file
                fig, file_stats = create_file_comparison_figure(
                    df, filename, subplot_width, subplot_height, dpi, point_size
                )
                
                # Save figure to PDF
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                total_figures += 1
                
                all_stats.append(file_stats)
                logger.info(f"Generated figure for {filename} with {file_stats['columns_analyzed']} columns")
                
            except Exception as e:
                logger.error(f"Failed to process {filename}: {e}", exc_info=True)
                error_stats = {"filename": filename, "error": str(e), "pairs_found": 0, "columns_analyzed": 0, "correlations": []}
                all_stats.append(error_stats)
    
    # Aggregate statistics
    aggregated_stats = {
        "files_processed": len(sorted_files),
        "figures_generated": total_figures,
        "file_statistics": all_stats,
        "all_correlations": []
    }
    
    # Collect all correlations
    for file_stats in all_stats:
        if "correlations" in file_stats:
            aggregated_stats["all_correlations"].extend(file_stats["correlations"])
    
    logger.info(f"Processed {len(sorted_files)} files, generated {total_figures} figures")
    return aggregated_stats

def save_summary_statistics(stats: Dict[str, Any], output_path: Path) -> None:
    """Save summary statistics to a TSV file.

    Args:
        stats (Dict[str, Any]): Analysis statistics.
        output_path (Path): Output PDF path (used to derive TSV filename).
    """
    if not stats.get("all_correlations"):
        logger.info("No correlation data to save.")
        return
    
    # Create DataFrame from correlations
    correlations_df = pd.DataFrame(stats["all_correlations"])
    
    # Construct output filename
    tsv_filename = output_path.stem + "_strand_correlations.tsv"
    output_tsv_path = output_path.parent / tsv_filename
    
    try:
        correlations_df.to_csv(output_tsv_path, sep='\t', index=False)
        logger.info(f"Correlation statistics saved to: {output_tsv_path}")
    except Exception as e:
        logger.error(f"Failed to save statistics to {output_tsv_path}: {e}")

def display_summary_table(stats: Dict[str, Any]) -> None:
    """Display summary statistics in a formatted table.

    Args:
        stats (Dict[str, Any]): Analysis statistics.
    """
    logger.info("\n--- Multi-File Analysis Summary ---")
    
    # Overall statistics
    summary_data = [
        ["Files processed", f"{stats['files_processed']:,}"],
        ["Figures generated", f"{stats['figures_generated']:,}"]
    ]
    
    try:
        table_str = tabulate(summary_data, headers=["Metric", "Value"], 
                           tablefmt="grid", stralign="left")
        for line in table_str.split('\n'):
            logger.info(line)
    except Exception as e:
        logger.error(f"Could not generate summary table: {e}")
        for metric, value in summary_data:
            logger.info(f"{metric}: {value}")
    
    # Per-file statistics
    logger.info("\n--- Per-File Statistics ---")
    file_table_data = []
    for file_stats in stats["file_statistics"]:
        if "error" in file_stats:
            file_table_data.append([
                file_stats["filename"],
                "ERROR",
                file_stats.get("error", "Unknown error")[:50]
            ])
        else:
            file_table_data.append([
                file_stats["filename"],
                f"{file_stats['pairs_found']:,}",
                f"{file_stats['columns_analyzed']:,}"
            ])
    
    try:
        file_table_str = tabulate(file_table_data, 
                                headers=["File", "Strand Pairs", "Columns Analyzed"], 
                                tablefmt="grid", stralign="left")
        for line in file_table_str.split('\n'):
            logger.info(line)
    except Exception as e:
        logger.error(f"Could not generate per-file table: {e}")
        for row in file_table_data:
            logger.info(f"  {row[0]}: {row[1]} pairs, {row[2]} columns")
    
    # Correlation statistics
    if stats.get("all_correlations"):
        correlations = [c["correlation"] for c in stats["all_correlations"] if not np.isnan(c["correlation"])]
        if correlations:
            logger.info(f"\nOverall Correlation Statistics:")
            logger.info(f"  Total correlations calculated: {len(correlations):,}")
            logger.info(f"  Mean correlation: {np.mean(correlations):.3f}")
            logger.info(f"  Median correlation: {np.median(correlations):.3f}")
            logger.info(f"  Min correlation: {np.min(correlations):.3f}")
            logger.info(f"  Max correlation: {np.max(correlations):.3f}")

def main() -> None:
    """Main execution function for the insertion orientation analysis."""
    args = parse_arguments()
    setup_logging(args.verbose)
    
    start_time = time.time()
    logger.info("Starting Multi-File Insertion Orientation Analysis script...")
    
    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Input files: {[str(f) for f in args.input]}")
    logger.info(f"Output PDF: {args.output}")
    logger.info(f"Subplot size: {args.subplot_width}x{args.subplot_height}")
    
    try:
        # Perform multi-file analysis
        stats = analyze_multiple_files(
            args.input, args.output, 
            args.subplot_width, args.subplot_height,
            args.dpi, args.point_size
        )
        
        # Display and save results
        display_summary_table(stats)
        save_summary_statistics(stats, args.output)
        
        end_time = time.time()
        total_time = end_time - start_time
        
        logger.info(f"--- Analysis complete. PDF saved to {args.output} ---")
        logger.info(f"Total processing time: {total_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}", exc_info=True)
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
