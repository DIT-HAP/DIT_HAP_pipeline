"""
Insertion Orientation Analysis Script for the DIT-HAP project.

This script analyzes strand orientation (+/-) pairs from multiple TSV files with multi-level indexing.
For each file, it creates a single figure with subplots arranged in 1 row × n columns (where n = number 
of numeric columns). Each subplot shows all +/- strand pairs for that column using log-scale scatter plots. 
Results are saved as a multi-page PDF report with correlation statistics.

Typical Usage:
    python insertion_orientation_analysis.py --input file1.tsv file2.tsv --output orientation_analysis.pdf

Input: One or more TSV files with multi-level indexing where level 2 represents strand orientation (+/-)
Output: Multi-page PDF report with strand orientation analysis plots and correlation statistics TSV file
Other information: The script extracts +/- strand pairs and creates log-scale scatter plots with correlation analysis.
"""

# =============================== Imports ===============================
import sys
import argparse
import time
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pydantic import BaseModel, Field, field_validator


# =============================== Constants ===============================
plt.style.use('/data/c/yangyusheng_optimized/DIT_HAP_pipeline/config/DIT_HAP.mplstyle')
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']


# =============================== Configuration & Models ===============================
class InsertionOrientationAnalysisConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths for insertion orientation analysis."""
    input_files: List[Path] = Field(..., description="List of input TSV files with multi-level indexing")
    output_path: Path = Field(..., description="Path for output PDF file")

    @field_validator('input_files')
    def validate_input_files(cls, v):
        if not v:
            raise ValueError("At least one input file must be provided")
        for file_path in v:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
            if not file_path.suffix.lower() in ['.tsv', '.txt']:
                raise ValueError(f"Input file must be a TSV file: {file_path}")
        return v
    
    @field_validator('output_path')
    def validate_output_path(cls, v):
        if not v.suffix.lower() == '.pdf':
            raise ValueError(f"Output file must be a PDF: {v}")
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True

class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    files_processed: int = Field(..., ge=0, description="Number of files processed")
    figures_generated: int = Field(..., ge=0, description="Number of figures generated")
    file_statistics: List[Dict[str, Any]] = Field(default_factory=list, description="Statistics for each file")
    all_correlations: List[Dict[str, Any]] = Field(default_factory=list, description="All correlation data")


# =============================== Setup Logging ===============================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# =============================== Core Functions ===============================
@logger.catch
def read_tsv_with_multiindex(file_path: Path) -> pd.DataFrame:
    """Read TSV file with multi-level indexing [0,1,2,3]."""
    logger.info(f"Reading TSV file: {file_path}")
    
    df = pd.read_csv(file_path, sep='\t', index_col=[0, 1, 2, 3], header=[0,1])
    logger.info(f"Successfully read {len(df)} rows with {len(df.columns)} columns")
    logger.debug(f"Index levels: {df.index.names}")
    logger.debug(f"Columns: {list(df.columns)}")
    return df

@logger.catch
def extract_strand_pairs(df: pd.DataFrame) -> List[Tuple[pd.DataFrame, pd.DataFrame, str]]:
    """Extract +/- strand pairs from the DataFrame."""
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

@logger.catch
def create_file_comparison_figure(
    df: pd.DataFrame,
    filename: str,
) -> Tuple[plt.Figure, Dict[str, Any]]:
    """Create a figure with subplots comparing +/- strand values for all columns."""
    
    # Extract strand pairs
    strand_pairs = extract_strand_pairs(df)
    
    if not strand_pairs:
        logger.warning(f"No +/- strand pairs found in {filename}")
        # Create empty figure
        fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
        ax.text(0.5, 0.5, f'No +/- strand pairs found\nin {filename}', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(filename)
        return fig, {"filename": filename, "pairs_found": 0, "columns_analyzed": 0, "correlations": []}
    
    # Get all numeric columns
    numeric_columns = df.select_dtypes(include=[np.number]).columns.tolist()
    
    if not numeric_columns:
        logger.warning(f"No numeric columns found in {filename}")
        # Create empty figure
        fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
        ax.text(0.5, 0.5, f'No numeric columns found\nin {filename}', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(filename)
        return fig, {"filename": filename, "pairs_found": len(strand_pairs), "columns_analyzed": 0, "correlations": []}
    
    logger.info(f"Creating figure for {filename} with {len(numeric_columns)} columns: {numeric_columns}")
    
    # Create figure with subplots: 1 row × n columns
    n_rows = len(numeric_columns)
    fig_width = AX_WIDTH
    fig_height = n_rows * AX_HEIGHT
    
    fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height))
    if n_rows == 1:
        axes = [axes]  # Ensure axes is always a list
    
    # Set overall title
    fig.suptitle(f'Strand Orientation Analysis: {filename}', y=1.05)
    
    stats = {
        "filename": filename,
        "pairs_found": len(strand_pairs),
        "columns_analyzed": len(numeric_columns),
        "correlations": []
    }
    
    # Process each column
    for col_idx, column in enumerate(numeric_columns):
        ax = axes[col_idx]
        
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
            ax.scatter(pos_array, neg_array, s=10, alpha=0.5, facecolor="none", edgecolor="gray", rasterized=True)
            
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
                ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment='top')
        else:
            # No data to plot
            ax.text(0.5, 0.5, 'No positive data\nfor log scale', ha='center', va='center',
                   transform=ax.transAxes)
        
        # Customize subplot
        ax.set_xlabel('Positive Strand (+)')
        if col_idx == 0:  # Only label y-axis on leftmost subplot
            ax.set_ylabel('Negative Strand (-)')
        ax.set_title(column)
        ax.grid(True)
    
    plt.tight_layout()
    return fig, stats

@logger.catch
def analyze_multiple_files(
    input_files: List[Path],
    output_path: Path,
) -> Dict[str, Any]:
    """Analyze strand orientations across multiple files and generate PDF report."""
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
                    df, filename
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

@logger.catch
def save_summary_statistics(stats: Dict[str, Any], output_path: Path) -> None:
    """Save summary statistics to a TSV file."""
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


# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(description="Analyze insertion orientation (+/-) strand pairs from multiple TSV files.")
    parser.add_argument("-i", "--input", nargs='+', type=Path, required=True, help="One or more input TSV files with multi-level indexing.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output PDF file path for the plots.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point of the script."""
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    try:
        config = InsertionOrientationAnalysisConfig(
            input_files=args.input,
            output_path=args.output
        )

        logger.info("=== Insertion Orientation Analysis ===")
        logger.info(f"Processing {len(config.input_files)} input files...")
        
        # Sort files by name
        sorted_files = sorted(config.input_files, key=lambda x: x.name)
        logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")
        
        start_time = time.time()
        
        # Perform multi-file analysis
        results = analyze_multiple_files(config.input_files, config.output_path)
        
        # Save results
        save_summary_statistics(results, config.output_path)
        
        end_time = time.time()
        total_time = end_time - start_time
        
        logger.success(f"Analysis complete! Output saved to: {config.output_path}")
        logger.info(f"Generated {results['figures_generated']} figures in PDF")
        logger.info(f"Total processing time: {total_time:.2f} seconds")
        
        # Print summary statistics
        logger.info("\n=== Summary Statistics ===")
        logger.info(f"Files processed: {results['files_processed']}")
        logger.info(f"Figures generated: {results['figures_generated']}")
        logger.info(f"Total correlations: {len(results['all_correlations'])}")
        
        if results['all_correlations']:
            all_correlations = [c['correlation'] for c in results['all_correlations']]
            logger.info(f"Mean correlation: {np.mean(all_correlations):.4f}")
            logger.info(f"Correlation range: {np.min(all_correlations):.4f} to {np.max(all_correlations):.4f}")
        
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()