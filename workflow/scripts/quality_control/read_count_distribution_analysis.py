"""
Read Count Distribution Analysis Script

This script takes multiple TSV files as input, visualizes the distribution of
values in each column using log-transformed data with original scale labels,
applies a cutoff to a specified initial time point column, and calculates
data retention. Results are saved as a multi-page PDF, and a summary
table of statistics is printed to the console.

Typical usage:
  python read_count_distribution_analysis.py -i sample*.tsv -t T0 -c 10 -o report.pdf

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
    # Suppress specific matplotlib warnings if necessary, and other general warnings
    warnings.filterwarnings('ignore', category=UserWarning) # Matplotlib may warn about fixed formatter
    warnings.filterwarnings('ignore', category=RuntimeWarning) # For log(0) issues if any slip through

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: An object holding the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Analyze read count distributions from TSV files and apply cutoffs.",
        formatter_class=argparse.RawTextHelpFormatter,
        # epilog is dynamically set by default based on usage in recent argparse versions
    )
    
    io_group = parser.add_argument_group('Input/Output Arguments')
    io_group.add_argument(
        '-i', '--input',
        nargs='+',
        type=Path,
        required=True,
        help='One or more input TSV files.'
    )
    io_group.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output PDF file path for the plots.'
    )
    
    analysis_group = parser.add_argument_group('Analysis Parameters')
    analysis_group.add_argument(
        '-t', '--initial_time_point',
        required=True,
        type=str,
        help='Name of the column representing the initial time point for cutoff application.'
    )
    analysis_group.add_argument(
        '-c', '--cutoff',
        required=True,
        type=float,
        help='Cutoff value to apply to the initial time point column (values >= cutoff are kept).'
    )
    
    plot_group = parser.add_argument_group('Plotting Parameters')
    plot_group.add_argument(
        '--bins',
        type=int,
        default=50,
        help='Number of bins for histograms (default: %(default)s).'
    )
    plot_group.add_argument(
        '--figsize_per_subplot_width',
        type=float,
        default=4.5,
        help='Width of each subplot in inches (default: %(default)s).'
    )
    plot_group.add_argument(
        '--figsize_subplot_height',
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
    
    control_group = parser.add_argument_group('Logging and Control')
    control_group.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose (DEBUG level) logging.'
    )
    
    # Update epilog with an example
    parser.epilog = f"Example: python {parser.prog} -i sampleA.tsv sampleB.tsv -t T0_reads -c 100 -o analysis.pdf"
    
    return parser.parse_args()

def plot_distributions_and_calculate_stats(
    df: pd.DataFrame,
    filename: str,
    initial_time_point_col: str,
    cutoff_val: float,
    bins: int,
    figsize_per_subplot_width: float,
    figsize_subplot_height: float,
    dpi: int
) -> Tuple[plt.Figure | None, Dict[str, Any]]:
    """
    Plots distributions for numeric columns and calculates filtering statistics.

    For each numeric column in the DataFrame, a histogram of its log10-transformed
    values is plotted. The x-axis tick labels show values in the original
    scale. A cutoff line is drawn on the plot for the initial time point column.
    Key statistics for the current file are displayed on the plot.

    Args:
        df (pd.DataFrame): DataFrame containing the data.
        filename (str): Name of the input file, used for titles and logging.
        initial_time_point_col (str): Column name for applying the cutoff.
        cutoff_val (float): Cutoff value for filtering.
        bins (int): Number of bins for histograms.
        figsize_per_subplot_width (float): Width for each subplot.
        figsize_subplot_height (float): Height for each subplot.
        dpi (int): Dots per inch for the figure.

    Returns:
        Tuple[plt.Figure | None, Dict[str, Any]]:
            - plt.Figure | None: The generated matplotlib Figure, or None if no
              numeric data is available for plotting.
            - Dict[str, Any]: A dictionary containing statistics about data
              retention after applying the cutoff.
    """
    numeric_cols = df.select_dtypes(include=np.number).columns
    if not numeric_cols.any():
        logger.warning(f"No numeric columns found in {filename}. Skipping plot generation.")
        return None, {}

    # Plotting style parameters from "python-plotting" rule
    plt.rcParams['font.family'] = 'sans-serif' # Rule: Font family
    # plt.rcParams['font.sans-serif'] = ['Arial'] # Rule: Arial preferred - often needs setup

    primary_hist_color = '#7fb775' # Rule: Green Family - Medium green
    axis_linewidth = 1.2
    data_linewidth = 2.0 # For cutoff line
    grid_alpha = 0.3
    
    # Typography
    title_fontsize = 14
    axis_label_fontsize = 12
    tick_label_fontsize = 10
    legend_fontsize = 10
    annotation_fontsize = 9


    num_subplots = len(numeric_cols)
    fig, axes = plt.subplots(
        1, num_subplots,
        figsize=(num_subplots * figsize_per_subplot_width, figsize_subplot_height),
        sharey=True,
        dpi=dpi
    )
    if num_subplots == 1:
        axes = [axes]

    fig.suptitle(
        f"Value Distributions: {filename}\\nInitial Time Point: '{initial_time_point_col}', Cutoff >= {cutoff_val}",
        fontsize=title_fontsize, fontweight='bold', y=0.98 
    )

    max_y_val = 0

    for i, col_name in enumerate(numeric_cols):
        ax = axes[i]
        col_data = df[col_name].dropna()
        positive_col_data = col_data[col_data > 0]

        # Remove top and right spines - Rule
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines['left'].set_linewidth(axis_linewidth)
        ax.spines['bottom'].set_linewidth(axis_linewidth)
        ax.tick_params(width=axis_linewidth, labelsize=tick_label_fontsize, labelleft=True, labelbottom=True)

        if positive_col_data.empty:
            ax.text(0.5, 0.5, 'No positive data', ha='center', va='center', transform=ax.transAxes, fontsize=tick_label_fontsize)
            log_col_data_for_plot = pd.Series(dtype=float)
        else:
            log_col_data_for_plot = np.log10(positive_col_data)
            hist_counts, _, _ = ax.hist(
                log_col_data_for_plot, bins=bins, color=primary_hist_color,
                edgecolor='black', alpha=0.9, rwidth=0.9 # Rule: alpha 0.9 for lines/bars
            )
            if hist_counts.size > 0:
                 max_y_val = max(max_y_val, hist_counts.max())

        ax.set_title(col_name, fontsize=axis_label_fontsize, fontweight='semibold') # Rule: axis labels (here used for subplot title)

        if not log_col_data_for_plot.empty:
            tick_locator = mticker.MaxNLocator(nbins=5, prune='both')
            log_ticks = tick_locator.tick_values(log_col_data_for_plot.min(), log_col_data_for_plot.max())
            
            ax.set_xticks(log_ticks)
            ax.set_xticklabels([f"{10**val:.1e}" for val in log_ticks], rotation=45, ha="right", fontsize=tick_label_fontsize)
        else:
            ax.set_xticks([])
            ax.set_xticklabels([])


        if col_name == initial_time_point_col:
            if cutoff_val > 0:
                log_cutoff = np.log10(cutoff_val)
                current_xlim = ax.get_xlim()
                # Check if current_xlim is valid (sometimes can be (1.0,0.0) if no data initially)
                if current_xlim[0] < current_xlim[1] and log_cutoff >= current_xlim[0] and log_cutoff <= current_xlim[1]:
                    ax.axvline(log_cutoff, color='red', linestyle='--', linewidth=data_linewidth, label=f'Cutoff = {cutoff_val:.2g}')
                elif current_xlim[0] >= current_xlim[1]: # xlims are inverted or invalid
                     ax.axvline(log_cutoff, color='red', linestyle='--', linewidth=data_linewidth, label=f'Cutoff = {cutoff_val:.2g}')
                     logger.debug(f"Cutoff line for {col_name} plotted, but xlims {current_xlim} were potentially invalid before plotting.")
                else:
                    logger.info(f"Cutoff value {cutoff_val} (log10: {log_cutoff:.2f}) for {col_name} is outside plot x-limits ({current_xlim[0]:.2f}, {current_xlim[1]:.2f}). Line not drawn.")
                
                if ax.has_data(): # Only add legend if there's something to plot
                    ax.legend(fontsize=legend_fontsize, frameon=False) # Rule: legend 10pt, no frame
            else:
                logger.warning(f"Cutoff value ({cutoff_val}) for {initial_time_point_col} is not positive. Cutoff line not plotted on log scale.")
        
        if i > 0 :
            ax.set_yticklabels([])


    for ax_idx, ax in enumerate(axes):
        ax.set_ylim(0, max_y_val * 1.15 if max_y_val > 0 else 1)
        if ax_idx == 0:
            ax.set_ylabel('Frequency', fontsize=axis_label_fontsize, fontweight='semibold')
        ax.set_xlabel('Value (Original Scale)', fontsize=axis_label_fontsize, fontweight='semibold')
        ax.grid(True, linestyle=':', alpha=grid_alpha, linewidth=0.8) # Rule: grid alpha 0.3, line 0.8
        ax.tick_params(labelsize=tick_label_fontsize)


    # Calculate statistics
    stats: Dict[str, Any] = {}
    original_rows = len(df)
    original_counts = df[initial_time_point_col].sum() 

    stats['original_rows'] = original_rows
    stats['original_counts'] = original_counts

    if initial_time_point_col not in df.columns:
        logger.warning(f"Initial time point column '{initial_time_point_col}' not found in {filename}. Cannot calculate cutoff stats.")
        stats.update({'rows_kept': 'N/A', 'percentage_rows_kept': 'N/A',
                      'count_kept': 'N/A', 'percentage_count_kept': 'N/A'})
    elif not pd.api.types.is_numeric_dtype(df[initial_time_point_col]):
        logger.warning(f"Initial time point column '{initial_time_point_col}' in {filename} is not numeric. Cannot apply numeric cutoff.")
        stats.update({'rows_kept': 'N/A (col not numeric)', 'percentage_rows_kept': 'N/A',
                      'count_kept': 'N/A', 'percentage_count_kept': 'N/A'})
    else:
        # Ensure cutoff_val is float for comparison, already done by argparse
        kept_df = df[df[initial_time_point_col] >= cutoff_val]
        stats['rows_kept'] = len(kept_df)
        stats['percentage_rows_kept'] = (stats['rows_kept'] / original_rows) * 100.0 if original_rows > 0 else 0.0
        stats['count_kept'] = kept_df[initial_time_point_col].sum()
        stats['percentage_count_kept'] = (stats['count_kept'] / original_counts) * 100.0 if original_counts > 0 else 0.0
            
    # Add on-plot statistics text box (Rule: Statistical Annotations)
    stat_text_lines = [
        f"File: {filename}",
        f"Initial Time Point: '{initial_time_point_col}'",
        f"Cutoff Applied: >= {cutoff_val:.2g}",
        f"Original Rows: {stats['original_rows']:,}",
        f"Rows Kept: {stats['rows_kept'] if isinstance(stats['rows_kept'], str) else f'{stats['rows_kept']:,}'} ({stats['percentage_rows_kept'] if isinstance(stats['percentage_rows_kept'], str) else f'{stats['percentage_rows_kept']:.1f}%'})",
        f"Original Counts: {stats['original_counts']:,}",
        f"Counts Kept: {stats['count_kept'] if isinstance(stats['count_kept'], str) else f'{stats['count_kept']:,}'} ({stats['percentage_count_kept'] if isinstance(stats['percentage_count_kept'], str) else f'{stats['percentage_count_kept']:.1f}%'})"
    ]
    stat_text = "\n".join(stat_text_lines)
    
    # Position: Top-left corner (0.05, 0.95 in axes coordinates) - using figure coords for overall box
    # Anchored to top-left of the figure, slightly offset.
    fig.text(0.01, 0.97, stat_text, transform=fig.transFigure, fontsize=annotation_fontsize,
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(boxstyle='round,pad=0.4', fc='white', alpha=0.7, ec='gray'))

    plt.tight_layout(rect=[0.03, 0.03, 0.97, 0.88]) # Adjust rect for suptitle, common labels, and new stats box

    return fig, stats

def display_summary_table(aggregated_stats: List[Dict[str, Any]], total_time: float, headers: Dict[str, str]) -> None:
    """Displays a summary table of aggregated statistics using tabulate.

    Args:
        aggregated_stats (List[Dict[str, Any]]): A list of statistics dictionaries,
            one for each processed file. Each dictionary should include a 'filename' key.
        total_time (float): Total processing time for the script.
        headers (Dict[str, str]): A dictionary mapping statistic keys to header names for the table.
    """
    if not aggregated_stats:
        logger.info("No statistics to display.")
        return

    table_data = []
    for file_stats in aggregated_stats:
        row = []
        # If there was an error, we might just have filename and error message
        if "error" in file_stats and len(file_stats) <= 2:
            row.append(file_stats.get("filename", "Unknown File"))
            # Fill other columns with the error message or N/A if error is the only other field
            error_message = file_stats.get("error", "Processing Error")
            for i, key in enumerate(headers.keys()):
                if i == 0: # Already added filename
                    continue
                # Show error in one prominent column, e.g., original_rows, and N/A elsewhere
                if key == "original_rows": 
                    row.append(error_message[:50]) # Truncate long error messages
                else:
                    row.append("N/A") 
        else:
            for key in headers: 
                value = file_stats.get(key, 'N/A')
                if isinstance(value, float) and ('percentage' in key or '%' in headers[key]):
                    row.append(f"{value:.2f}%")
                elif isinstance(value, (int)) and not ('percentage' in key or '%' in headers[key]): # Check for int specifically
                    row.append(f"{value:,}") 
                elif isinstance(value, float) and not ('percentage' in key or '%' in headers[key]): # float values not for percentage
                    row.append(f"{value:.2f}") # format float to 2 decimal places
                else:
                    row.append(str(value)) # Ensure all values are strings for tabulate
        table_data.append(row)
    
    logger.info("\n--- Processing Summary ---")
    try:
        table_str = tabulate(table_data, headers=list(headers.values()), tablefmt="grid", stralign="left", numalign="right")
        for line in table_str.split('\n'):
            logger.info(line)
    except Exception as e:
        logger.error(f"Could not generate summary table with tabulate: {e}", exc_info=True)
        logger.info("Raw aggregated stats:")
        for stat_item in aggregated_stats:
            logger.info(str(stat_item))

def main() -> None:
    """Main execution function for the read count distribution analysis."""
    args = parse_arguments()
    setup_logging(args.verbose)

    start_time = time.time()
    logger.info("Starting Read Count Distribution Analysis script...")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    input_files = sorted(args.input, key=lambda p: p.name)

    logger.info(f"Processing {len(input_files)} input files.")
    logger.info(f"Initial time point column: '{args.initial_time_point}'")
    logger.info(f"Cutoff value: {args.cutoff}")
    logger.info(f"Histogram bins: {args.bins}")
    logger.info(f"Output PDF: {args.output}")

    all_file_stats: List[Dict[str, Any]] = []

    with PdfPages(args.output) as pdf:
        for file_path in input_files:
            filename = file_path.name
            logger.info(f"--- Processing file: {filename} ---")

            try:
                df_reader_kwargs = {'sep':'\t', 'index_col':[0,1,2,3]}
                # Simplified CSV reading based on user's diff (removed c-engine fallback try-except)
                df = pd.read_csv(file_path, engine='python', **df_reader_kwargs)
                
                if df.empty:
                    logger.warning(f"File {filename} resulted in an empty DataFrame after reading. Skipping.")
                    all_file_stats.append({"filename": filename, "error": "Empty DataFrame after read", **{k: 'N/A' for k in headers if k != 'filename'}}) # Add N/A for other stats fields
                    continue

                fig, current_stats = plot_distributions_and_calculate_stats(
                    df,
                    filename,
                    args.initial_time_point,
                    args.cutoff,
                    args.bins,
                    args.figsize_per_subplot_width,
                    args.figsize_subplot_height,
                    args.dpi
                )
                
                current_stats["filename"] = filename 
                all_file_stats.append(current_stats)

                if fig:
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close(fig)
                    logger.info(f"Plot generated for {filename}.")
                    logger.debug(f"Detailed stats for {filename}: {current_stats}")
                else:
                    logger.info(f"Plotting skipped for {filename} (e.g., no numeric data or other issue).")

            except FileNotFoundError:
                logger.error(f"File {filename} not found. Skipping.")
                all_file_stats.append({"filename": filename, "error": "File not found", **{k: 'N/A' for k in headers if k != 'filename'}})
                continue
            except pd.errors.EmptyDataError:
                logger.warning(f"File {filename} is empty. Skipping.")
                all_file_stats.append({"filename": filename, "error": "Empty file", **{k: 'N/A' for k in headers if k != 'filename'}})
                continue
            except pd.errors.ParserError as pe:
                logger.error(f"ParserError for {filename}: {pe}. Skipping.")
                all_file_stats.append({"filename": filename, "error": f"Parsing failed: {pe}", **{k: 'N/A' for k in headers if k != 'filename'}})
                continue
            except Exception as e:
                logger.error(f"An unexpected error occurred while processing {filename}: {e}", exc_info=True)
                all_file_stats.append({"filename": filename, "error": str(e), **{k: 'N/A' for k in headers if k != 'filename'}})
    
    end_time = time.time()
    total_processing_time = end_time - start_time

    logger.info(f"--- Analysis complete for all files. PDF saved to {args.output} ---")
    # Headers dict needed for display_summary_table and for filling N/A stats on error
    headers = {
        "filename": "File Name",
        "original_rows": "Original Rows",
        "original_counts": "Original Counts",
        "rows_kept": "Rows Kept",
        "percentage_rows_kept": "% Rows Kept",
        "count_kept": "Counts Kept",
        "percentage_count_kept": "% Counts Kept"
    } # This was inside display_summary_table, moving out for broader use
    display_summary_table(all_file_stats, total_processing_time, headers)

if __name__ == "__main__":
    main()
