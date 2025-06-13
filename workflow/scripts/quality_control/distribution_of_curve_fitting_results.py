"""
Distribution Analysis of Curve Fitting Results

This script reads a CSV file with statistical data and generates histogram plots
for all numeric columns, arranged in a 4-column subplot layout and saved as PDF.

Typical usage: python distribution_of_curve_fitting_results.py -i input.csv -o output.pdf
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from collections import defaultdict
import time

# Custom color palette following cursor rules
CUSTOM_COLORS = [
    '#962955',  # Deep pink-purple (primary)
    '#7fb775',  # Medium green (primary) 
    '#6479cc',  # Medium blue (primary)
    '#ad933c',  # Golden brown (primary)
    '#26b1fd',  # Bright blue
    '#8c397b',  # Medium purple
    '#9ab25d',  # Yellow-green
    '#48b7cd',  # Blue-cyan
    '#be6940',  # Orange-brown
    '#d3616e'   # Rose red
]

def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate histogram plots for numeric columns in CSV data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input',
        type=Path,
        required=True,
        help='Input CSV file with statistical data'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output PDF file for histogram plots'
    )
    
    parser.add_argument(
        '--bins',
        type=int,
        default=30,
        help='Number of bins for histograms'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    return parser.parse_args()

def load_and_analyze_data(input_file: Path) -> Tuple[pd.DataFrame, List[str], Dict[str, any]]:
    """
    Load CSV data and identify numeric columns.
    
    Args:
        input_file: Path to input CSV file
        
    Returns:
        Tuple of (dataframe, numeric_columns, statistics)
    """
    logging.info(f"Loading data from {input_file}")
    
    # Load data
    df = pd.read_csv(input_file)
    logging.info(f"Loaded {len(df)} rows and {len(df.columns)} columns")
    
    # Identify numeric columns
    numeric_columns = df.select_dtypes(include=[np.number]).columns.tolist()
    logging.info(f"Found {len(numeric_columns)} numeric columns: {numeric_columns}")
    
    # Generate basic statistics
    stats = {
        'total_rows': len(df),
        'total_columns': len(df.columns),
        'numeric_columns': len(numeric_columns),
        'missing_values': df.isnull().sum().sum()
    }
    
    return df, numeric_columns, stats

def setup_matplotlib_style() -> None:
    """Configure matplotlib with professional styling."""
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
        'font.size': 10,
        'axes.linewidth': 1.2,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linewidth': 0.8,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight'
    })

def create_histogram_plots(df: pd.DataFrame, numeric_columns: List[str], 
                         output_file: Path, bins: int = 30) -> Dict[str, any]:
    """
    Create histogram plots for all numeric columns.
    
    Args:
        df: Input dataframe
        numeric_columns: List of numeric column names
        output_file: Output PDF file path
        bins: Number of bins for histograms
        
    Returns:
        Dictionary with plotting statistics
    """
    setup_matplotlib_style()
    
    # Calculate subplot layout (4 columns)
    n_cols = 4
    n_rows = (len(numeric_columns) + n_cols - 1) // n_cols  # Ceiling division
    
    logging.info(f"Creating {n_rows}x{n_cols} subplot layout for {len(numeric_columns)} plots")
    
    # Create figure
    fig_width = 16  # 4 inches per column
    fig_height = n_rows * 4  # 4 inches per row
    
    with PdfPages(output_file) as pdf:
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
        
        # Handle case where we have only one row or column
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)
        elif len(numeric_columns) == 1:
            axes = np.array([[axes]])
        
        plot_stats = defaultdict(dict)
        
        for idx, column in enumerate(numeric_columns):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row, col]
            
            # Get data and remove NaN values
            data = df[column].dropna()
            
            if len(data) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(column, fontsize=12, fontweight='bold')
                continue
            
            # Create histogram
            color = CUSTOM_COLORS[idx % len(CUSTOM_COLORS)]
            n, bins_edges, patches = ax.hist(
                data, 
                bins=bins, 
                color=color, 
                alpha=0.8, 
                edgecolor='white', 
                linewidth=0.5
            )
            
            # Customize plot
            ax.set_title(column, fontsize=12, fontweight='bold', pad=10)
            ax.set_xlabel('Value', fontsize=10)
            ax.set_ylabel('Frequency', fontsize=10)
            
            # Add statistics text box
            stats_text = f'n = {len(data):,}\nMean = {data.mean():.3f}\nStd = {data.std():.3f}'
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', fontsize=9,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
            
            # Store statistics
            plot_stats[column] = {
                'count': len(data),
                'mean': data.mean(),
                'std': data.std(),
                'min': data.min(),
                'max': data.max(),
                'median': data.median()
            }
            
            logging.debug(f"Plotted {column}: {len(data)} values, range [{data.min():.3f}, {data.max():.3f}]")
        
        # Hide empty subplots
        for idx in range(len(numeric_columns), n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].set_visible(False)
        
        # Add overall title
        fig.suptitle('Distribution of Numeric Variables', fontsize=16, fontweight='bold', y=0.98)
        
        # Adjust layout
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # Save to PDF
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    logging.info(f"Histogram plots saved to {output_file}")
    
    return dict(plot_stats)

def display_summary_table(stats: Dict[str, any], plot_stats: Dict[str, Dict]) -> None:
    """Display summary statistics in table format."""
    try:
        from tabulate import tabulate
        
        # Overall statistics
        print("\n" + "="*60)
        print("SUMMARY STATISTICS")
        print("="*60)
        
        overall_data = [
            ["Total rows", f"{stats['total_rows']:,}"],
            ["Total columns", f"{stats['total_columns']:,}"],
            ["Numeric columns", f"{stats['numeric_columns']:,}"],
            ["Missing values", f"{stats['missing_values']:,}"],
            ["Plots generated", f"{len(plot_stats):,}"]
        ]
        
        print(tabulate(overall_data, headers=["Metric", "Value"], tablefmt="grid"))
        
        # Column-wise statistics
        if plot_stats:
            print(f"\nCOLUMN-WISE STATISTICS")
            print("="*60)
            
            column_data = []
            for col, col_stats in plot_stats.items():
                column_data.append([
                    col,
                    f"{col_stats['count']:,}",
                    f"{col_stats['mean']:.3f}",
                    f"{col_stats['std']:.3f}",
                    f"{col_stats['min']:.3f}",
                    f"{col_stats['max']:.3f}"
                ])
            
            headers = ["Column", "Count", "Mean", "Std", "Min", "Max"]
            print(tabulate(column_data, headers=headers, tablefmt="grid"))
    
    except ImportError:
        # Fallback to simple print if tabulate not available
        print(f"\nSummary: Generated {len(plot_stats)} histogram plots")
        print(f"Data: {stats['total_rows']:,} rows, {stats['numeric_columns']} numeric columns")

def main() -> None:
    """Main execution function."""
    start_time = time.time()
    
    # Parse arguments and setup logging
    args = parse_arguments()
    setup_logging(args.verbose)
    
    logging.info("Starting histogram generation script")
    
    try:
        # Validate input file
        if not args.input.exists():
            raise FileNotFoundError(f"Input file not found: {args.input}")
        
        # Create output directory if needed
        args.output.parent.mkdir(parents=True, exist_ok=True)
        
        # Load and analyze data
        df, numeric_columns, stats = load_and_analyze_data(args.input)
        
        if not numeric_columns:
            logging.warning("No numeric columns found in the input data")
            return
        
        # Create histogram plots
        plot_stats = create_histogram_plots(df, numeric_columns, args.output, args.bins)
        
        # Display summary
        execution_time = time.time() - start_time
        stats['execution_time'] = execution_time
        
        display_summary_table(stats, plot_stats)
        
        logging.info(f"Script completed successfully in {execution_time:.2f} seconds")
        
    except Exception as e:
        logging.error(f"Script failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()
