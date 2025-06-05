"""
PBL-PBR Correlation Analysis Script

This script analyzes correlation between PBL and PBR values from multiple TSV files.
It generates scatter plots with diagonal reference lines, calculates statistical
metrics (PCC, R²), and saves results to a PDF with rasterization.
"""

import argparse
import os
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze correlation between PBL and PBR from multiple TSV files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i file1.tsv file2.tsv file3.tsv -o correlation_analysis.pdf
  %(prog)s -i *.tsv -o results/correlation_plot.pdf
        """
    )
    
    parser.add_argument(
        '-i', '--input', 
        nargs='+', 
        required=True,
        help='Input TSV files (space-separated)'
    )
    
    parser.add_argument(
        '-o', '--output', 
        required=True,
        help='Output PDF file path'
    )
    
    parser.add_argument(
        '--dpi', 
        type=int, 
        default=300,
        help='Resolution for rasterization (default: 300)'
    )
    
    parser.add_argument(
        '--figsize', 
        nargs=2, 
        type=float, 
        default=[8, 8],
        help='Figure size in inches (width height, default: 10 8)'
    )
    
    return parser.parse_args()

def read_tsv_file(file_path):
    """Read TSV file and validate required columns."""
    try:
        df = pd.read_csv(file_path, sep='\t', index_col=[0,1,2])
        
        # Check if PBL and PBR columns exist
        required_cols = ['PBL', 'PBR']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"Warning: Missing columns {missing_cols} in {file_path}")
            return None
        
        # Remove rows with missing values in PBL or PBR
        df_clean = df[['PBL', 'PBR']].dropna()
        
        # Remove zero or negative values for log scaling
        df_clean = df_clean[(df_clean['PBL'] > 0) & (df_clean['PBR'] > 0)]
        
        if df_clean.empty:
            print(f"Warning: No valid data points in {file_path}")
            return None
        
        return df_clean
        
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return None

def calculate_statistics(x, y):
    """Calculate correlation statistics."""
    # Pearson correlation coefficient
    pcc = np.corrcoef(x, y)[0, 1]
    p_value = 0
    
    # R-squared
    r_squared = pcc ** 2
    
    # Linear regression for trend line
    slope, intercept = np.polyfit(np.log10(x), np.log10(y), 1)
    
    return {
        'pcc': pcc,
        'p_value': p_value,
        'r_squared': r_squared,
        'slope': slope,
        'intercept': intercept,
    }

def create_correlation_plot(filename, df, stats_info, figsize, dpi):
    """Create correlation plot for a single file with statistics."""
    # Set up the plot style
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Set white background with grid
    ax.set_facecolor('white')
    
    # Plot data points
    ax.scatter(df['PBL'], df['PBR'], 
              alpha=0.7, 
              s=40, 
              color='steelblue', 
              rasterized=True,
              edgecolors='navy',
              linewidth=0.5)
    
    # Set log scale for both axes
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Get axis limits for diagonal line
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Plot diagonal reference line (y=x)
    min_val = max(min(xlim), min(ylim))
    max_val = min(max(xlim), max(ylim))
    ax.plot([min_val, max_val], [min_val, max_val], 
            'k--', alpha=0.8, linewidth=2, label='y=x diagonal')
    
    # Add regression line
    x_line = np.logspace(np.log10(min_val), np.log10(max_val), 100)
    y_line = 10**(stats_info['intercept'] + stats_info['slope'] * np.log10(x_line))
    ax.plot(x_line, y_line, 'r-', alpha=0.8, linewidth=2, label='Regression line')
    
    # Customize the plot
    ax.set_xlabel('PBL (log scale)', fontsize=12, fontweight='bold')
    ax.set_ylabel('PBR (log scale)', fontsize=12, fontweight='bold')
    ax.set_title(f'PBL vs PBR Correlation Analysis\n{filename}', fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add legend
    ax.legend(loc='lower right')
    
    # Add statistics text box
    stats_text = []
    stats_text.append(f"Dataset: {filename}")
    stats_text.append(f"Data points: {len(df):,}")
    stats_text.append(f"PCC: {stats_info['pcc']:.4f}")
    stats_text.append(f"R²: {stats_info['r_squared']:.4f}")
    stats_text.append(f"Log-slope: {stats_info['slope']:.3f}")
    stats_text.append(f"Log-intercept: {stats_info['intercept']:.3f}")
    
    # Add text box with statistics
    textstr = '\n'.join(stats_text)
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=props, fontfamily='monospace')
    
    plt.tight_layout()
    return fig

def main():
    """Main function to run the correlation analysis."""
    # Parse arguments
    args = parse_arguments()
    
    print("=== PBL-PBR Correlation Analysis ===")
    print(f"Processing {len(args.input)} input files...")
    
    input_files = [str(input_file) for input_file in args.input]
    # Sort files by name
    input_files.sort(key=lambda x: Path(x).name)
    print(f"Processing files in order: {[Path(f).name for f in input_files]}")
    
    # Read and process files
    data_dict = {}
    stats_dict = {}
    
    for file_path in input_files:
        filename = os.path.basename(file_path)
        print(f"Reading {filename}...")
        
        df = read_tsv_file(file_path)
        if df is not None:
            data_dict[filename] = df
            
            # Calculate statistics
            stats_info = calculate_statistics(df['PBL'], df['PBR'])
            stats_dict[filename] = stats_info
            
            print(f"  - {len(df)} valid data points")
            print(f"  - PCC: {stats_info['pcc']:.3f}")
            print(f"  - R²: {stats_info['r_squared']:.3f}")
    
    if not data_dict:
        print("Error: No valid data found in any input file!")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    # Create and save plots
    print(f"Creating correlation plots...")
    
    # Save to PDF with rasterization
    print(f"Saving plots to {args.output}...")
    try:
        with PdfPages(args.output) as pdf:
            for filename, df in data_dict.items():
                print(f"  - Processing {filename}...")
                stats_info = stats_dict[filename]
                fig = create_correlation_plot(filename, df, stats_info, args.figsize, args.dpi)
                pdf.savefig(fig, dpi=args.dpi, bbox_inches='tight', 
                           facecolor='white', edgecolor='none')
                plt.close(fig)  # Close figure to free memory
        
        print(f"✓ Analysis complete! Output saved to: {args.output}")
        print(f"✓ Generated {len(data_dict)} correlation plots in PDF")
        
        # Print summary statistics
        print("\n=== Summary Statistics ===")
        total_points = sum(len(df) for df in data_dict.values())
        print(f"Total data points analyzed: {total_points}")
        print(f"Files processed: {len(data_dict)}")
        
        # Print individual file statistics
        for filename, stats_info in stats_dict.items():
            df = data_dict[filename]
            print(f"\n{filename}:")
            print(f"  Data points: {len(df):,}")
            print(f"  PCC: {stats_info['pcc']:.4f}")
            print(f"  R²: {stats_info['r_squared']:.4f}")
        
        # Overall statistics if multiple files
        if len(data_dict) > 1:
            all_pbl = pd.concat([df['PBL'] for df in data_dict.values()])
            all_pbr = pd.concat([df['PBR'] for df in data_dict.values()])
            overall_stats = calculate_statistics(all_pbl, all_pbr)
            print(f"\nOverall (combined data):")
            print(f"  PCC: {overall_stats['pcc']:.4f}")
            print(f"  R²: {overall_stats['r_squared']:.4f}")
    
    except Exception as e:
        print(f"Error saving plots: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
