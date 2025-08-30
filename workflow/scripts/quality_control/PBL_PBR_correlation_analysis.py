"""
PBL-PBR correlation analysis script for the DIT-HAP project.

This script analyzes the correlation between PBL and PBR values from multiple TSV files,
generating scatter plots with regression lines and statistical summaries.

Typical Usage:
    python PBL_PBR_correlation_analysis.py --input file1.tsv file2.tsv --output results.pdf

Input: TSV files containing PBL and PBR columns with multi-index structure
Output: PDF file containing correlation plots and statistical analysis
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Dict, Any, Optional
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# =============================== Constants ===============================
plt.style.use('/data/c/yangyusheng_optimized/DIT_HAP_pipeline/config/DIT_HAP.mplstyle')
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']


# =============================== Configuration & Models ===============================
class PBL_PBR_CorrelationAnalysisConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths for PBL-PBR correlation analysis."""
    input_files: List[Path] = Field(..., description="List of input TSV files")
    output_path: Path = Field(..., description="Path for output PDF file")

    @field_validator('input_files')
    def validate_input_files(cls, v):
        """Validate that all input files exist."""
        for file_path in v:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        return v
    
    @field_validator('output_path')
    def validate_output_path(cls, v):
        """Validate output path and create directory if needed."""
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True

class CorrelationStatistics(BaseModel):
    """Pydantic model to hold and validate correlation statistics results."""
    pcc: float = Field(..., ge=-1.0, le=1.0, description="Pearson correlation coefficient")
    p_value: float = Field(..., ge=0.0, le=1.0, description="P-value for correlation")
    r_squared: float = Field(..., ge=0.0, le=1.0, description="R-squared value")
    slope: float = Field(..., description="Log-scale slope")
    intercept: float = Field(..., description="Log-scale intercept")


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
def read_tsv_file(file_path: Path) -> Optional[pd.DataFrame]:
    """Read TSV file and validate required columns."""
    logger.info(f"Reading TSV file: {file_path}")
    
    try:
        df = pd.read_csv(file_path, sep='\t', index_col=[0,1,2])
        
        # Check if PBL and PBR columns exist
        required_cols = ['PBL', 'PBR']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            logger.warning(f"Warning: Missing columns {missing_cols} in {file_path}")
            return None
        
        # Remove rows with missing values in PBL or PBR
        df_clean = df[['PBL', 'PBR']].dropna()
        
        # Remove zero or negative values for log scaling
        df_clean = df_clean[(df_clean['PBL'] > 0) & (df_clean['PBR'] > 0)]
        
        if df_clean.empty:
            logger.warning(f"Warning: No valid data points in {file_path}")
            return None
        
        return df_clean
        
    except Exception as e:
        logger.error(f"Error reading {file_path}: {str(e)}")
        return None

@logger.catch
def calculate_correlation_statistics(x: pd.Series, y: pd.Series) -> CorrelationStatistics:
    """Calculate correlation statistics between two data series."""
    # Pearson correlation coefficient
    pcc = np.corrcoef(x, y)[0, 1]
    p_value = 0
    
    # R-squared
    r_squared = pcc ** 2
    
    # Linear regression for trend line
    slope, intercept = np.polyfit(np.log10(x), np.log10(y), 1)
    
    return CorrelationStatistics(
        pcc=pcc,
        p_value=p_value,
        r_squared=r_squared,
        slope=slope,
        intercept=intercept
    )

@logger.catch
def create_correlation_plot(filename: str, df: pd.DataFrame, stats_info: CorrelationStatistics) -> plt.Figure:
    """Create correlation plot for a single file with statistics."""
    fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))
    
    # Plot data points
    ax.scatter(df['PBL'], df['PBR'], 
              alpha=0.5, s=10, facecolor="none", edgecolor=COLORS[7],
              rasterized=True)
    
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
            'k--', alpha=0.8, linewidth=2)
    
    # Add regression line
    x_line = np.logspace(np.log10(min_val), np.log10(max_val), 100)
    y_line = 10**(stats_info.intercept + stats_info.slope * np.log10(x_line))
    ax.plot(x_line, y_line, 'r-', alpha=0.8, linewidth=2)
    
    # Customize the plot
    ax.set_xlabel('PBL (log scale)')
    ax.set_ylabel('PBR (log scale)')
    ax.set_title(f'PBL vs PBR Correlation Analysis\n{filename}')
    
    # Add grid
    ax.grid(True)
    
    # Add statistics text box
    stats_text = []
    stats_text.append(f"Dataset: {filename}")
    stats_text.append(f"Data points: {len(df):,}")
    stats_text.append(f"PCC: {stats_info.pcc:.4f}")
    stats_text.append(f"R²: {stats_info.r_squared:.4f}")
    stats_text.append(f"Log-slope: {stats_info.slope:.3f}")
    stats_text.append(f"Log-intercept: {stats_info.intercept:.3f}")
    
    # Add text box with statistics
    textstr = '\n'.join(stats_text)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, verticalalignment='top')
    
    return fig

# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(description="Analyze correlation between PBL and PBR from multiple TSV files")
    parser.add_argument("-i", "--input", nargs='+', type=Path, required=True, help="Input TSV files (space-separated)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output PDF file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point of the script."""
    
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    # Validate input and output paths using the Pydantic model
    try:
        config = PBL_PBR_CorrelationAnalysisConfig(
            input_files=args.input,
            output_path=args.output
        )

        logger.info("=== PBL-PBR Correlation Analysis ===")
        logger.info(f"Processing {len(config.input_files)} input files...")
        
        # Sort files by name
        sorted_files = sorted(config.input_files, key=lambda x: x.name)
        logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")
        
        # Read and process files
        data_dict = {}
        stats_dict = {}
        
        for file_path in sorted_files:
            filename = file_path.name
            logger.info(f"Reading {filename}...")
            
            df = read_tsv_file(file_path)
            if df is not None:
                data_dict[filename] = df
                
                # Calculate statistics
                stats_info = calculate_correlation_statistics(df['PBL'], df['PBR'])
                stats_dict[filename] = stats_info
                
                logger.info(f"  - {len(df)} valid data points")
                logger.info(f"  - PCC: {stats_info.pcc:.3f}")
                logger.info(f"  - R²: {stats_info.r_squared:.3f}")
        
        if not data_dict:
            logger.error("Error: No valid data found in any input file!")
            sys.exit(1)
        
        # Create and save plots
        logger.info("Creating correlation plots...")
        
        # Save to PDF with rasterization
        logger.info(f"Saving plots to {config.output_path}...")
        try:
            with PdfPages(config.output_path) as pdf:
                for filename, df in data_dict.items():
                    logger.info(f"  - Processing {filename}...")
                    stats_info = stats_dict[filename]
                    fig = create_correlation_plot(filename, df, stats_info)
                    pdf.savefig(fig)
                    plt.close(fig)  # Close figure to free memory
            
            logger.success(f"Analysis complete! Output saved to: {config.output_path}")
            logger.info(f"Generated {len(data_dict)} correlation plots in PDF")
            
            # Print summary statistics
            logger.info("\n=== Summary Statistics ===")
            total_points = sum(len(df) for df in data_dict.values())
            logger.info(f"Total data points analyzed: {total_points}")
            logger.info(f"Files processed: {len(data_dict)}")
            
            # Print individual file statistics
            for filename, stats_info in stats_dict.items():
                df = data_dict[filename]
                logger.info(f"\n{filename}:")
                logger.info(f"  Data points: {len(df):,}")
                logger.info(f"  PCC: {stats_info.pcc:.4f}")
                logger.info(f"  R²: {stats_info.r_squared:.4f}")
            
            # Overall statistics if multiple files
            if len(data_dict) > 1:
                all_pbl = pd.concat([df['PBL'] for df in data_dict.values()])
                all_pbr = pd.concat([df['PBR'] for df in data_dict.values()])
                overall_stats = calculate_correlation_statistics(all_pbl, all_pbr)
                logger.info("\nOverall (combined data):")
                logger.info(f"  PCC: {overall_stats.pcc:.4f}")
                logger.info(f"  R²: {overall_stats.r_squared:.4f}")
        
        except Exception as e:
            logger.error(f"Error saving plots: {str(e)}")
            sys.exit(1)
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
