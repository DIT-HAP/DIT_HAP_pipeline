#!/usr/bin/env python3
"""
Gompertz Curve Fitting for Depletion Analysis

This script fits Gompertz growth curves to depletion time-series data from 
transposon insertion sequencing experiments. It processes multiple datasets
simultaneously and generates publication-quality plots with fitted parameters.

The Gompertz function models depletion curves with the equation:
y = A * exp(-exp((um * e / A) * (lam - x) + 1))

Where:
- A: maximum depletion level (asymptote)
- um: maximum depletion rate
- lam: lag time parameter

Typical usage: 
    python curve_fitting.py -i data.csv -t 0 2 4 6 8 10 12 14 -o results.csv

Input: CSV file with gene/insertion data as rows and time points as columns
Output: CSV file with fitted parameters and PDF with visualization plots

Author: Bioinformatics Pipeline
Date: 2024
"""

import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from scipy.optimize import minimize
from typing import List, Dict, Tuple, Optional, Union
import argparse
from collections import defaultdict
import time
from tqdm import tqdm

# Configure matplotlib for publication quality
plt.rcParams.update({
    'figure.max_open_warning': 0,
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
    'fit_line': '#962955',
    'data_points': '#6479cc',
    'constraint_lines': '#333333'
}


def setup_logging(verbose: bool = False) -> None:
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def gompertz_function(x: np.ndarray, A: float, um: float, lam: float) -> np.ndarray:
    """
    Calculate Gompertz function values with numerical stability.
    
    Args:
        x: Input time points
        A: Maximum depletion level (asymptote)
        um: Maximum depletion rate
        lam: Lag time parameter
        
    Returns:
        Gompertz function values at input points
    """
    # Check for division by zero
    if A == 0:
        return np.zeros_like(x)
    
    exponent = np.clip((um * np.e / A) * (lam - x) + 1, -700, 700)
    return A * np.exp(-np.exp(exponent))


def gompertz_derivative(x: np.ndarray, A: float, um: float, lam: float) -> np.ndarray:
    """
    Calculate derivative of Gompertz function.
    
    Args:
        x: Input time points
        A: Maximum depletion level
        um: Maximum depletion rate  
        lam: Lag time parameter
        
    Returns:
        Derivative values at input points
    """
    alpha = (um * np.e) / A
    u = alpha * (lam - x) + 1
    return A * alpha * np.exp(u - np.exp(u))


def objective_function(params: List[float], x: np.ndarray, y: np.ndarray) -> float:
    """
    Objective function for curve fitting using Huber loss.
    
    Args:
        params: Gompertz parameters [A, um, lam]
        x: Time point values
        y: Depletion measurements
        
    Returns:
        Huber loss value
    """
    A, um, lam = params
    y_fit = gompertz_function(x, A, um, lam)
    residuals = y - y_fit
    z = residuals**2
    # Huber loss for robustness to outliers
    rho_z = np.where(z <= 1, z, 2*np.sqrt(z) - 1)
    return np.sum(rho_z)


def constraint_function1(params: List[float], t_last: float) -> float:
    """Constraint to ensure reasonable parameter bounds."""
    A, um, lam = params
    return t_last + 3 - abs(A)/abs(um) - lam


def constraint_function2(params: List[float]) -> float:
    """Constraint to ensure smooth curve behavior."""
    A, um, lam = params
    x0 = lam + A/um/np.e
    return (abs(gompertz_derivative(x0-1, A, um, lam)) + 
            abs(gompertz_derivative(x0+1, A, um, lam)) - 1.8 * abs(um))


def fit_single_curve(x_values: np.ndarray, y_values: np.ndarray, 
                    gene_name: str, t_last: float) -> Dict[str, Union[str, float]]:
    """
    Fit Gompertz curve to a single dataset.
    
    Args:
        x_values: Time points
        y_values: Depletion measurements
        gene_name: Gene/insertion identifier
        t_last: Last time point for constraints
        
    Returns:
        Dictionary with fitting results and statistics
    """
    constraints = (
        {'type': 'ineq', 'fun': constraint_function1, 'args': (t_last,)},
        {'type': 'ineq', 'fun': constraint_function2}
    )
    
    try:
        result = minimize(
            objective_function,
            x0=[1, 1, 1],
            args=(x_values, y_values),
            bounds=((-1, t_last), (-1, np.inf), (-1e-6, t_last)),
            constraints=constraints,
            options={'maxiter': 3000, 'disp': False}
        )

        if result.success:
            A, um, lam = result.x
            residuals = y_values - gompertz_function(x_values, A, um, lam)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y_values - np.mean(y_values))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(ss_res / len(y_values))

            return {
                'gene': gene_name,
                'Status': 'Success',
                'A': A, 'um': um, 'lam': lam,
                'R2': r_squared, 'RMSE': rmse
            }
        else:
            logging.warning(f"Optimization failed for {gene_name}")
            return {
                'gene': gene_name,
                'Status': 'Optimization failed',
                'A': np.nan, 'um': np.nan, 'lam': np.nan,
                'R2': np.nan, 'RMSE': np.nan
            }
    
    except Exception as e:
        logging.error(f"Error fitting {gene_name}: {e}")
        return {
            'gene': gene_name,
            'Status': 'Fitting error',
            'A': np.nan, 'um': np.nan, 'lam': np.nan,
            'R2': np.nan, 'RMSE': np.nan
        }


def create_fitted_plot(ax: plt.Axes, x_values: np.ndarray, y_values: np.ndarray, 
                      params: Dict[str, Union[str, float]], gene_name: str) -> None:
    """
    Create a publication-quality plot for fitted curve.
    
    Args:
        ax: Matplotlib axes object
        x_values: Time points
        y_values: Depletion measurements
        params: Fitted parameters
        gene_name: Gene identifier
    """
    ax.grid(True, alpha=0.3)
    
    if params['Status'] == 'Success':
        A, um, lam = params['A'], params['um'], params['lam']
        
        # Plot data points
        ax.scatter(x_values, y_values, 
                  color=COLOR_PALETTE['data_points'], 
                  s=50, alpha=0.8, 
                  edgecolors='white', linewidth=0.5,
                  label='Data')
        
        # Plot fitted curve
        x_smooth = np.linspace(min(x_values), max(x_values), 100)
        y_fit = gompertz_function(x_smooth, A, um, lam)
        ax.plot(x_smooth, y_fit, 
               color=COLOR_PALETTE['fit_line'], 
               linewidth=2.0, label='Fitted')
        
        # Add constraint lines with subtle styling
        ax.axhline(y=A, color=COLOR_PALETTE['constraint_lines'], 
                  linestyle='--', alpha=0.3, linewidth=1.0)
        ax.axvline(x=lam, color=COLOR_PALETTE['constraint_lines'], 
                  linestyle='--', alpha=0.3, linewidth=1.0)
        
        # Add parameter text
        param_text = f'A={A:.2f}, μₘ={um:.2f}, λ={lam:.2f}\nR²={params["R2"]:.3f}'
        ax.text(0.05, 0.95, param_text, 
               transform=ax.transAxes, fontsize=8,
               verticalalignment='top', 
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    else:
        # Plot failed fit
        ax.scatter(x_values, y_values, 
                  color='gray', s=30, alpha=0.6)
        ax.text(0.5, 0.5, 'Fit Failed', 
               transform=ax.transAxes, fontsize=10,
               horizontalalignment='center', color='red')
    
    ax.set_ylim(-2, 10)
    ax.set_title(gene_name, fontsize=9, pad=5)
    ax.tick_params(labelsize=8)


def generate_fitting_plots(results_df: pd.DataFrame, x_values: np.ndarray, 
                          y_values: np.ndarray, output_plot: Path) -> None:
    """
    Generate multi-page PDF with fitting plots.
    
    Args:
        results_df: DataFrame with fitting results
        x_values: Time points
        y_values: Depletion data matrix
        output_plot: Output PDF path
    """
    plots_per_page = 32
    num_pages = int(np.ceil(len(results_df) / plots_per_page))
    
    logging.info(f"Generating {num_pages} pages of plots...")
    
    with PdfPages(output_plot) as pdf:
        for page in range(num_pages):
            fig, axes = plt.subplots(8, 4, figsize=(11.69, 8.27))  # A4 landscape
            axes = axes.flatten()
            
            start_idx = page * plots_per_page
            end_idx = min((page + 1) * plots_per_page, len(results_df))
            
            for idx in range(start_idx, end_idx):
                ax_idx = idx % plots_per_page
                row = results_df.iloc[idx]
                
                create_fitted_plot(
                    axes[ax_idx], 
                    x_values, 
                    y_values[idx], 
                    row.to_dict(),
                    row['gene']
                )
            
            # Hide unused subplots
            for ax_idx in range(end_idx - start_idx, plots_per_page):
                axes[ax_idx].set_visible(False)
            
            plt.tight_layout(pad=1.0)
            pdf.savefig(fig, bbox_inches='tight', dpi=300)
            plt.close(fig)


def process_depletion_data(input_file: Path, time_points: List[float]) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Load and process depletion data from CSV file.
    
    Args:
        input_file: Path to input CSV file
        time_points: List of time point values
        
    Returns:
        Tuple of (x_values, y_values, gene_names)
    """
    logging.info(f"Loading data from {input_file}")
    
    # Load data with multi-level index for insertions
    data = pd.read_csv(input_file, header=0)
    len_columns = len(data.columns)
    index_column_num = len_columns - len(time_points)
    index_columns = data.columns.tolist()[:index_column_num]
    timepoint_columns = data.columns.tolist()[index_column_num:]
    data.set_index(index_columns, inplace=True)

    # Create gene identifiers
    gene_names = ["_".join(map(str, idx)) for idx in data.index.tolist()]
    
    x_values = time_points
    y_values = data.values
    
    logging.info(f"Loaded {len(gene_names)} datasets with {len(x_values)} time points")
    
    return x_values, y_values, gene_names


def generate_summary_statistics(results_df: pd.DataFrame) -> Dict[str, Union[int, float]]:
    """
    Generate comprehensive summary statistics.
    
    Args:
        results_df: DataFrame with fitting results
        
    Returns:
        Dictionary with summary statistics
    """
    total_count = len(results_df)
    success_count = len(results_df[results_df['Status'] == 'Success'])
    success_rate = (success_count / total_count * 100) if total_count > 0 else 0
    
    # Statistics for successful fits only
    successful_fits = results_df[results_df['Status'] == 'Success']
    
    stats = {
        'Total datasets': total_count,
        'Successful fits': success_count,
        'Success rate (%)': success_rate,
        'Failed fits': total_count - success_count
    }
    
    if len(successful_fits) > 0:
        stats.update({
            'Mean R²': successful_fits['R2'].mean(),
            'Median R²': successful_fits['R2'].median(),
            'Mean RMSE': successful_fits['RMSE'].mean(),
            'Mean A parameter': successful_fits['A'].mean(),
            'Mean um parameter': successful_fits['um'].mean()
        })
    
    return stats


def display_summary_table(stats: Dict[str, Union[int, float]]) -> None:
    """
    Display summary statistics in formatted table.
    
    Args:
        stats: Dictionary with summary statistics
    """
    logging.info("\n" + "="*50)
    logging.info("CURVE FITTING SUMMARY STATISTICS")
    logging.info("="*50)
    
    for key, value in stats.items():
        if isinstance(value, float):
            logging.info(f"{key:<25}: {value:.3f}")
        else:
            logging.info(f"{key:<25}: {value}")
    
    logging.info("="*50)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='''
        Fit Gompertz curves to depletion time-series data.
        
        This script processes transposon insertion depletion data and fits
        Gompertz growth curves to model depletion kinetics. Generates both
        parameter estimates and visualization plots.
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--input', type=Path, required=True,
                       help='Input CSV file with depletion data')
    parser.add_argument('-t', '--time_points', required=True, nargs='+', 
                       type=float, help='Time points for the experiment')
    parser.add_argument('-o', '--output', type=Path, required=True,
                       help='Output CSV file for fitted parameters')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose logging')
    
    return parser.parse_args()


def main() -> None:
    """Main execution function."""
    start_time = time.time()
    
    # Parse arguments and setup
    args = parse_arguments()
    setup_logging(args.verbose)
    
    logging.info("Starting Gompertz curve fitting analysis")
    logging.info(f"Input file: {args.input}")
    logging.info(f"Time points: {args.time_points}")
    
    # Create output directories
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_plot = args.output.with_suffix('.pdf').with_name(args.output.stem + '_fitted_curves.pdf')
    
    # Process data
    x_values, y_values, gene_names = process_depletion_data(args.input, args.time_points)
    t_last = x_values[-1]
    
    # Fit curves with progress tracking
    logging.info("Fitting Gompertz curves...")
    all_results = []
    
    with tqdm(total=len(y_values), desc="Fitting progress") as pbar:
        for i, (y_data, gene_name) in enumerate(zip(y_values, gene_names)):
            result = fit_single_curve(x_values, y_data, gene_name, t_last)
            
            # Add time series data to result
            for j, time_val in enumerate(x_values):
                result[f't{j}'] = y_data[j]
            
            all_results.append(result)
            pbar.update(1)
    
    # Create results DataFrame
    results_df = pd.DataFrame(all_results)
    results_df.insert(1, 'time_points', [list(x_values)] * len(results_df))
    
    # Round numeric columns
    numeric_columns = ['A', 'um', 'lam', 'R2', 'RMSE']
    results_df[numeric_columns] = results_df[numeric_columns].round(3)

    # Save results
    results_df.to_csv(args.output, index=False)
    
    # Generate plots
    generate_fitting_plots(results_df, x_values, y_values, output_plot)
    
    # Calculate and display statistics
    stats = generate_summary_statistics(results_df)
    display_summary_table(stats)
    
    # Final summary
    elapsed_time = time.time() - start_time
    logging.info(f"\nAnalysis completed in {elapsed_time:.1f} seconds")
    logging.info(f"Results saved to: {args.output}")
    logging.info(f"Plots saved to: {output_plot}")


if __name__ == "__main__":
    main()