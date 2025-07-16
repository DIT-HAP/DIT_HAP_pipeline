"""
Piecewise Linear Function Curve Fitting for Depletion Analysis

This script fits Piecewise Linear Function growth curves to depletion time-series data from 
transposon insertion sequencing experiments. It processes multiple datasets
simultaneously and generates publication-quality plots with fitted parameters.

The Piecewise Linear Function models depletion curves with the equation:
y = 0, x < lam
y = um * (x - lam), lam <= x <= lam + A/um
y = A, x > lam + A/um

Where:
- A: maximum depletion level (asymptote)
- um: depletion rate
- lam: lag time parameter

Typical usage: 
    python Piecewise_Linear_Function_curve_fitting.py -i data.csv -t 0 2 4 6 8 10 12 14 -o results.csv

Input: CSV file with gene/insertion data as rows and time points as columns
Output: CSV file with fitted parameters and PDF with visualization plots
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


def piecewise_linear_function(x: np.ndarray, A: float, um: float, lam: float) -> np.ndarray:
    """
    Piecewise Linear Function
    """
    # Avoid division by zero when um is zero or very close to zero
    y = np.piecewise(x, 
                        [
                            x < lam, 
                            (x >= lam) & (x <= lam + abs(A)/(abs(um)+1e-10)),
                            x > lam + abs(A)/(abs(um)+1e-10)
                        ], 
                        [
                            0, 
                            lambda x: um * (x - lam),
                            A
                        ]
                    )
    return y


def objective_function(params: List[float], x: np.ndarray, y: np.ndarray, weight_values: np.ndarray) -> float:
    """
    Objective function for curve fitting using Huber loss.
    
    Args:
        params: Piecewise Linear Function parameters [A, um, lam]
        x: Time point values
        y: Depletion measurements
        weight_values: Weight values
        
    Returns:
        Huber loss value
    """
    A, um, lam = params
    y_fit = piecewise_linear_function(x, A, um, lam)
    residuals = y - y_fit
    z = (residuals*weight_values)**2
    # Huber loss for robustness to outliers
    rho_z = np.where(z <= 1, z, 2*np.sqrt(z) - 1)
    return np.sum(rho_z)

def constraint_function(params: List[float]) -> float:
    """
    Constraint function for curve fitting.
    """
    A, um, lam = params
    return A*um

def constraint_function_2(params: List[float]) -> float:
    """
    Constraint function for curve fitting.
    """
    A, um, lam = params

    return abs(A)/(abs(um)+1e-10) - 3

def fit_single_curve(x_values: np.ndarray, y_values: np.ndarray, 
                    weight_values: np.ndarray, ID: str, t_last: float) -> Dict[str, Union[str, float]]:
    """
    Fit Piecewise Linear Function curve to a single dataset.
    
    Args:
        x_values: Time points
        y_values: Depletion measurements
        weight_values: Weight values
        ID: Gene/insertion identifier
        t_last: Last time point for constraints
        
    Returns:
        Dictionary with fitting results and statistics
    """
    
    try:

        cons = ({'type': 'ineq', 'fun': constraint_function},
                {'type': 'ineq', 'fun': constraint_function_2})

        result = minimize(
            objective_function,
            x0=[1, 1, 1],
            args=(x_values, y_values, weight_values),
            bounds=((-1, t_last), (-0.5, 2), (-1e-6, t_last)),
            options={'maxiter': 3000, 'disp': False},
            constraints=cons
        )

        if result.success:
            A, um, lam = result.x
            residuals = y_values - piecewise_linear_function(x_values, A, um, lam)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y_values - np.mean(y_values))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(ss_res / len(y_values))
            normalized_rmse = rmse / (y_values.max() - y_values.min())

            return {
                'ID': ID,
                'Status': 'Success',
                'A': A,
                'um': um, 'lam': lam,
                'R2': r_squared, 'RMSE': rmse, 'normalized_RMSE': normalized_rmse
            }
        else:
            logging.warning(f"Optimization failed for {ID}")
            return {
                'ID': ID,
                'Status': 'Optimization failed',
                'A': np.nan,
                'um': np.nan, 'lam': np.nan,
                'R2': np.nan, 'RMSE': np.nan, 'normalized_RMSE': np.nan
            }
    
    except Exception as e:
        logging.error(f"Error fitting {ID}: {e}")
        return {
            'ID': ID,
            'Status': 'Fitting error',
            'A': np.nan,
            'um': np.nan, 'lam': np.nan,
            'R2': np.nan, 'RMSE': np.nan, 'normalized_RMSE': np.nan
        }


def create_fitted_plot(ax: plt.Axes, x_values: np.ndarray, y_values: np.ndarray, 
                      params: Dict[str, Union[str, float]], ID: str) -> None:
    """
    Create a publication-quality plot for fitted curve.
    
    Args:
        ax: Matplotlib axes object
        x_values: Time points
        y_values: Depletion measurements
        params: Fitted parameters
        ID: Gene/insertion identifier
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
        y_fit = piecewise_linear_function(x_smooth, A, um, lam)
        ax.plot(x_smooth, y_fit, 
               color=COLOR_PALETTE['fit_line'], 
               linewidth=2.0, label='Fitted')
        
        # Add constraint lines with subtle styling
        ax.axhline(y=A, color=COLOR_PALETTE['constraint_lines'], 
                  linestyle='--', alpha=0.3, linewidth=1.0)
        ax.axvline(x=lam, color=COLOR_PALETTE['constraint_lines'], 
                  linestyle='--', alpha=0.3, linewidth=1.0)
        # add the straight line that crosses with the x axis at x=lam, and with the slope of um
        x_line = np.linspace(lam-1, 14, 100)
        y_line = um * (x_line - lam)
        ax.plot(x_line, y_line, color=COLOR_PALETTE['constraint_lines'], linestyle='--', alpha=0.3, linewidth=1.0)
        
        # Add parameter text
        param_text = f'A={A:.2f}    R²={params["R2"]:.3f}\num={um:.2f}  RMSE={params["RMSE"]:.3f}\nlam={lam:.2f}    NRMSE={params["normalized_RMSE"]:.3f}'
        ax.text(0.05, 0.95, param_text, 
               transform=ax.transAxes, fontsize=8,
               verticalalignment='top')
    else:
        # Plot failed fit
        ax.scatter(x_values, y_values, 
                  color='gray', s=30, alpha=0.6)
        ax.text(0.5, 0.5, 'Fit Failed', 
               transform=ax.transAxes, fontsize=10,
               horizontalalignment='center', color='red')
    
    ax.set_ylim(-1.5, 8.5)
    ax.set_title(" ".join(ID.split("=")), fontsize=9, pad=5)
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
            fig, axes = plt.subplots(8, 4, figsize=(8.27, 11.69))  # A4 landscape
            axes = axes.flatten()

            # progress bar
            if page % 10 == 0:
                logging.info(f"Generating page {page+1} of {num_pages}...")
            
            start_idx = page * plots_per_page
            end_idx = min((page + 1) * plots_per_page, len(results_df))
            
            for idx in range(start_idx, end_idx):
                ax_idx = idx % plots_per_page
                row = results_df.iloc[idx]
                ID = " ".join(map(str, row.name))
                
                create_fitted_plot(
                    axes[ax_idx], 
                    x_values, 
                    y_values[idx], 
                    row.to_dict(),
                    ID
                )
            
            # Hide unused subplots
            for ax_idx in range(end_idx - start_idx, plots_per_page):
                axes[ax_idx].set_visible(False)
            
            plt.tight_layout(pad=1.0)
            pdf.savefig(fig, bbox_inches='tight', dpi=300)
            plt.close(fig)


def process_depletion_data(input_file: Path, time_points: List[float], weight_file: Optional[Path] = None) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Load and process depletion data from CSV file.
    
    Args:
        input_file: Path to input CSV file
        time_points: List of time point values
        weight_file: Path to weight CSV file
        
    Returns:
        Tuple of (x_values, y_values, gene_names)
    """
    logging.info(f"Loading data from {input_file}")
    
    # Load data with multi-level index for insertions
    data = pd.read_csv(input_file, header=0, sep="\t")
    len_columns = len(data.columns)
    index_column_num = len_columns - len(time_points)
    index_columns = data.columns.tolist()[:index_column_num]
    timepoint_columns = data.columns.tolist()[index_column_num:]
    data.set_index(index_columns, inplace=True)

    index_names = data.index.names

    # Create gene identifiers
    IDs = ["=".join(map(str, idx)) for idx in data.index.tolist()]
    
    x_values = time_points
    y_values = data.values

    if weight_file is not None:
        weight_data = pd.read_csv(weight_file, header=0)
        weight_data.set_index(index_columns, inplace=True)
        weight_data = weight_data.loc[data.index].fillna(0.01)
        # if weights are lfcSE, then convert to 1/lfcSE^2
        if "lfcSE" in Path(weight_file).stem:
            weight_data = 1/(weight_data**2)
        weight_values = weight_data.values
    else:
        weight_values = np.ones(shape=(len(IDs), len(x_values)))
    
    logging.info(f"Loaded {len(IDs)} datasets with {len(x_values)} time points")
    
    return x_values, y_values, weight_values, IDs, index_names


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
            'Mean um parameter': successful_fits['um'].mean(),
            'Mean normalized RMSE': successful_fits['normalized_RMSE'].mean()
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


def main() -> None:
    """Main execution function."""
    start_time = time.time()
    
    input_file = Path("../results/HD_DIT_HAP/17_gene_level_depletion_analysis/LFC.tsv")
    time_points = [0, 2.352, 5.588, 9.104, 12.480]
    output_file = Path("../tmp/Piecewise_Linear_Function_curve_fitting/Gene_level_LFC_fitted.tsv")

    logging.info("Starting Piecewise Linear Function curve fitting analysis")
    logging.info(f"Input file: {input_file}")
    logging.info(f"Time points: {time_points}")
    
    # Create output directories
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_plot = output_file.with_suffix('.pdf').with_name(output_file.stem + '_fitted_curves.pdf')
    
    x_values, y_values, weight_values, IDs, index_names = process_depletion_data(input_file, time_points)
    t_last = x_values[-1]
    
    # Fit curves with progress tracking
    logging.info("Fitting Piecewise Linear Function curves...")
    all_results = []
    
    with tqdm(total=len(y_values), desc="Fitting progress") as pbar:
        for i, (y_data, ID) in enumerate(zip(y_values, IDs)):
            result = fit_single_curve(x_values, y_data, weight_values[i], ID, t_last)
            
            # Add time series data to result
            for j, time_val in enumerate(x_values):
                result[f't{j}'] = y_data[j]
            # Add the fitted time series data to the result
            try:
                for j, time_val in enumerate(x_values):
                    result[f't{j}_fitted'] = np.round(piecewise_linear_function(time_val, result['A'], result['um'], result['lam']), 3)
                for j, time_val in enumerate(x_values):
                    result[f't{j}_residual'] = np.round(result[f't{j}'] - result[f't{j}_fitted'], 3)
            except:
                for j, time_val in enumerate(x_values):
                    result[f't{j}_fitted'] = np.nan
                for j, time_val in enumerate(x_values):
                    result[f't{j}_residual'] = np.nan

            all_results.append(result)
            pbar.update(1)
    
    # Create results DataFrame
    results_df = pd.DataFrame(all_results)
    results_df.insert(1, 'time_points', [list(x_values)] * len(results_df))
    
    # Round numeric columns
    numeric_columns = ['A', 'um', 'lam', 'R2', 'RMSE', 'normalized_RMSE']
    results_df[numeric_columns] = results_df[numeric_columns].round(3)

    results_df.set_index("ID", inplace=True)
    multiple_index = pd.MultiIndex.from_tuples([ idx.split("=") for idx in results_df.index.tolist()])
    results_df.index = multiple_index
    results_df.rename_axis(index_names, inplace=True)

    # Save results
    results_df.to_csv(output_file, index=True, float_format='%.3f', sep="\t")
    
    # Generate plots
    generate_fitting_plots(results_df, x_values, y_values, output_plot)
    
    # Calculate and display statistics
    stats = generate_summary_statistics(results_df)
    display_summary_table(stats)
    # Final summary
    elapsed_time = time.time() - start_time
    logging.info(f"\nAnalysis completed in {elapsed_time:.1f} seconds")
    logging.info(f"Results saved to: {output_file}")
    logging.info(f"Plots saved to: {output_plot}")


if __name__ == "__main__":
    main()