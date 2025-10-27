from xmlrpc.server import resolve_dotted_attribute
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ================================ Constants =================================
SCRIPT_DIR = Path(__file__).parent.resolve()
TARGET_path = str((SCRIPT_DIR / "../../config/DIT_HAP.mplstyle").resolve())
plt.style.use(TARGET_path)
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# ================================ Functions =================================
def create_scatter_correlation_plot(
x: pd.Series | np.ndarray | list,
y: pd.Series | np.ndarray | list,
ax: plt.Axes,
xscale: None | str = None,
yscale: None | str = None,
) -> plt.Axes:
    """Create correlation plot for a single file with statistics."""    
    # Plot data points
    ax.scatter(
        x, y,
        alpha=0.5,
        s=10,
        facecolor="none",
        edgecolor=COLORS[1],
        rasterized=True
    )

    # Get axis limits for diagonal line
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Plot diagonal reference line (y=x)
    min_val = min(min(xlim), min(ylim))
    max_val = max(max(xlim), max(ylim))
    ax.plot([min_val, max_val], [min_val, max_val], 
            'k--', alpha=0.8, linewidth=2)
    
    # Set log scale for both axes
    if xscale == 'log':
        ax.set_xscale('log')
        x_for_fitting = np.log10(x)
    else:
        x_for_fitting = x

    if yscale == 'log':
        ax.set_yscale('log')
        y_for_fitting = np.log10(y)
    else:
        y_for_fitting = y

    # Calculate correlation statistics
    # Pearson correlation coefficient
    pcc = np.corrcoef(x_for_fitting, y_for_fitting)[0, 1]
    # R-squared
    r_squared = pcc**2
    # Linear regression
    slope, intercept = np.polyfit(x_for_fitting, y_for_fitting, 1)
    # RMSE
    y_pred = intercept + slope * x_for_fitting
    rmse = np.sqrt(np.mean((y_for_fitting - y_pred)**2))
    
    # Add statistics text box
    stats_text = []
    stats_text.append(f"Data points: {len(x):,}")
    stats_text.append(f"PCC: {pcc:.4f}")
    stats_text.append(f"R²: {r_squared:.4f}")
    stats_text.append(f"Slope: {slope:.4f}")
    stats_text.append(f"Intercept: {intercept:.4f}")
    stats_text.append(f"RMSE: {rmse:.4f}")

    # Add text box with statistics
    textstr = '\n'.join(stats_text)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, verticalalignment='top')
    
    return ax