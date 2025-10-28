# ================================ Imports =================================
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

# ================================ Constants =================================
SCRIPT_DIR = Path(__file__).parent.resolve()
TARGET_path = str((SCRIPT_DIR / "../../config/DIT_HAP.mplstyle").resolve())
plt.style.use(TARGET_path)
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# ================================ Functions =================================
def plot_depletion_curves_for_given_genes(
    ax: Axes,
    data_df: pd.DataFrame,
    time_points: list[float],
    genes: list[str],
    gene_column: str,
    title: str,
    use_fitted: bool = False,
    **kwargs
) -> Axes:
    """Plot depletion curves for the given genes."""
    
    # Choose columns based on fitted or raw
    if use_fitted:
        value_cols = ['YES0_fitted', 'YES1_fitted', 'YES2_fitted', 'YES3_fitted', 'YES4_fitted']
    else:
        value_cols = ['YES0', 'YES1', 'YES2', 'YES3', 'YES4']
    
    # Plot individual curves in gray with transparency
    for gene in genes:
        values = data_df.query(f"{gene_column} == @gene")[value_cols].values.flatten().tolist()
        ax.plot(time_points, values, color='gray', alpha=0.5, linewidth=0.5, **kwargs)
    
    # Calculate and plot centroid in red
    subset_df = data_df.query(f"{gene_column} in @genes")
    if subset_df.shape[0] > 0:
        centroid = [subset_df[col].mean() for col in value_cols]
        ax.plot(time_points, centroid, color='red', linewidth=2, 
            label=f'Centroid (n={len(subset_df)})', **kwargs)
    
    ax.set_title(f'{title}\n(n={len(subset_df)})')
    ax.set_xlabel('Generations')
    ax.set_ylabel('LFC Value')
    ax.grid(True)

    return ax

def plot_depletion_curves_for_groups(
    data_df: pd.DataFrame,
    time_points: list[float],
    group_column: str,
    gene_column: str,
    col_num: int = 4,
    use_fitted: bool = False,
    **kwargs
) -> Figure:
    """Plot depletion curves for the groups based on the group column. """

    # Determine number of rows needed
    groups = sorted(data_df[group_column].unique().tolist())
    row_num = int(np.ceil(len(groups) / col_num))

    # Create subplots
    fig, axes = plt.subplots(row_num, col_num, figsize=(AX_WIDTH*col_num, AX_HEIGHT * row_num), sharex=True, sharey=True)
    if row_num == 1 and col_num == 1:
        axes = axes.reshape(1, -1)
    else:
        axes = axes.flatten()

    # Plot for each group
    for idx, group in enumerate(groups):
        # Get genes for the group
        group_df = data_df.query(f"{group_column} == @group")
        genes = group_df[gene_column].unique().tolist()
        
        # Get the corresponding axis
        ax = axes[idx]
        ax = plot_depletion_curves_for_given_genes(
            ax=ax,
            data_df=data_df,
            time_points=time_points,
            genes=genes,
            gene_column=gene_column,
            title=f'{group_column} {group}',
            use_fitted=use_fitted,
            **kwargs
        )
        ax.tick_params(axis="both", which="major", labelleft=True, labelbottom=True)

    # Remove unused subplots
    for j in range(idx + 1, len(axes)):
        fig.delaxes(axes[j])

    return fig

def plot_given_genes_on_feature_space(
    ax: Axes,
    data_df: pd.DataFrame,
    genes: list[str],
    gene_column: str,
    title: str,
    x_feature: str = "um",
    y_feature: str = "lam",
    **kwargs
) -> Axes:
    """ Plot given genes on feature space. """

    # all points in light gray
    x_all = data_df[x_feature]
    y_all = data_df[y_feature]
    ax.scatter(x_all, y_all, color='lightgray', alpha=0.4, **kwargs)

    # points for given genes
    subset_df = data_df.query(f"{gene_column} in @genes")
    x_subset = subset_df[x_feature]
    y_subset = subset_df[y_feature]
    try:
        xy_subset = np.vstack([x_subset, y_subset])
        z = gaussian_kde(xy_subset)(xy_subset)
        ax.scatter(x_subset, y_subset, c=z, cmap='viridis', **kwargs)
    except Exception:
        ax.scatter(x_subset, y_subset, color='red', **kwargs)

    ax.set_title(f'{title}\n(n={len(subset_df)})')
    ax.set_xlabel(x_feature)
    ax.set_ylabel(y_feature)
    ax.grid(True)

    return ax

def plot_groups_on_feature_space(
    data_df: pd.DataFrame,
    group_column: str,
    gene_column: str,
    col_num: int = 4,
    x_feature: str = "um",
    y_feature: str = "lam",
    **kwargs
) -> Figure:
    """ Plot groups on feature space based on the group column. """

    # Determine number of rows needed
    groups = sorted(data_df[group_column].unique().tolist())
    row_num = int(np.ceil(len(groups) / col_num))

    # Create subplots
    fig, axes = plt.subplots(row_num, col_num, figsize=(AX_WIDTH*col_num, AX_HEIGHT * row_num), sharex=True, sharey=True)
    if row_num == 1 and col_num == 1:
        axes = axes.reshape(1, -1)
    else:
        axes = axes.flatten()

    # Plot for each group
    for idx, group in enumerate(groups):
        # Get genes for the group
        group_df = data_df.query(f"{group_column} == @group")
        genes = group_df[gene_column].unique().tolist()
        
        # Get the corresponding axis
        ax = axes[idx]
        ax = plot_given_genes_on_feature_space(
            ax=ax,
            data_df=data_df,
            genes=genes,
            gene_column=gene_column,
            title=f'{group_column} {group}',
            x_feature=x_feature,
            y_feature=y_feature,
            **kwargs
        )
        ax.tick_params(axis="both", which="major", labelleft=True, labelbottom=True)

    # Remove unused subplots
    for j in range(idx + 1, len(axes)):
        fig.delaxes(axes[j])

    return fig