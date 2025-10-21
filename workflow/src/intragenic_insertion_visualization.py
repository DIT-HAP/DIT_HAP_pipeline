# ================================ Imports =================================
from pathlib import Path
from typing import List, Optional
import pandas as pd
import matplotlib.pyplot as plt

# ================================= Constants =================================
SCRIPT_DIR = Path(__file__).parent.resolve()
plt.style.use(SCRIPT_DIR / "../../config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# ================================= Functions =================================
def intragenic_insertion_visualization(
    insertions_with_annotations: pd.DataFrame,
    gene_statistics: pd.DataFrame,
    gene: str,
    color: str,
    ax: plt.Axes,
    feature: str = "um",
) -> plt.Axes:
    """Visualize intragenic insertions for a given gene."""
    gene_insertions = insertions_with_annotations.query("Name == @gene")
    gene_metrics = gene_statistics.loc[gene]
    ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
    ax.axhline(gene_metrics[feature], color="darkred", linestyle="--", alpha=0.5)
    x = gene_insertions['Residue_affected']
    y = gene_insertions[feature]
    ax.scatter(x, y, marker='o', color=color, alpha=0.7)
    ax.set_title(gene, fontstyle='italic')
    ax.set_xlabel('Protein Position (residues)')

    match feature:
        case "um":
            feature_name = "DR"
        case "lam":
            feature_name = "DL"
        case _:
            raise ValueError(f"Invalid feature: {feature}")
    ax.set_ylabel(f"insertion-level {feature_name}")
    return ax