# Machine Learning Analysis Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Rewrite `machine_learning_analysis.ipynb` to run mljar-supervised AutoML regression for `um` and `lam` targets on the `um_gt_p35` subset, with full feature importance and SHAP analysis.

**Architecture:** Load pre-transformed data from `um_gt_p35_transformed_features_and_targets.csv`, run AutoML in Explain then Perform modes for each target, visualize predictions and residuals, extract and compare SHAP-based feature importance between `um` and `lam`.

**Tech Stack:** mljar-supervised AutoML, scikit-learn, matplotlib, pandas, numpy, SHAP (via mljar explain_level=2)

---

## Context

- **Notebook path:** `workflow/notebooks/machine_learning_analysis.ipynb`
- **Input data:** `results/HD_DIT_HAP_generationRAW/22_machine_learning_modeling/um_gt_p35_transformed_features_and_targets.csv`
  - Shape: (1421, 87) — 83 feature columns + `um`, `lam`, `A`, `DIT_HAP_cluster`
  - Features already transformed (PowerTransformer, StandardScaler, OneHotEncoder)
- **Output base dir:** `results/HD_DIT_HAP_generationRAW/22_machine_learning_modeling/regression/`
- **Target variables:** `um` and `lam` (regression only, not `A` or `DIT_HAP_cluster`)
- **Subset:** `um_gt_p35` only (already pre-filtered in the CSV)
- **AutoML modes:** First `Explain`, then `Perform` (separate cells, run independently)

---

### Task 1: Clear existing cells and write Cell 1 — Imports

**Files:**
- Modify: `workflow/notebooks/machine_learning_analysis.ipynb` (replace all existing cells)

**Step 1: Replace all notebook cells with a single imports cell**

Delete all existing cells and write the following as cell 0 (code cell):

```python
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr

from sklearn.model_selection import train_test_split
from supervised.automl import AutoML
```

**Step 2: Verify cell is correct by reading notebook**

Check the cell was written correctly.

---

### Task 2: Write Cell 2 — Config dataclass

**Files:**
- Modify: `workflow/notebooks/machine_learning_analysis.ipynb` (add cell)

**Step 1: Add the Config cell**

```python
@dataclass
class Config:
    target: str   # 'um' or 'lam'
    mode: str     # 'Explain' or 'Perform'

    data_path: Path = Path(
        "../../results/HD_DIT_HAP_generationRAW/"
        "22_machine_learning_modeling/"
        "um_gt_p35_transformed_features_and_targets.csv"
    )
    base_output_dir: Path = Path(
        "../../results/HD_DIT_HAP_generationRAW/"
        "22_machine_learning_modeling/regression"
    )
    test_size: float = 0.2
    random_state: int = 42
    index_column: str = "Systematic_ID"
    all_target_columns: list = field(
        default_factory=lambda: ['um', 'lam', 'A', 'DIT_HAP_cluster']
    )

    @property
    def output_dir(self) -> Path:
        path = self.base_output_dir / f"{self.target}_{self.mode.lower()}"
        path.mkdir(parents=True, exist_ok=True)
        return path

    def __post_init__(self):
        df = pd.read_csv(self.data_path, index_col=self.index_column)
        self.feature_columns = [
            c for c in df.columns if c not in self.all_target_columns
        ]
        self.data = df[self.feature_columns + [self.target]].dropna()
        print(f"Config: target={self.target}, mode={self.mode}")
        print(f"  Data shape: {self.data.shape}")
        print(f"  Features: {len(self.feature_columns)}")
        print(f"  Output dir: {self.output_dir}")
```

**Step 2: Verify by instantiating one Config in next cell temporarily**

```python
# quick smoke test - delete after confirming
cfg = Config(target='um', mode='Explain')
```

Expected output:
```
Config: target=um, mode=Explain
  Data shape: (1421, 84)
  Features: 83
  Output dir: .../regression/um_explain
```

---

### Task 3: Write Cell 3 — run_automl() function

**Files:**
- Modify: `workflow/notebooks/machine_learning_analysis.ipynb` (add cell)

**Step 1: Add the training function cell**

```python
def run_automl(cfg: Config):
    """Train AutoML and return model, test split, and predictions."""
    X = cfg.data[cfg.feature_columns]
    y = cfg.data[cfg.target]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=cfg.test_size,
        random_state=cfg.random_state
    )

    automl = AutoML(
        mode=cfg.mode,
        ml_task="regression",
        results_path=str(cfg.output_dir),
        explain_level=2,
    )
    automl.fit(X_train, y_train)
    y_pred = automl.predict(X_test)

    return automl, X_train, X_test, y_train, y_test, y_pred
```

---

### Task 4: Write Cell 4 — evaluate_and_plot() function

**Files:**
- Modify: `workflow/notebooks/machine_learning_analysis.ipynb` (add cell)

**Step 1: Add evaluation and plotting function**

```python
def evaluate_and_plot(cfg: Config, y_test, y_pred, automl):
    """Compute metrics, save prediction scatter + residual plots."""
    from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

    r2   = r2_score(y_test, y_pred)
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    mae  = mean_absolute_error(y_test, y_pred)
    r, p = pearsonr(y_test, y_pred)

    metrics = {"R2": r2, "RMSE": rmse, "MAE": mae, "Pearson_r": r, "Pearson_p": p}
    print(f"\n=== {cfg.target} [{cfg.mode}] ===")
    for k, v in metrics.items():
        print(f"  {k}: {v:.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Prediction vs actual
    ax = axes[0]
    ax.scatter(y_test, y_pred, s=20, alpha=0.5, color='#962955', edgecolor='none')
    lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
    ax.plot(lims, lims, 'k--', linewidth=1, alpha=0.6)
    ax.set_xlabel(f"True {cfg.target}", fontsize=13)
    ax.set_ylabel(f"Predicted {cfg.target}", fontsize=13)
    ax.set_title(f"Predicted vs True ({cfg.mode})", fontsize=14, fontweight='bold')
    ax.text(0.05, 0.92,
            f"R²={r2:.3f}\nr={r:.3f}\nRMSE={rmse:.3f}",
            transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8),
            fontsize=11)
    ax.grid(True, alpha=0.3)

    # Residuals
    residuals = y_test.values - y_pred
    ax = axes[1]
    ax.scatter(y_pred, residuals, s=20, alpha=0.5, color='#6479cc', edgecolor='none')
    ax.axhline(0, color='k', linestyle='--', linewidth=1, alpha=0.6)
    ax.set_xlabel(f"Predicted {cfg.target}", fontsize=13)
    ax.set_ylabel("Residual (True − Predicted)", fontsize=13)
    ax.set_title(f"Residuals ({cfg.mode})", fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.suptitle(f"Target: {cfg.target}  |  Mode: {cfg.mode}", fontsize=15, fontweight='bold')
    plt.tight_layout()
    out_path = cfg.output_dir / "prediction_and_residuals.pdf"
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"  Saved: {out_path}")

    return metrics
```

---

### Task 5: Write Cell 5 — plot_feature_importance() function

**Files:**
- Modify: `workflow/notebooks/machine_learning_analysis.ipynb` (add cell)

**Step 1: Add feature importance plotting function**

mljar AutoML writes a `features_importance.csv` to the output directory. Read it and plot top 20. Also look for SHAP images generated by mljar.

```python
def plot_feature_importance(cfg: Config, top_n: int = 20):
    """Read mljar feature importance CSV and plot top N features."""
    fi_path = cfg.output_dir / "features_importance.csv"
    if not fi_path.exists():
        print(f"  Feature importance file not found: {fi_path}")
        return None

    fi = pd.read_csv(fi_path).sort_values("mean_importance", ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.35)))
    colors_list = ['#962955' if i < 5 else '#8c397b' if i < 10 else '#b684d5'
                   for i in range(len(fi))]
    ax.barh(fi["feature"][::-1], fi["mean_importance"][::-1],
            color=colors_list[::-1], edgecolor='white', linewidth=0.5)
    ax.set_xlabel("Mean Importance", fontsize=13)
    ax.set_title(
        f"Top {top_n} Feature Importances — {cfg.target} [{cfg.mode}]",
        fontsize=14, fontweight='bold'
    )
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()

    out_path = cfg.output_dir / f"feature_importance_top{top_n}.pdf"
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"  Saved: {out_path}")

    return fi
```

---

### Task 6: Write Cell 6 — compare_feature_importance() function

**Files:**
- Modify: `workflow/notebooks/machine_learning_analysis.ipynb` (add cell)

**Step 1: Add comparison function (reads saved CSV, does not require models in memory)**

```python
def compare_feature_importance(
    cfg_um: Config,
    cfg_lam: Config,
    top_n: int = 20,
    base_output_dir: Path = Path(
        "../../results/HD_DIT_HAP_generationRAW/"
        "22_machine_learning_modeling/regression"
    )
):
    """Side-by-side feature importance for um vs lam, highlighting shared top features."""
    fi_um_path  = cfg_um.output_dir  / "features_importance.csv"
    fi_lam_path = cfg_lam.output_dir / "features_importance.csv"

    if not fi_um_path.exists() or not fi_lam_path.exists():
        print("Feature importance files missing — run AutoML for both targets first.")
        return

    fi_um  = pd.read_csv(fi_um_path).sort_values("mean_importance", ascending=False).head(top_n)
    fi_lam = pd.read_csv(fi_lam_path).sort_values("mean_importance", ascending=False).head(top_n)

    shared = set(fi_um["feature"]) & set(fi_lam["feature"])
    print(f"Shared top-{top_n} features: {len(shared)}")
    print(f"  {sorted(shared)}")

    fig, axes = plt.subplots(1, 2, figsize=(18, max(6, top_n * 0.35)))

    for ax, fi, label, color_hi, color_lo in [
        (axes[0], fi_um,  "um",  '#962955', '#b684d5'),
        (axes[1], fi_lam, "lam", '#6479cc', '#89afba'),
    ]:
        colors_bar = [color_hi if f in shared else color_lo for f in fi["feature"][::-1]]
        ax.barh(fi["feature"][::-1], fi["mean_importance"][::-1],
                color=colors_bar, edgecolor='white', linewidth=0.5)
        ax.set_xlabel("Mean Importance", fontsize=12)
        ax.set_title(f"Top {top_n} — {label} [{cfg_um.mode}]", fontsize=13, fontweight='bold')
        ax.grid(True, axis='x', alpha=0.3)

    # legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#962955', label='Shared (in both um & lam)'),
        Patch(facecolor='#b684d5', label='um only'),
        Patch(facecolor='#89afba', label='lam only'),
    ]
    axes[1].legend(handles=legend_elements, loc='lower right', fontsize=10)

    plt.suptitle(f"Feature Importance: um vs lam  [{cfg_um.mode}]",
                 fontsize=15, fontweight='bold')
    plt.tight_layout()

    out_path = base_output_dir / f"um_vs_lam_feature_importance_{cfg_um.mode.lower()}.pdf"
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved: {out_path}")
```

---

### Task 7: Write Cell 7 — Section header (markdown) + Explain mode execution

**Files:**
- Modify: `workflow/notebooks/machine_learning_analysis.ipynb` (add cells)

**Step 1: Add markdown cell**

```markdown
# AutoML — Explain Mode

Runs AutoML in `Explain` mode for `um` and `lam` targets.
Results saved to `regression/um_explain/` and `regression/lam_explain/`.

Re-run this section to retrain. Results from a previous run will be overwritten.
```

**Step 2: Add execution cell**

```python
# --- Explain mode ---
explain_metrics = {}
explain_fi = {}

for target in ['um', 'lam']:
    cfg = Config(target=target, mode='Explain')

    automl, X_train, X_test, y_train, y_test, y_pred = run_automl(cfg)
    explain_metrics[target] = evaluate_and_plot(cfg, y_test, y_pred, automl)
    explain_fi[target] = plot_feature_importance(cfg, top_n=20)

# Metrics summary table
print("\n=== Explain Mode — Summary ===")
summary = pd.DataFrame(explain_metrics).T
print(summary.round(4).to_string())
```

**Step 3: Add comparison cell**

```python
# Feature importance comparison (Explain)
cfg_um_explain  = Config(target='um',  mode='Explain')
cfg_lam_explain = Config(target='lam', mode='Explain')
compare_feature_importance(cfg_um_explain, cfg_lam_explain, top_n=20)
```

**Step 4: Run the execution cell and verify**

Expected: AutoML trains, prints metrics, saves PDFs to `regression/um_explain/` and `regression/lam_explain/`. Metrics table printed.

---

### Task 8: Write Cell 8 — Section header (markdown) + Perform mode execution

**Files:**
- Modify: `workflow/notebooks/machine_learning_analysis.ipynb` (add cells)

**Step 1: Add markdown cell**

```markdown
# AutoML — Perform Mode

Runs AutoML in `Perform` mode for `um` and `lam` targets.
**Note:** This section takes ~30–60 minutes to complete.

Results saved to `regression/um_perform/` and `regression/lam_perform/`.
```

**Step 2: Add execution cell (identical structure to Explain, different mode)**

```python
# --- Perform mode ---
perform_metrics = {}
perform_fi = {}

for target in ['um', 'lam']:
    cfg = Config(target=target, mode='Perform')

    automl, X_train, X_test, y_train, y_test, y_pred = run_automl(cfg)
    perform_metrics[target] = evaluate_and_plot(cfg, y_test, y_pred, automl)
    perform_fi[target] = plot_feature_importance(cfg, top_n=20)

# Metrics summary table
print("\n=== Perform Mode — Summary ===")
summary = pd.DataFrame(perform_metrics).T
print(summary.round(4).to_string())
```

**Step 3: Add comparison cell**

```python
# Feature importance comparison (Perform)
cfg_um_perform  = Config(target='um',  mode='Perform')
cfg_lam_perform = Config(target='lam', mode='Perform')
compare_feature_importance(cfg_um_perform, cfg_lam_perform, top_n=20)
```

---

### Task 9: Final verification

**Step 1: Verify notebook structure**

Notebook should have these cells in order:
1. Imports
2. Config dataclass
3. `run_automl()` function
4. `evaluate_and_plot()` function
5. `plot_feature_importance()` function
6. `compare_feature_importance()` function
7. Markdown: "AutoML — Explain Mode"
8. Explain execution cell
9. Explain comparison cell
10. Markdown: "AutoML — Perform Mode"
11. Perform execution cell
12. Perform comparison cell

**Step 2: Verify output files exist after Explain run**

```bash
ls results/HD_DIT_HAP_generationRAW/22_machine_learning_modeling/regression/um_explain/
ls results/HD_DIT_HAP_generationRAW/22_machine_learning_modeling/regression/lam_explain/
```

Expected files in each: `features_importance.csv`, `prediction_and_residuals.pdf`, `feature_importance_top20.pdf`

Expected in parent `regression/` dir: `um_vs_lam_feature_importance_explain.pdf`
