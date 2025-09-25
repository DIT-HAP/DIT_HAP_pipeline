"""
Test the effect of lam_penalty and tol on the curve fitting results.
"""

# =============================== Imports ===============================
from pathlib import Path
import sys
import importlib.util
import pandas as pd
import matplotlib.pyplot as plt
import mlflow
import argparse
from loguru import logger
from multiprocessing import Pool
from types import ModuleType
from typing import List

# =============================== Constants ===============================
TARGET_SCRIPT_PATH = "../../workflow/scripts/depletion_analysis/curve_fitting.py"
NUM_PROCESSES = 32

plt.style.use("/data/c/yangyusheng_optimized/DIT_HAP_pipeline/config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# =============================== Configuration ===============================
# lam_penalty_values = [1e-3, 6e-3, 1e-2, 5e-2]
# tol_values = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4]
lam_penalty_values = [5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2]
tol_values = [1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 1e-5]
input_files = {
    "Previous": Path("../../results/HD_DIT_HAP/17_gene_level_depletion_analysis/LFC.tsv"),
    "Now": Path("../../results/HD_DIT_HAP_generationPLUS1/16_gene_level_depletion_analysis/LFC.tsv"),
}

time_points = {
    "Previous": [0, 2.352, 5.588, 9.104, 12.480],
    "Now": [0, 3.352, 6.588, 10.104, 13.480]
}

# setting up mlflow
mlflow.set_tracking_uri(f"file://{(Path(__file__).parent / 'mlruns').resolve()}")
mlflow.set_experiment("lam_penalty_tol_sweep") 

# =============================== Setup Logging ===============================
def setup_logging():
    """Setup logging."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level="INFO", colorize=True)
    return logger

# =============================== Functions ===============================
@logger.catch
def load_script_as_module(script_path: Path) -> ModuleType:
    """Load a Python script as a module."""
    module_name = Path(script_path).stem
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    module = importlib.util.module_from_spec(spec)
    # add the module path to sys.path, to handle relative imports in the script
    sys.path.insert(0, Path(script_path).parent)
    spec.loader.exec_module(module)
    sys.path.pop(0) # restore sys.path
    return module

@logger.catch
def run_with_params(lam_penalty: float, tol: float, target_script_module: ModuleType, input_file: Path, input_file_label: str, time_points: List[float]) -> None:
    """Run the script with the given parameters."""
    run_name = f"lam_penalty_{lam_penalty}_tol_{tol}_{input_file_label}"
    with mlflow.start_run(run_name=run_name):
        logger.info(f"--- Running with LAM_PENALTY={lam_penalty}, TOL={tol}, INPUT_FILE={input_file_label} ---")
        
        # log parameters
        mlflow.log_param("LAM_PENALTY", lam_penalty)
        mlflow.log_param("TOL", tol)

        # create the output directory for this run
        output_dir = Path("results")/run_name
        output_dir.mkdir(parents=True, exist_ok=True)
        output_tsv_path = output_dir/"gene_level_fitting_statistics.tsv"
        
        original_parser = target_script_module.parse_arguments

        # define a fake, replacement function
        def fake_parser():
            # create a mock Namespace object, like the one returned by argparse
            args = argparse.Namespace(
                input=input_file,
                time_points=time_points,
                output=output_tsv_path,
                weight=None,
                verbose=False
            )
            return args

        # --- core steps: modify the constants in the loaded module ---
        setattr(target_script_module, 'parse_arguments', fake_parser)
        setattr(target_script_module, 'LAM_PENALTY', lam_penalty)
        setattr(target_script_module, 'TOL', tol)

        try:
            target_script_module.main()
        except Exception as e:
            logger.error(f"Error: executing the main function of the script: {e}")
            mlflow.log_param("error", str(e))
        finally:
            setattr(target_script_module, 'parse_arguments', original_parser)
            setattr(target_script_module, 'LAM_PENALTY', 6e-3)
            setattr(target_script_module, 'TOL', 1e-6)
    
        # analyze and plot (same as before)
        try:
            df = pd.read_csv(output_tsv_path, sep='\t')

            success_rate = df["Status"].value_counts()["Success"] / len(df)
            mlflow.log_metric("success_rate", success_rate)

            fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
            df["um"].hist(bins=50, ax=axes[0], alpha=0.7)
            axes[0].set_title("um")
            df["lam"].hist(bins=50, ax=axes[1], alpha=0.7)
            axes[1].set_title("lam")
            plt.tight_layout()
            plt.savefig(output_dir/"hist_um_lam.png")
            plt.close()
            
            # log all artifacts
            mlflow.log_artifacts(output_dir)
            logger.info("Logged artifacts to MLflow.")
        except Exception as e:
            logger.error(f"Error: analyzing the results: {e}")
            mlflow.log_param("analysis_error", str(e))

# =============================== Main ===============================
@logger.catch
def main():
    """Main function."""
    logger = setup_logging()
    try:
        target_script_module = load_script_as_module(TARGET_SCRIPT_PATH)
    except FileNotFoundError:
        logger.error(f"Error: target script '{TARGET_SCRIPT_PATH}' not found. Please check the path.")
        exit(1)

    # iterate over all parameter combinations
    for lam_penalty in lam_penalty_values:
        for tol in tol_values:
            for input_file_label, input_file in input_files.items():
                run_with_params(lam_penalty, tol, target_script_module, input_file, input_file_label, time_points[input_file_label])

    logger.info("\nAll parameter combinations tested!")
    logger.info("Please run 'mlflow ui' in the 'experiments' folder to view the results.")


if __name__ == "__main__":
    main()