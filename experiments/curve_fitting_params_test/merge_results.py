"""
Merge the results of the curve fitting parameters test.
"""

# =============================== Imports ===============================
from pathlib import Path
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from loguru import logger
# =============================== Constants ===============================
RESULT_DIR = Path("results")

plt.style.use("/data/c/yangyusheng_optimized/DIT_HAP_pipeline/config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# =============================== Setup Logging ===============================
def setup_logging():
    """Setup logging."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level="INFO", colorize=False)
    return logger

# =============================== Functions ===============================
@logger.catch
def parse_files(dir_path: Path) -> pd.DataFrame:
    """Parse the files in the directory."""
    sub_dirs = list(dir_path.glob("*"))
    logger.info(f"Found {len(sub_dirs)} subdirectories.")
    file_dataframes = pd.DataFrame()
    for idx, sub_dir in enumerate(sub_dirs):
        reg_pattern = re.compile(r"lam_penalty_(\S+)_tol_(\S+)_(Previous|Now)")
        match = reg_pattern.search(sub_dir.name)
        if match:
            logger.info(f"Found match for {sub_dir.name}")
            file_dataframes.loc[idx, "lam_penalty"] = float(match.group(1))
            file_dataframes.loc[idx, "tol"] = float(match.group(2))
            file_dataframes.loc[idx, "input_file_label"] = match.group(3)
            file_dataframes.loc[idx, "input_file"] = str(sub_dir/"gene_level_fitting_statistics.tsv")
    return file_dataframes

@logger.catch
def plot_results(file_dataframes: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    """Plot the results."""
    with PdfPages(output_path) as pdf:
        for idx in file_dataframes.index:
            input_file_label = file_dataframes.loc[idx, "input_file_label"]
            lam_penalty = file_dataframes.loc[idx, "lam_penalty"]
            tol = file_dataframes.loc[idx, "tol"]
            fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
            df = pd.read_csv(file_dataframes.loc[idx, "input_file"], sep="\t")
            success_rate = df["Status"].value_counts()["Success"] / len(df)
            file_dataframes.loc[idx, "success_rate"] = success_rate
            failed_num = sum(df["Status"]!="Success")
            file_dataframes.loc[idx, "failed_num"] = failed_num
            logger.info(f"input_file_label: {input_file_label}, lam_penalty: {lam_penalty}, tol: {tol}\nsuccess_rate: {success_rate:.2f}, failed_num: {failed_num}")
            df["um"].hist(bins=50, ax=axes[0], alpha=0.7)
            axes[0].set_title("um")
            df["lam"].hist(bins=50, ax=axes[1], alpha=0.7)
            axes[1].set_title("lam")
            plt.title(f"input_file_label: {input_file_label}, lam_penalty: {lam_penalty}, tol: {tol}\nsuccess_rate: {success_rate:.2f}, failed_num: {failed_num}")
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    return file_dataframes
# =============================== Main ===============================
@logger.catch
def main():
    """Main function."""
    logger = setup_logging()
    logger.info("Merging the results of the curve fitting parameters test.")
    file_dataframes = parse_files(RESULT_DIR)
    file_dataframes = plot_results(file_dataframes, RESULT_DIR/"results.pdf")
    file_dataframes.to_csv(RESULT_DIR/"results.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()