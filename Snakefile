# Main entrypoint of the workflow. 
# Please follow the best practices: 
# https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html,
# in particular regarding the standardized folder structure mentioned there. 

from snakemake.utils import min_version
from pathlib import Path
import pandas as pd

min_version("8.0")

# load configuration
# -----------------------------------------------------
snakemake_config_file="config/config_HD_generationPLUS1.yaml"
# snakemake_config_file="config/config_LD_generationPLUS1.yaml"
# snakemake_config_file="config/config_1328_spore2YES6.yaml"
# snakemake_config_file="config/config_LD_haploid.yaml"
# snakemake_config_file="config/config_HD_diploid.yaml"
configfile: snakemake_config_file
workdir: "/data/c/yangyusheng_optimized/DIT_HAP_pipeline"

# optional messages, log and error handling
# -----------------------------------------------------
onstart:
    print("\n--- Analysis started ---\n")


onsuccess:
    print("\n--- Workflow finished! ---\n")


onerror:
    print("\n--- An error occurred! ---\n")

# load project name and project directory
# -----------------------------------------------------
project_name = config["project_name"]

# load sample sheet and rename the fastq files, and update the sample sheet
# -----------------------------------------------------
sample_sheet = pd.read_csv(config["sample_sheet"], sep="\t")
samples = sample_sheet["Sample"].unique().tolist()
timepoints = sample_sheet["Timepoint"].unique().tolist()
conditions = sample_sheet["Condition"].unique().tolist()

sample_sheet_dict = { s: {t: {c: {"fq1": None, "fq2": None} for c in conditions} for t in timepoints} for s in samples }

for index, row in sample_sheet.iterrows():
    # Original source paths from the sample sheet
    source_fq1_path = Path(row["read1"]) # Original file path
    source_fq2_path = Path(row["read2"]) # Original file path

    sample_sheet_dict[row["Sample"]][row["Timepoint"]][row["Condition"]]["fq1"] = source_fq1_path
    sample_sheet_dict[row["Sample"]][row["Timepoint"]][row["Condition"]]["fq2"] = source_fq2_path

# wildcard constraint
# -----------------------------------------------------
wildcard_constraints:
    sample = "|".join(samples),
    timepoint = "|".join(timepoints),
    condition = "|".join(conditions)

# target rules
# -----------------------------------------------------
rule all:
    input:
        f"results/{project_name}/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
        f"results/{project_name}/17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
        f"reports/{project_name}/multiqc/quality_control_multiqc_report.html",
        f"reports/{project_name}/mapping_filtering_statistics/mapping_filtering_statistics.tsv",
        f"reports/{project_name}/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
        f"reports/{project_name}/read_count_distribution_analysis/read_count_distribution_analysis.pdf",
        f"reports/{project_name}/insertion_orientation_analysis/insertion_orientation_analysis.pdf",
        f"reports/{project_name}/insertion_density_analysis/insertion_density_analysis.tsv",
        f"reports/{project_name}/depletion_LFC_and_curve_features_analysis/insertion_level_depletion_LFC_and_curve_features_analysis.pdf",
        f"reports/{project_name}/depletion_LFC_and_curve_features_analysis/gene_level_depletion_and_curve_features_analysis.pdf",
        f"reports/{project_name}/gene_coverage_analysis"

# load rules
# -----------------------------------------------------
include: "workflow/rules/preparation.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/depletion_analysis.smk"
include: "workflow/rules/quality_control.smk"
# include: "workflow/rules/misc.smk"


