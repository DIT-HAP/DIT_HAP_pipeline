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
configfile: "config/config_1328_spore2YES6.yaml"
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
        # "resources/pombase_data/{release_version}/genome_region/genome_intervals.bed".format(release_version=config["Pombase_release_version"])
        f"results/{project_name}/13_filtered/control_insertions.tsv"
        # f"results/{project_name}/16_insertion_level_curve_fitting/insertions_LFC_fitted.tsv",
        # f"results/{project_name}/18_gene_level_curve_fitting/Gene_level_statistics_fitted.tsv",
        # f"results/{project_name}/19_insertion_in_non_coding_genes/Non_coding_genes_Gene_level_statistics_fitted.tsv",
        # f"results/{project_name}/12_concatenated/annotations.tsv"
        # directory(f"reports/{project_name}/gene_coverage_analysis")

# load rules
# -----------------------------------------------------
include: "workflow/rules/preparation.smk"
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/depletion_analysis.smk"
# include: "workflow/rules/quality_control.smk"
# include: "workflow/rules/misc.smk"


