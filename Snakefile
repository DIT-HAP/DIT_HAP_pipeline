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
configfile: "config/config.yaml"
workdir: "/data/c/yangyusheng_optimized/DIT_HAP_pipeline"
multiqc_config: "config/multiqc_config.yml"

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
file_to_rename_dict = {}
for index, row in sample_sheet.iterrows():
    # Original source paths from the sample sheet
    source_fq1_path = Path(row["read1"]) # Original file path
    source_fq2_path = Path(row["read2"]) # Original file path
    # This is the base name for the new symlinks and the key for the dictionary
    renamed_name_key = row["Sample"] + "_" + row["Timepoint"] + "_" + row["Condition"]
    # Using renamed_name_key which matches the original intent for the key.
    file_to_rename_dict[renamed_name_key] = [source_fq1_path, source_fq2_path]

# wildcard constraint
# -----------------------------------------------------
wildcard_constraints:
    sample = "|".join(file_to_rename_dict.keys())


# target rules
# -----------------------------------------------------
rule all:
    input:
        # expand(f"results/{project_name}/5_sorted/{{sample}}.PBL.sorted.bam", sample=list(file_to_rename_dict.keys())),
        # expand(f"results/{project_name}/5_sorted/{{sample}}.PBR.sorted.bam", sample=list(file_to_rename_dict.keys())),
        # f"resources/pombase_data/{config['Pombase_release_version']}/fasta/Schizosaccharomyces_pombe_all_chromosomes.fa",
        # expand(f"results/{project_name}/4_filtered/{{sample_id}}.{{name}}.unique.bam.bai", sample_id=list(file_to_rename_dict.keys()), name=["PBL", "PBR"]),
        # report(f"reports/{project_name}/multiqc/preprocessing_multiqc_report.html")
        # directory(f"reports/{project_name}/fastqc")
        expand(f"reports/{project_name}/samtools_mapping_statistics_after_filtering/{{sample}}.PBL.stats.txt", sample=list(file_to_rename_dict.keys())),
        expand(f"reports/{project_name}/samtools_mapping_statistics_after_filtering/{{sample}}.PBR.stats.txt", sample=list(file_to_rename_dict.keys()))

# load rules
# -----------------------------------------------------
include: "workflow/rules/preparation.smk"
include: "workflow/rules/preprocessing.smk"


