# Control insertion selection
rule control_insertion_selection:
    input:
        counts_df = rules.hard_filtering.output,
        annotation_df = rules.concat_counts_and_annotations.output.annotations
    output:
        f"results/{project_name}/13_filtered/control_insertions.tsv"
    log:
        f"logs/{project_name}/depletion_analysis/control_insertion_selection.log"
    conda:
        "../envs/tabular_operations.yml"
    message:
        "*** Selecting control insertions..."
    shell:
        """
        python workflow/scripts/depletion_analysis/def_ctr_insertions.py \
        -i {input.counts_df} \
        -a {input.annotation_df} \
        -o {output} &> {log}
        """


# Imputation of the insertions in the coding genes using Forward/Reverse insertions for better coverage in the coding genes
# The forward/reverse insertions tend to have the similar disruption effect on the fitness
# -----------------------------------------------------
rule impute_missing_values_using_FR:
    input:
        filtered_reads = rules.hard_filtering.output,
        annotation = rules.concat_counts_and_annotations.output.annotations
    output:
        f"results/{project_name}/14_imputed_missing_values_using_FR/imputed_raw_reads.tsv"
    log:
        "logs/depletion_analysis/impute_missing_values_using_FR.log"
    conda:
        "../envs/tabular_operations.yml"
    message:
        "*** Imputing missing values using FR..."
    shell:
        """
        python workflow/scripts/depletion_analysis/impute_missing_values_using_FR.py \
        -i {input.filtered_reads} \
        -a {input.annotation} \
        -o {output} &> {log}
        """

# insertion-level depletion analysis
# -----------------------------------------------------
rule insertion_level_depletion_analysis:
    input:
        counts_df = rules.impute_missing_values_using_FR.output,
        annotations_df = rules.concat_counts_and_annotations.output.annotations
    output:
        all_statistics = f"results/{project_name}/15_insertion_level_depletion_analysis/insertions_LFC.tsv",
        LFC = f"results/{project_name}/15_insertion_level_depletion_analysis/LFC.tsv",
        lfcSE = f"results/{project_name}/15_insertion_level_depletion_analysis/lfcSE.tsv",
        padj = f"results/{project_name}/15_insertion_level_depletion_analysis/padj.tsv"
    log:
        "logs/depletion_analysis/insertion_level_depletion_analysis.log"
    params:
        initial_time_point=config["initial_time_point"]
    conda:
        "../envs/pydeseq2.yml"
    message:
        "*** Running insertion-level depletion analysis..."
    shell:
        """
        python workflow/scripts/depletion_analysis/insertion_level_depletion_analysis.py \
        -i {input.counts_df} \
        -a {input.annotations_df} \
        -t {params.initial_time_point} \
        -o {output.all_statistics} &> {log}
        """

# insertion-level curve fitting
# -----------------------------------------------------
rule insertion_level_curve_fitting:
    input:
        LFC = rules.insertion_level_depletion_analysis.output.LFC
    output:
        f"results/{project_name}/16_insertion_level_curve_fitting/insertions_LFC_fitted.tsv"
    log:
        "logs/depletion_analysis/insertion_level_curve_fitting.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Running insertion-level curve fitting..."
    params:
        time_points = " ".join(map(str, config["time_points"]))
    shell:
        """
        python workflow/scripts/depletion_analysis/curve_fitting.py \
        -i {input.LFC} \
        -t {params.time_points} \
        -o {output} &> {log}
        """

# gene-level depletion analysis
# -----------------------------------------------------
rule gene_level_depletion_analysis:
    input:
        lfc_path = rules.insertion_level_depletion_analysis.output.all_statistics,
        annotations_path = rules.concat_counts_and_annotations.output.annotations
    output:
        all_statistics = f"results/{project_name}/17_gene_level_depletion_analysis/Gene_level_statistics.tsv",
        LFC = f"results/{project_name}/17_gene_level_depletion_analysis/LFC.tsv"
    log:
        "logs/depletion_analysis/gene_level_depletion_analysis.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Running gene-level depletion analysis..."
    params:
        method = "fisher"
    shell:
        """
        python workflow/scripts/depletion_analysis/gene_level_depletion_analysis.py \
        -l {input.lfc_path} \
        -a {input.annotations_path} \
        -o {output.all_statistics} \
        -m {params.method} &> {log}
        """

rule gene_level_curve_fitting:
    input:
        LFC = rules.gene_level_depletion_analysis.output.LFC
    output:
        f"results/{project_name}/18_gene_level_curve_fitting/Gene_level_statistics_fitted.tsv"
    log:
        "logs/depletion_analysis/gene_level_curve_fitting.log"
    params:
        time_points = " ".join(map(str, config["time_points"]))
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Running gene-level curve fitting..."
    shell:
        """
        python workflow/scripts/depletion_analysis/curve_fitting.py \
        -i {input.LFC} \
        -t {params.time_points} \
        -o {output} &> {log}
        """