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

# distinguish the replicates and non-replicates
# -----------------------------------------------------
if True:
    # insertion-level depletion analysis
    # -----------------------------------------------------
    rule insertion_level_depletion_analysis_no_replicates:
        input:
            counts_df = rules.hard_filtering.output,
            control_insertions_df = rules.control_insertion_selection.output
        output:
            all_statistics = f"results/{project_name}/14_insertion_level_depletion_analysis/insertions_LFC.tsv",
            LFC = f"results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv"
        log:
            f"logs/{project_name}/depletion_analysis/insertion_level_depletion_analysis_no_replicates.log"
        params:
            initial_time_point=config["initial_time_point"]
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** Running insertion-level depletion analysis..."
        shell:
            """
            python workflow/scripts/depletion_analysis/insertion_level_depletion_analysis_no_replicates.py \
            -i {input.counts_df} \
            -c {input.control_insertions_df} \
            -t {params.initial_time_point} \
            -a {output.all_statistics} \
            -l {output.LFC} &> {log}
            """

else:
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
            f"logs/{project_name}/depletion_analysis/impute_missing_values_using_FR.log"
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
    rule insertion_level_depletion_analysis_has_replicates:
        input:
            counts_df = rules.impute_missing_values_using_FR.output,
            annotations_df = rules.concat_counts_and_annotations.output.annotations
        output:
            all_statistics = f"results/{project_name}/14_insertion_level_depletion_analysis/insertions_LFC.tsv",
            LFC = f"results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
            lfcSE = f"results/{project_name}/14_insertion_level_depletion_analysis/lfcSE.tsv",
            padj = f"results/{project_name}/14_insertion_level_depletion_analysis/padj.tsv"
        log:
            f"logs/{project_name}/depletion_analysis/insertion_level_depletion_analysis_has_replicates.log"
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
        branch(
            config["has_replicates"],
            then=f"results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
            otherwise=f"results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv"
        )
    output:
        f"results/{project_name}/15_insertion_level_curve_fitting/insertions_LFC_fitted.tsv"
    log:
        f"logs/{project_name}/depletion_analysis/insertion_level_curve_fitting.log"
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

if True:
    rule r_square_as_weights:
        input:
            rules.insertion_level_curve_fitting.output
        output:
            f"results/{project_name}/15_insertion_level_curve_fitting/insertions_LFC_fitted_with_r_square_as_weights.tsv"
        log:
            f"logs/{project_name}/depletion_analysis/r_square_as_weights.log"
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** Running r-square as weights..."
        run:
            import pandas as pd
            fitting_res = pd.read_csv(input[0], sep="\t", index_col=[0,1,2,3])

            fitting_res["confidence"] = fitting_res["R2"].clip(lower=1e-3, upper=1-1e-3)

            log2FoldChange = fitting_res.filter(regex=r"^t(\d+)$").astype(float)
            log2FoldChange.columns = log2FoldChange.columns.str.replace("t", "YES")

            padj = pd.DataFrame(index=fitting_res.index, columns=log2FoldChange.columns)

            for col in log2FoldChange.columns:
                padj[col] = 1 - fitting_res["confidence"]

            concated_statistics = pd.concat([log2FoldChange, padj], axis=1, keys=["log2FoldChange", "padj"])
            concated_statistics = concated_statistics.rename_axis(["Statistic", "Timepoint"], axis=1)
            concated_statistics.reorder_levels(["Timepoint", "Statistic"], axis=1).to_csv(output[0], sep="\t")

# gene-level depletion analysis
# -----------------------------------------------------
rule gene_level_depletion_analysis:
    input:
        lfc_path = branch(
            config["has_replicates"],
            then=f"results/{project_name}/14_insertion_level_depletion_analysis/insertions_LFC.tsv",
            otherwise=f"results/{project_name}/15_insertion_level_curve_fitting/insertions_LFC_fitted_with_r_square_as_weights.tsv"
        ),
        annotations_path = rules.concat_counts_and_annotations.output.annotations
    output:
        all_statistics = f"results/{project_name}/17_gene_level_depletion_analysis/Gene_level_statistics.tsv",
        LFC = f"results/{project_name}/17_gene_level_depletion_analysis/LFC.tsv"
    log:
        f"logs/{project_name}/depletion_analysis/gene_level_depletion_analysis.log"
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

# rule gene_level_curve_fitting:
#     input:
#         LFC = rules.gene_level_depletion_analysis.output.LFC
#     output:
#         f"results/{project_name}/18_gene_level_curve_fitting/Gene_level_statistics_fitted.tsv"
#     log:
#         "logs/depletion_analysis/gene_level_curve_fitting.log"
#     params:
#         time_points = " ".join(map(str, config["time_points"]))
#     conda:
#         "../envs/statistics_and_figure_plotting.yml"
#     message:
#         "*** Running gene-level curve fitting..."
#     shell:
#         """
#         python workflow/scripts/depletion_analysis/curve_fitting.py \
#         -i {input.LFC} \
#         -t {params.time_points} \
#         -o {output} &> {log}
#         """