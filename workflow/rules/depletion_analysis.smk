# Imputation of the insertions in the coding genes using Forward/Reverse insertions for better coverage in the coding genes
# The forward/reverse insertions tend to have the similar disruption effect on the fitness
# -----------------------------------------------------
rule impute_missing_values_using_FR:
    input:
        raw_reads = rules.concat_counts_and_annotations.output.counts,
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
        -i {input.raw_reads} \
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
        all_statistics = f"results/{project_name}/15_insertion_level_depletion_analysis/insertions_LFC.csv",
        LFC = f"results/{project_name}/15_insertion_level_depletion_analysis/LFC.csv",
        padj = f"results/{project_name}/15_insertion_level_depletion_analysis/padj.csv"
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

rule insertion_level_curve_fitting:
    input:
        rules.insertion_level_depletion_analysis.output.LFC,
    output:
        f"results/{project_name}/16_insertion_level_curve_fitting/insertions_LFC_fitted.csv"
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
        -i {input} \
        -t {params.time_points} \
        -o {output} &> {log}
        """

# 

# rule calculate_gene_level_statistics:
#     input:
#         lfc_path = rules.DEseq2_analysis.output,
#         annotations_path = rules.merge_annotation.output
#     output:
#         project_dir/"13_gene_level_statistics/GWMs.csv"
#     shell:
#         """
#         python workflow/scripts/annotating_and_DEseq2_analysis/calculate_gene_level_statistic.py \
#         -l {input.lfc_path} \
#         -a {input.annotations_path} \
#         -o {output}
#         """