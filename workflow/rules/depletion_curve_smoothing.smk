rule concat_samples:
    input:
        Ms = expand(project_dir/"11_normalization/M_values/{sample}.M.csv", sample=samples),
        As = expand(project_dir/"11_normalization/A_values/{sample}.A.csv", sample=samples),
        reads = expand(project_dir/"11_normalization/normalized_values/{sample}.normalized.csv", sample=samples)
    output:
        normalized=project_dir/f"12_concated_values/{project_name}.normalized.csv",
        M_value=project_dir/f"12_concated_values/{project_name}.M.csv",
        A_value=project_dir/f"12_concated_values/{project_name}.A.csv",
        Confidence_score=project_dir/f"12_concated_values/{project_name}.Confidence_score.csv"
    params:
        initial_timepoint = config["initial_timepoint"],
        unused_timepoints = config["unused_timepoints"]
    shell:
        """
        python workflow/scripts/depletion_curve_smoothing/concat_samples.py \
        -M {input.Ms} \
        -A {input.As} \
        -N {input.reads} \
        -uTP {params.unused_timepoints} \
        -itp {params.initial_timepoint} \
        -oM {output.M_value} \
        -oA {output.A_value} \
        -oN {output.normalized} \
        -oCS {output.Confidence_score}
        """

# rule calculate_Site_level_weighted_M:
#     input:
#         M_value=rules.concat_samples.output.M_value,
#         Confidence_score=rules.concat_samples.output.Confidence_score
#     output:
#         Weighted_M="data/2_weighted_M/Site/{Group}.Site.csv",
#         Statistic="data/2_weighted_M/Site/{Group}.Site_statistic.csv"
#     params:
#         initial_timepoint=config["initial_timepoint"]
#     shell:
#         """
#         python src/site_level_weighted_M.py \
#             -M {input.M_value} \
#             -CS {input.Confidence_score} \
#             -itp {params.initial_timepoint} \
#             -o {output.Weighted_M} \
#             -os {output.Statistic}
#         """