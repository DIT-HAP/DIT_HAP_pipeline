rule concat_samples:
    input:
        Ms = expand(project_dir/"12_FWHM_filtered_reads/{sample}.M.csv", sample=samples),
        As = expand(project_dir/"12_FWHM_filtered_reads/{sample}.A.csv", sample=samples),
        reads = expand(project_dir/"12_FWHM_filtered_reads/{sample}.normalized.csv", sample=samples)
    output:
        normalized=project_dir/f"13_concated_values/{project_name}.normalized.csv",
        M_value=project_dir/f"13_concated_values/{project_name}.M.csv",
        A_value=project_dir/f"13_concated_values/{project_name}.A.csv",
        Confidence_score=project_dir/f"13_concated_values/{project_name}.Confidence_score.csv"
    params:
        initial_timepoint = config["initial_time_point"],
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

rule data_smoothing:
    input:
        M_value=rules.concat_samples.output.M_value,
        Confidence_score=rules.concat_samples.output.Confidence_score,
        generation = config["samples_with_timepoints"]
    output:
        before_smoothing=project_dir/f"14_smoothed_values/{project_name}.before_smoothing.csv",
        after_smoothing=project_dir/f"14_smoothed_values/{project_name}.after_smoothing.csv"
    params:
        initial_timepoint=config["initial_time_point"]
    threads: 12
    shell:
        """
        python workflow/scripts/depletion_curve_smoothing/data_smoothing.py \
            -M {input.M_value} \
            -CS {input.Confidence_score} \
            -g {input.generation} \
            -C {threads} \
            -itp {params.initial_timepoint} \
            -ob {output.before_smoothing} \
            -oa {output.after_smoothing}
        """

rule site_level_curve_fitting:
    input:
        project_dir/f"14_smoothed_values/{project_name}.{{smooth}}.csv"
    output:
        project_dir/f"15_curve_fitting/{project_name}.{{smooth}}.site.csv"
    params:
        fatol = 5e-6
    threads: 12
    shell:
        """
        python workflow/scripts/depletion_curve_smoothing/site_level_curve_fitting.py \
            -GM {input} \
            -C {threads} \
            -fatol {params.fatol} \
            -o {output}
        """

rule domain_level_curve_fitting:
    input:
        annotations = project_dir/f"9_annotated_insertions/{project_name}.annotated.csv",
        GMs = project_dir/f"14_smoothed_values/{project_name}.{{smooth}}.csv",
        domains = config["domain_file"]
    output:
        project_dir/f"15_curve_fitting/{project_name}.{{smooth}}.domain.csv"
    params:
        fatol = 5e-3
    shell:
        """
        python workflow/scripts/depletion_curve_smoothing/domain_level_curve_fitting.py \
            -anno {input.annotations} \
            -GM {input.GMs} \
            -d {input.domains} \
            -o {output} \
            -fatol {params.fatol}
        """

rule gene_level_curve_fitting:
    input:
        annotations = project_dir/f"9_annotated_insertions/{project_name}.annotated.csv",
        GMs = project_dir/f"14_smoothed_values/{project_name}.{{smooth}}.csv",
        domains = config["domain_file"],
        DDR = rules.domain_level_curve_fitting.output
    output:
        project_dir/f"15_curve_fitting/{project_name}.{{smooth}}.gene.csv"
    params:
        fatol = 5e-3
    shell:
        """
        python workflow/scripts/depletion_curve_smoothing/gene_level_curve_fitting.py \
            -anno {input.annotations} \
            -GM {input.GMs} \
            -d {input.domains} \
            -ddr {input.DDR} \
            -o {output} \
            -fatol {params.fatol}
        """