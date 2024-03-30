rule reads_hard_filtering:
    input:
        project_dir/"7_concated/{sample}/{sample}.Reads"
    output:
        project_dir/"8_filtered_reads/{sample}.filtered.Reads"
    params:
        initial_time_point=config["initial_time_point"],
        cutoff=config["first_round_hard_filtering_cutoff"]
    shell:
        """
        python workflow/scripts/annotating_and_normalizing/reads_hard_filtering.py \
        -i {input} \
        -itp {params.initial_time_point} \
        -c {params.cutoff} \
        -o {output}
        """

rule annotate_insertions:
    input:
        insertion=rules.reads_hard_filtering.output,
        genome_region=output_dir / "Genome_regions_CDS_intron_IGR_annotated.bed"
    output:
        annotated_insertions=project_dir/"9_annotated_insertions/{sample}.annotated.csv"
    conda:
        "pybedtools"
    shell:
        """
        python workflow/scripts/annotating_and_normalizing/annotate_insertions.py \
        -i {input.insertion} \
        -g {input.genome_region} \
        -o {output.annotated_insertions}
        """

rule merge_annotation:
    input:
        expand(rules.annotate_insertions.output, sample=samples)
    output:
        project_dir/f"9_annotated_insertions/{config['project_name']}.annotated.csv"
    run:
        import pandas as pd
        from pathlib import Path
        insertions = [pd.read_csv(f, header=0) for f in input]

        merged_annotated_insertions = insertions[0]
        for i in insertions[1:]:
            merged_annotated_insertions = merged_annotated_insertions.merge(i, how="outer")
        
        merged_annotated_insertions.rename(columns={"Start": "Coordinate"}, inplace=True)
        merged_annotated_insertions.drop(columns=["End"], inplace=True)

        merged_annotated_insertions.sort_values(["#Chr", "Coordinate", "Strand", "Target"]).to_csv(output[0], header=True, index=False, float_format="%.3f")
    
rule reformat_insertion:
    input:
        rules.reads_hard_filtering.output
    output:
        project_dir/"10_reformatted_insertions/{sample}.csv"
    shell:
        """
        python workflow/scripts/annotating_and_normalizing/reformat_insertion.py \
        -i {input} \
        -o {output}
        """

checkpoint normalization:
    input:
        insertions = expand(rules.reformat_insertion.output, sample=samples),
        annotations = rules.merge_annotation.output
    output:
        directory(project_dir/"11_normalization/")
    params:
        initial_time_point=config["initial_time_point"]
    shell:
        """
        python workflow/scripts/annotating_and_normalizing/normalization_multi_in_multi_out.py \
        -i {input.insertions} \
        -a {input.annotations} \
        -t {params.initial_time_point} \
        -o {output}
        """

# def FWHM_analysis_and_filtering_input(wildcards):
#     checkpoint_output = checkpoints.normalization.get(**wildcards).output.ouput
#     return f"data/4_second_round_normalization/{wildcards.sample}.normalized.csv"

# checkpoint FWHM_analysis_and_filtering:
#     input:
#         normalized_reads = FWHM_analysis_and_filtering_input,
#     output:
#         FWHM_report = report("reports/FWHM_analysis_and_filtering/FWHM_report/{sample}.FWHM_report.csv"),
#         FWHM_plot = report("reports/FWHM_analysis_and_filtering/FWHM_plot/{sample}.FWHM_plot.pdf"),
#         MA_plot = report("reports/FWHM_analysis_and_filtering/MA_plot/{sample}.MA_plot.pdf"),
#         filtered_reads = "data/5_FWHM_filtered_reads/{sample}.normalized.csv",
#         Ms = "data/5_FWHM_filtered_reads/{sample}.M.csv",
#         As = "data/5_FWHM_filtered_reads/{sample}.A.csv"
#     params:
#         initial_time_point=config["initial_time_point"]
#     shell:
#         """
#         python src/FWHM_analysis_and_filtering.py \
#         -i {input.normalized_reads} \
#         -itp {params.initial_time_point} \
#         -or {output.filtered_reads} \
#         -om {output.Ms} \
#         -oa {output.As} \
#         -r {output.FWHM_report} \
#         -f {output.FWHM_plot} \
#         -ma {output.MA_plot}
#         """


