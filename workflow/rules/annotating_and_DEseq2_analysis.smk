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
        python workflow/scripts/annotating_and_DEseq2_analysis/reads_hard_filtering.py \
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
        python workflow/scripts/annotating_and_DEseq2_analysis/annotate_insertions.py \
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
        python workflow/scripts/annotating_and_DEseq2_analysis/reformat_insertion.py \
        -i {input} \
        -o {output}
        """

rule impute_missing_values_using_FR:
    input:
        raw_reads = expand(rules.reformat_insertion.output, sample=samples),
        annotation = rules.merge_annotation.output
    output:
        project_dir/"11_imputed_missing_values_using_FR/imputed_raw_reads.csv"
    shell:
        """
        python workflow/scripts/annotating_and_DEseq2_analysis/impute_missing_values_using_FR.py \
        -i {input.raw_reads} \
        -a {input.annotation} \
        -o {output}
        """

rule DEseq2_analysis:
    input:
        counts_df = rules.impute_missing_values_using_FR.output,
        annotations_df = rules.merge_annotation.output
    output:
        project_dir/"12_DEseq2_analysis/insertions_LFC.csv"
    params:
        initial_time_point=config["initial_time_point"]
    conda:
        "pydeseq2"
    shell:
        """
        python workflow/scripts/annotating_and_DEseq2_analysis/DEseq2_analysis.py \
        -i {input.counts_df} \
        -a {input.annotations_df} \
        -o {output} \
        -t {params.initial_time_point}
        """

rule calculate_gene_level_statistics:
    input:
        lfc_path = rules.DEseq2_analysis.output,
        annotations_path = rules.merge_annotation.output
    output:
        project_dir/"13_gene_level_statistics/GWMs.csv"
    shell:
        """
        python workflow/scripts/annotating_and_DEseq2_analysis/calculate_gene_level_statistic.py \
        -l {input.lfc_path} \
        -a {input.annotations_path} \
        -o {output}
        """