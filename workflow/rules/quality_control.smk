# MultiQC for preprocessing reports
# -----------------------------------------------------
rule multiqc_preprocessing:
    input:
        fastp_json=expand(rules.fastp_preprocessing.output.json, sample=samples, timepoint=timepoints, condition=conditions),
        demux_json=expand(rules.demultiplexing.output.json, sample=samples, timepoint=timepoints, condition=conditions),
        fastqc=expand(f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.{{read}}_fastqc.zip", sample=samples, timepoint=timepoints, condition=conditions, read=["PBL_1", "PBL_2", "PBR_1", "PBR_2"]),
        bam_stats=expand(f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.{{type}}.txt", sample=samples, timepoint=timepoints, condition=conditions, fragment=["PBL", "PBR"], type=["stats", "flagstat", "idxstats"])
    output:
        report(f"reports/{project_name}/multiqc/{project_name}_quality_control_multiqc_report.html")
    log:
        f"logs/quality_control/multiqc_quality_control.log"
    conda:
        "../envs/multiqc.yml"
    params:
        outdir=f"reports/{project_name}/multiqc",
        title=f"Quality Control Report for {project_name}"
        # cl_config=config["multiqc_config"]
    message:
        "*** Generating MultiQC report for quality control..."
    shell:
        """
        multiqc \
            --title "{params.title}" \
            --filename {project_name}_quality_control_multiqc_report.html \
            --outdir {params.outdir} \
            --force \
            --verbose \
            --fn_as_s_name \
            {input.fastp_json} \
            {input.demux_json} \
            {input.fastqc} \
            {input.bam_stats} \
            &> {log}
        """

# Mapping filtering statistics
# -----------------------------------------------------
rule mapping_filtering_statistics:
    input:
        expand(rules.read_pair_filtering.log, sample=samples, timepoint=timepoints, condition=conditions)
    output:
        report(f"reports/{project_name}/mapping_filtering_statistics/{project_name}_mapping_filtering_statistics.tsv")
    log:
        f"logs/quality_control/mapping_filtering_statistics.log"
    conda:
        "../envs/tabular_operations.yml"
    message:
        "*** Generating mapping filtering statistics..."
    shell:
        """
        python workflow/scripts/quality_control/extract_mapping_filtering_statistics.py \
        -i {input} \
        -o {output} &> {log}
        """
        
# PBL-PBR correlation analysis
# -----------------------------------------------------
rule PBL_PBR_correlation_analysis:
    input:
        expand(rules.merge_insertions.output, sample=samples, timepoint=timepoints, condition=conditions)
    output:
        report(f"reports/{project_name}/PBL_PBR_correlation_analysis/{project_name}_PBL_PBR_correlation_analysis.pdf")
    log:
        f"logs/quality_control/PBL_PBR_correlation_analysis.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing PBL-PBR correlation analysis..."
    shell:
        """
        python workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py -i {input} -o {output} &> {log}
        """

# Read count distribution analysis
# -----------------------------------------------------
rule read_count_distribution_analysis:
    input:
        expand(rules.merge_similar_timepoints.output, sample=samples, condition=conditions)
    output:
        report(f"reports/{project_name}/read_count_distribution_analysis/{project_name}_read_count_distribution_analysis.pdf")
    log:
        f"logs/quality_control/read_count_distribution_analysis.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing read count distribution analysis..."
    params:
        initial_time_point = config["initial_time_point"],
        hard_filtering_cutoff = config["hard_filtering_cutoff"]
    shell:
        """
        python workflow/scripts/quality_control/read_count_distribution_analysis.py -i {input} -t {params.initial_time_point} -c {params.hard_filtering_cutoff} -o {output} &> {log}
        """

# Insertion orientation analysis
# -----------------------------------------------------
rule insertion_orientation_analysis:
    input:
        expand(rules.hard_filtering.output, sample=samples, condition=conditions)
    output:
        report(f"reports/{project_name}/insertion_orientation_analysis/{project_name}_insertion_orientation_analysis.pdf")
    log:
        f"logs/quality_control/insertion_orientation_analysis.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing insertion orientation analysis..."
    shell:
        """
        python workflow/scripts/quality_control/insertion_orientation_analysis.py -i {input} -o {output} &> {log}
        """