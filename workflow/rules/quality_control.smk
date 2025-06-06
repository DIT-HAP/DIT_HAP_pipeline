# Quality Control Rules for Transposon Insertion Sequencing
# ==========================================================
#
# This file contains rules for comprehensive quality control analysis
# of transposon insertion sequencing data. Each rule generates reports
# with enhanced metadata for beautiful HTML report generation.
#
# To generate the comprehensive HTML report, run:
#   snakemake --report reports/{project_name}/{project_name}_quality_control_report.html
#
# The report includes:
# - MultiQC preprocessing summary
# - Mapping and filtering statistics  
# - PBL-PBR correlation analysis
# - Read count distribution analysis
# - Insertion orientation analysis
# - Insertion density analysis
#
# All reports include detailed captions, categorization, and metadata
# for professional presentation and easy interpretation.

# MultiQC for preprocessing reports
# -----------------------------------------------------
rule multiqc_preprocessing:
    input:
        fastp_json=expand(rules.fastp_preprocessing.output.json, sample=samples, timepoint=timepoints, condition=conditions),
        demux_json=expand(rules.demultiplexing.output.json, sample=samples, timepoint=timepoints, condition=conditions),
        fastqc=expand(f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.{{read}}_fastqc.zip", sample=samples, timepoint=timepoints, condition=conditions, read=["PBL_1", "PBL_2", "PBR_1", "PBR_2"]),
        bam_stats=expand(f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.{{type}}.txt", sample=samples, timepoint=timepoints, condition=conditions, fragment=["PBL", "PBR"], type=["stats", "flagstat", "idxstats"])
    output:
        report(
            f"reports/{project_name}/multiqc/{project_name}_quality_control_multiqc_report.html",
            caption="../report/quality_control_multiqc.rst",
            category="Quality Control",
            subcategory="Preprocessing Summary",
            labels={
                "type": "MultiQC Report",
                "stage": "Preprocessing",
                "format": "HTML"
            }
        )
    log:
        f"logs/quality_control/multiqc_quality_control.log"
    conda:
        "../envs/multiqc.yml"
    params:
        outdir=f"reports/{project_name}/multiqc",
        title=f"Quality Control Report for {project_name}",
        multiqc_config=config["multiqc_config"]
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
            --config {params.multiqc_config} \
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
        report(
            f"reports/{project_name}/mapping_filtering_statistics/{project_name}_mapping_filtering_statistics.tsv",
            caption="../report/mapping_filtering_statistics.rst",
            category="Quality Control",
            subcategory="Mapping Statistics",
            labels={
                "type": "Statistics Table",
                "stage": "Mapping",
                "format": "TSV"
            }
        )
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
        report(
            f"reports/{project_name}/PBL_PBR_correlation_analysis/{project_name}_PBL_PBR_correlation_analysis.pdf",
            caption="../report/pbl_pbr_correlation.rst",
            category="Quality Control",
            subcategory="Correlation Analysis",
            labels={
                "type": "Correlation Plot",
                "stage": "Insertion Analysis",
                "format": "PDF",
                "comparison": "PBL vs PBR"
            }
        )
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
        report(
            f"reports/{project_name}/read_count_distribution_analysis/{project_name}_read_count_distribution_analysis.pdf",
            caption="../report/read_count_distribution.rst",
            category="Quality Control",
            subcategory="Distribution Analysis",
            labels={
                "type": "Distribution Plot",
                "stage": "Read Analysis",
                "format": "PDF",
                "metric": "Read Counts"
            }
        )
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
        report(
            f"reports/{project_name}/insertion_orientation_analysis/{project_name}_insertion_orientation_analysis.pdf",
            caption="../report/insertion_orientation.rst",
            category="Quality Control",
            subcategory="Insertion Analysis",
            labels={
                "type": "Orientation Plot",
                "stage": "Insertion Analysis",
                "format": "PDF",
                "metric": "Insertion Orientation"
            }
        )
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

# Insertion density analysis
# -----------------------------------------------------
rule insertion_density_analysis:
    input:
        insertion_data = rules.insertion_level_depletion_analysis.output.all_statistics,
        annotation = rules.concat_counts_and_annotations.output.annotations
    output:
        table = report(
            f"reports/{project_name}/insertion_density_analysis/{project_name}_insertion_density_analysis.csv",
            caption="../report/insertion_density_analysis.rst",
            category="Quality Control",
            subcategory="Density Analysis",
            labels={
                "type": "Statistics Table",
                "stage": "Insertion Analysis",
                "format": "CSV",
                "metric": "Insertion Density"
            }
        ),
        plot = report(
            f"reports/{project_name}/insertion_density_analysis/{project_name}_insertion_density_analysis_histograms.pdf",
            caption="../report/insertion_density_analysis.rst",
            category="Quality Control",
            subcategory="Density Analysis",
            labels={
                "type": "Density Plot",
                "stage": "Insertion Analysis",
                "format": "PDF",
                "metric": "Insertion Density"
            }
        )
    log:
        f"logs/quality_control/insertion_density_analysis.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing insertion density analysis..."
    params:
        initial_time_point = config["initial_time_point"]
    shell:
        """
        python workflow/scripts/quality_control/insertion_density_analysis.py -i {input.insertion_data} \
                                                                              -a {input.annotation} \
                                                                              -t {params.initial_time_point} \
                                                                              -o {output.table} &> {log}
        """
