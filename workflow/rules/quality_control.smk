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
            f"reports/{project_name}/multiqc/quality_control_multiqc_report.html",
            category="Quality Control",
            subcategory="Preprocessing Summary",
            labels={
                "type": "MultiQC Report",
                "stage": "Preprocessing",
                "format": "HTML"
            }
        )
    log:
        f"logs/{project_name}/quality_control/multiqc_quality_control.log"
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
            --filename quality_control_multiqc_report.html \
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
        expand(rules.filter_aligned_reads.log, sample=samples, timepoint=timepoints, condition=conditions)
    output:
        report(
            f"reports/{project_name}/mapping_filtering_statistics/mapping_filtering_statistics.tsv",
            category="Quality Control",
            subcategory="Mapping Statistics",
            labels={
                "type": "Statistics Table",
                "stage": "Mapping",
                "format": "TSV"
            }
        )
    log:
        f"logs/{project_name}/quality_control/mapping_filtering_statistics.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
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
        expand(rules.merge_strand_insertions.output, sample=samples, timepoint=timepoints, condition=conditions)
    output:
        report(
            f"reports/{project_name}/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
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
        f"logs/{project_name}/quality_control/PBL_PBR_correlation_analysis.log"
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
        branch(
            config["merge_similar_timepoints"],
            expand(f"results/{project_name}/11_merged/{{sample}}_{{condition}}.merged.tsv", sample=samples, condition=conditions),
            expand(rules.concat_timepoints.output.Reads, sample=samples, condition=conditions)
        )
    output:
        report(
            f"reports/{project_name}/read_count_distribution_analysis/read_count_distribution_analysis.pdf",
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
        f"logs/{project_name}/quality_control/read_count_distribution_analysis.log"
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
            f"reports/{project_name}/insertion_orientation_analysis/insertion_orientation_analysis.pdf",
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
        f"logs/{project_name}/quality_control/insertion_orientation_analysis.log"
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
        insertion_data = f"results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
        annotation = rules.concat_counts_and_annotations.output.annotations
    output:
        table = report(
            f"reports/{project_name}/insertion_density_analysis/insertion_density_analysis.tsv",
            category="Quality Control",
            subcategory="Density Analysis",
            labels={
                "type": "Statistics Table",
                "stage": "Insertion Analysis",
                "format": "TSV",
                "metric": "Insertion Density"
            }
        ),
        plot = report(
            f"reports/{project_name}/insertion_density_analysis/insertion_density_analysis_histograms.pdf",
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
        f"logs/{project_name}/quality_control/insertion_density_analysis.log"
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

# Insertion-level depletion LFC and curve features analysis
# -----------------------------------------------------
rule insertion_level_depletion_LFC_and_curve_features_analysis:
    input:
        rules.insertion_level_curve_fitting.output
    output:
        report(
            f"reports/{project_name}/depletion_LFC_and_curve_features_analysis/insertion_level_depletion_LFC_and_curve_features_analysis.pdf"
        )
    log:
        f"logs/{project_name}/quality_control/insertion_level_depletion_LFC_and_curve_features_analysis.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing insertion level depletion LFC and curve features analysis..."
    shell:
        """
        python workflow/scripts/quality_control/distribution_of_curve_fitting_results.py -i {input} -o {output} &> {log}
        """

# Gene-level depletion analysis and curve features analysis
# -----------------------------------------------------
use rule insertion_level_depletion_LFC_and_curve_features_analysis as gene_level_depletion_and_curve_features_analysis with:
    input:
        rules.gene_level_curve_fitting.output
    output:
        report(
            f"reports/{project_name}/depletion_LFC_and_curve_features_analysis/gene_level_depletion_and_curve_features_analysis.pdf"
        )
    log:
        f"logs/{project_name}/quality_control/gene_level_depletion_and_curve_features_analysis.log"

# Gene coverage analysis
# -----------------------------------------------------
rule gene_coverage_analysis:
    input:
        covered_genes = rules.gene_level_curve_fitting.output,
        all_genes = rules.download_pombase_data.output.gene_IDs_names_products.format(release_version=config["Pombase_release_version"])
    output:
        directory(f"reports/{project_name}/gene_coverage_analysis")
    log:
        f"logs/{project_name}/quality_control/gene_coverage_analysis.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing gene coverage analysis..."
    shell:
        """
        python workflow/scripts/quality_control/gene_coverage_analysis.py -c {input.covered_genes} -a {input.all_genes} -o {output} &> {log}
        """