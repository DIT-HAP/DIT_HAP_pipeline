# MultiQC for preprocessing reports
# -----------------------------------------------------
rule multiqc_preprocessing:
    input:
        fastp_json=expand(rules.fastp_preprocessing.output.json, sample=samples, timepoint=timepoints, condition=conditions),
        junction_classification_json=expand(rules.junction_classification.output.json, sample=samples, timepoint=timepoints, condition=conditions),
        fastqc=expand(f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.{{read}}_fastqc.zip", sample=samples, timepoint=timepoints, condition=conditions, read=["PBL_1", "PBL_2", "PBR_1", "PBR_2"]),
        bam_stats=expand(f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.{{type}}.txt", sample=samples, timepoint=timepoints, condition=conditions, fragment=["PBL", "PBR"], type=["stats", "flagstat", "idxstats"]),
        picard_metrics=expand(f"reports/{project_name}/picard_insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.txt", sample=samples, timepoint=timepoints, condition=conditions, fragment=["PBL", "PBR"])
    output:
        report(
            f"reports/{project_name}/multiqc/quality_control_multiqc_report.html",
            category="Quality Control",
            labels={
                "name": "1. Preprocessing MultiQC Report",
                "type": "MultiQC Report",
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
            {input.junction_classification_json} \
            {input.fastqc} \
            {input.bam_stats} \
            {input.picard_metrics} \
            &> {log}
        """

# Mapping filtering statistics
# -----------------------------------------------------
# Concat filtering statistics
# ---------------------------
rule mapping_filtering_statistics:
    input:
        expand(rules.filter_aligned_reads.log, sample=samples, timepoint=timepoints, condition=conditions)
    output:
        f"reports/{project_name}/mapping_filtering_statistics/mapping_filtering_statistics.tsv"
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

# Datavzrd report for mapping filtering statistics
# ---------------------------
rule datavzrd_mapping_filtering_statistics:
    input:
        config="workflow/reports/datavzrd/mapping_filtering_statistics.yaml",
        table=rules.mapping_filtering_statistics.output[0]
    params:
        extra="",
    output:
        report(
            directory(f"reports/{project_name}/mapping_filtering_statistics/datavzrd_mapping_filtering_statistics"),
            htmlindex="index.html",
            category="Quality Control",
            labels={
                "name": "2. Mapping Filtering Statistics",
                "type": "Datavzrd Report",
                "format": "Datavzrd HTML"
            }
        ),
    log:
        f"logs/{project_name}/quality_control/mapping_filtering_statistics_datavzrd.log"
    wrapper:
        f"{snakemake_wrapper_version}/utils/datavzrd"
        
# PBL-PBR correlation analysis
# -----------------------------------------------------
rule PBL_PBR_correlation_analysis:
    input:
        expand(rules.merge_strand_insertions.output, sample=samples, timepoint=timepoints, condition=conditions)
    output:
        report(
            f"reports/{project_name}/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
            category="Quality Control",
            labels={
                "name": "3. PBL-PBR Correlation Analysis",
                "type": "Correlation Plot",
                "format": "PDF",
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
            labels={
                "name": "4. Read Count Distribution Analysis",
                "type": "Distribution Plot",
                "format": "PDF",
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
            labels={
                "name": "5. Insertion Orientation Analysis",
                "type": "Correlation Plot",
                "format": "PDF",
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
        table = f"reports/{project_name}/insertion_density_analysis/insertion_density_analysis.tsv",
        plot = report(
            f"reports/{project_name}/insertion_density_analysis/insertion_density_analysis_histograms.pdf",
            category="Quality Control",
            labels={
                "name": "6. Insertion Density",
                "type": "Distribution Plot",
                "format": "PDF"
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
# Insertion density datavzrd report
# ---------------------------
rule datavzrd_insertion_density_analysis:
    input:
        config="workflow/reports/datavzrd/insertion_density_analysis.yaml",
        table=rules.insertion_density_analysis.output.table
    params:
        extra="",
    output:
        report(
            directory(f"reports/{project_name}/insertion_density_analysis/datavzrd_insertion_density_analysis"),
            htmlindex="index.html",
            category="Quality Control",
            labels={
                "name": "6. Insertion Density",
                "type": "Datavzrd Report",
                "format": "Datavzrd HTML"
            }
        ),
    log:
        f"logs/{project_name}/quality_control/insertion_density_analysis_datavzrd.log"
    wrapper:
        f"{snakemake_wrapper_version}/utils/datavzrd"

# Insertion-level depletion LFC and curve features analysis
# -----------------------------------------------------
rule insertion_level_depletion_LFC_and_curve_features_analysis:
    input:
        rules.insertion_level_curve_fitting.output
    output:
        report(
            f"reports/{project_name}/depletion_LFC_and_curve_features_analysis/insertion_level_depletion_LFC_and_curve_features_analysis.pdf",
            category="Insertion-level results",
            labels={
                "name": "Insertion-level Depletion LFC and Curve Features Analysis",
                "type": "Distribution Plot",
                "format": "PDF",
            }
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
            f"reports/{project_name}/depletion_LFC_and_curve_features_analysis/gene_level_depletion_and_curve_features_analysis.pdf",
            category="Gene-level results",
            labels={
                "name": "Gene-level Depletion LFC and Curve Features Analysis",
                "type": "Distribution Plot",
                "format": "PDF",
            }
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
        report(
            directory(f"reports/{project_name}/gene_coverage_analysis"),
            patterns=["{name}.png"],
            category="Gene-level results",
            labels={
                "name": "Gene Coverage Analysis",
                "type": "Donut Plot",
                "format": "PNG",
            }
        )
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