# Annotate the insertions with the non-coding genes
# -----------------------------------------------------
rule insertion_annotation_with_non_coding_genes:
    input:
        insertions = rules.insertion_level_depletion_analysis.output.LFC,
        non_coding_genes = rules.extract_genome_region.output.non_coding_rna_bed.format(release_version=config["Pombase_release_version"])
    output:
       f"results/{project_name}/19_insertion_in_non_coding_genes/annotated_insertions_in_non_coding_genes.tsv"
    log:
        "logs/misc/insertion_annotation_with_non_coding_genes.log"
    conda:
        "../envs/pybedtools.yml"
    message:
        "*** Annotating insertions with non-coding genes..."
    shell:
        """
        python workflow/scripts/preprocessing/annotate_insertions.py -i {input.insertions} -g {input.non_coding_genes} -o {output} &> {log}
        """

# Gene-level depletion analysis for insertions in non-coding genes
# -----------------------------------------------------
use rule gene_level_depletion_analysis as gene_level_depletion_analysis_non_coding_genes with:
    input:
        lfc_path = rules.insertion_level_depletion_analysis.output.all_statistics,
        annotations_path = rules.insertion_annotation_with_non_coding_genes.output
    output:
        all_statistics = f"results/{project_name}/19_insertion_in_non_coding_genes/Gene_level_statistics.csv",
        LFC = f"results/{project_name}/19_insertion_in_non_coding_genes/LFC.csv"
    log:
        "logs/misc/gene_level_depletion_analysis_non_coding_genes.log"
    message:
        "*** Running gene-level depletion analysis for insertions in non-coding genes..."

# Gene-level curve fitting for insertions in non-coding genes
# -----------------------------------------------------
use rule gene_level_curve_fitting as gene_level_curve_fitting_non_coding_genes with:
    input:
        LFC = rules.gene_level_depletion_analysis_non_coding_genes.output.LFC
    output:
        f"results/{project_name}/19_insertion_in_non_coding_genes/Non_coding_genes_Gene_level_statistics_fitted.csv"
    log:
        "logs/misc/gene_level_curve_fitting_non_coding_genes.log"