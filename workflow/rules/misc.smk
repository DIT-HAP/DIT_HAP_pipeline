# Annotate the insertions with the non-coding genes
# -----------------------------------------------------
rule insertion_annotation_with_non_coding_genes:
    input:
        insertions = rules.insertion_level_curve_fitting.output,
        non_coding_genes = rules.extract_genome_region.output.non_coding_rna_bed.format(release_version=config["Pombase_release_version"])
    output:
       f"results/{project_name}/19_insertion_in_non_coding_genes/annotated_insertions_in_non_coding_genes.tsv"
    log:
        f"logs/{project_name}/misc/insertion_annotation_with_non_coding_genes.log"
    conda:
        "../envs/pybedtools.yml"
    message:
        "*** Annotating insertions with non-coding genes..."
    shell:
        """
        python workflow/scripts/preprocessing/annotate_genomic_features.py -i {input.insertions} -g {input.non_coding_genes} -o {output} &> {log}
        """

# Gene-level depletion analysis for insertions in non-coding genes
# -----------------------------------------------------
use rule gene_level_depletion_analysis as gene_level_depletion_analysis_non_coding_genes with:
    input:
        lfc_path = f"results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
        weights_path = branch(
            config.get("use_DEseq2_for_biological_replicates", False),
            f"results/{project_name}/14_insertion_level_depletion_analysis/padj.tsv",
            f"results/{project_name}/15_insertion_level_curve_fitting/insertions_LFC_fitted_with_r_square_as_weights.tsv"
        ),
        annotations_path = rules.insertion_annotation_with_non_coding_genes.output
    output:
        all_statistics = f"results/{project_name}/19_insertion_in_non_coding_genes/Gene_level_statistics.tsv",
        LFC = f"results/{project_name}/19_insertion_in_non_coding_genes/LFC.tsv"
    log:
        f"logs/{project_name}/misc/gene_level_depletion_analysis_non_coding_genes.log"
    message:
        "*** Running gene-level depletion analysis for insertions in non-coding genes..."

# Gene-level curve fitting for insertions in non-coding genes
# -----------------------------------------------------
use rule gene_level_curve_fitting as gene_level_curve_fitting_non_coding_genes with:
    input:
        LFC = rules.gene_level_depletion_analysis_non_coding_genes.output.LFC
    output:
        f"results/{project_name}/19_insertion_in_non_coding_genes/Non_coding_genes_Gene_level_statistics_fitted.tsv"
    log:
        f"logs/{project_name}/misc/gene_level_curve_fitting_non_coding_genes.log"
    message:
        "*** Running gene-level curve fitting for insertions in non-coding genes..."