rule extract_coding_region_from_gff3:
    input:
        gff3=gff3_file,
        peptide_stats=peptide_stats
    output:
        output_dir/"Coding_genes.bed"
    shell:
        """
        python workflow/scripts/making_genome_intervals/extract_coding_region_from_gff3.py \
            -g {input.gff3} \
            -p {input.peptide_stats} \
            -o {output}
        """

rule find_overlapped_gene_regions:
    input:
        rules.extract_coding_region_from_gff3.output
    output:
        output_dir/"Overlapped_gene_regions.bed"
    conda:
        "pybedtools"
    shell:
        """
        python workflow/scripts/making_genome_intervals/find_overlapped_gene_regions.py \
            -i {input} \
            -o {output}
        """

rule annotate_intergenic_regions:
    input:
        coding_genes = rules.extract_coding_region_from_gff3.output,
        overlapped_gene_regions = rules.find_overlapped_gene_regions.output,
        fai = fai_file
    output:
        output_dir / "Intergenic_regions.bed"
    conda:
        "pybedtools"
    shell:
        """
        python workflow/scripts/making_genome_intervals/annotate_intergenic_regions.py \
            -ic {input.coding_genes} \
            -iol {input.overlapped_gene_regions} \
            -f {input.fai} \
            -o {output}
        """

rule merge_genome_regions:
    input:
        coding_genes = rules.extract_coding_region_from_gff3.output,
        overlapped_gene_regions = rules.find_overlapped_gene_regions.output,
        intergenic_regions = rules.annotate_intergenic_regions.output
    output:
        not_separated_output = output_dir / "Genome_regions_overlapped_regions_not_separated.bed",
        separated_output = output_dir / "Genome_regions_overlapped_regions_separated.bed"
    conda:
        "pybedtools"
    shell:
        """
        python workflow/scripts/making_genome_intervals/merge_genome_regions.py \
            -ic {input.coding_genes} \
            -iol {input.overlapped_gene_regions} \
            -iIGR {input.intergenic_regions} \
            -o1 {output.not_separated_output} \
            -o2 {output.separated_output}
        """

rule extract_CDS_intron_from_coding_region_gff3:
    input:
        gff3=gff3_file,
        peptide_stats=peptide_stats
    output:
        output_dir / "Coding_gene_CDS_and_intron.bed"
    shell:
        """
        python workflow/scripts/making_genome_intervals/extract_CDS_intron_from_coding_region_gff3.py \
            -g {input.gff3} \
            -p {input.peptide_stats} \
            -o {output}
        """

rule merge_CDS_intron_intergenic:
    input:
        CDS_intro = rules.extract_CDS_intron_from_coding_region_gff3.output,
        intergenic_regions = rules.annotate_intergenic_regions.output
    output:
        output_dir / "Genome_regions_CDS_intron_IGR.bed"
    shell:
        """
        python workflow/scripts/making_genome_intervals/merge_CDS_intron_intergenic.py \
            -ci {input.CDS_intro} \
            -i {input.intergenic_regions} \
            -o {output}
        """

rule annotate_the_genome_regions:
    input:
        genome_regions = rules.merge_CDS_intron_intergenic.output,
        gene_product = gene_product_file,
        gene_essentiality = gene_essentiality_file,
    output:
        output_dir / "Genome_regions_CDS_intron_IGR_annotated.bed"
    shell:
        """
        python workflow/scripts/making_genome_intervals/annotate_the_genome_regions.py \
            -i {input.genome_regions} \
            -n {input.gene_product} \
            -e {input.gene_essentiality} \
            -o {output}
        """

