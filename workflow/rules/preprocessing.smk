# Fastp preprocessing for QC and adapter trimming
# -----------------------------------------------------
rule fastp_preprocessing:
    input:
        fq1=lambda wildcards: sample_sheet_dict[wildcards.sample][wildcards.timepoint][wildcards.condition]["fq1"],
        fq2=lambda wildcards: sample_sheet_dict[wildcards.sample][wildcards.timepoint][wildcards.condition]["fq2"]
    output:
        fq1=temp(f"results/{project_name}/1_fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp_1.fq.gz"),
        fq2=temp(f"results/{project_name}/1_fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp_2.fq.gz"),
        html=f"reports/{project_name}/fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp.html",
        json=f"reports/{project_name}/fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp.json"
    log:
        f"logs/{project_name}/preprocessing/fastp/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/fastp.yml"
    params:
        adapter_sequence=config["adapter_sequence"],
        adapter_sequence_r2=config["adapter_sequence_r2"]
    threads: 6
    message:
        "*** Preprocessing fastp for {input.fq1} and {input.fq2}..."
    shell:
        """
        fastp --adapter_sequence {params.adapter_sequence} \
              --adapter_sequence_r2 {params.adapter_sequence_r2} \
              --disable_quality_filtering \
              --disable_length_filtering \
              --overrepresentation_analysis \
              --correction \
              -j {output.json} \
              -h {output.html} \
              --thread {threads} \
              --in1 {input.fq1} \
              --in2 {input.fq2} \
              --out1 {output.fq1} \
              --out2 {output.fq2} &> {log}
        """

# Cutadapt preprocessing for classification of PBL and PBR
# -----------------------------------------------------
rule junction_classification:
    input:
        fq1 = rules.fastp_preprocessing.output.fq1,
        fq2 = rules.fastp_preprocessing.output.fq2
    output:
        PBL_r1=f"results/{project_name}/2_junction_classification/{{sample}}_{{timepoint}}_{{condition}}.PBL_1.fq.gz",
        PBL_r2=f"results/{project_name}/2_junction_classification/{{sample}}_{{timepoint}}_{{condition}}.PBL_2.fq.gz",
        PBR_r1=f"results/{project_name}/2_junction_classification/{{sample}}_{{timepoint}}_{{condition}}.PBR_1.fq.gz",
        PBR_r2=f"results/{project_name}/2_junction_classification/{{sample}}_{{timepoint}}_{{condition}}.PBR_2.fq.gz",
        json=f"reports/{project_name}/junction_classification/{{sample}}_{{timepoint}}_{{condition}}.json",
    log:
        f"logs/{project_name}/preprocessing/junction_classification/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/cutadapt.yml"
    params:
        PBL_adapter=config["PBL_adapter"],
        PBR_adapter=config["PBR_adapter"],
        PBL_reverseComplement_adapter=config["PBL_reverseComplement_adapter"],
        PBR_reverseComplement_adapter=config["PBR_reverseComplement_adapter"],
        output_folder = f"results/{project_name}/2_junction_classification"
    threads: 6
    message:
        "*** Junction classification {input.fq1} and {input.fq2}..."
    shell:
        """
        cutadapt --cores {threads} \
                 -q 15 \
                 --overlap 15 \
                 -g PBL={params.PBL_adapter} \
                 -g PBR={params.PBR_adapter} \
                 -A PBL={params.PBL_reverseComplement_adapter} \
                 -A PBR={params.PBR_reverseComplement_adapter} \
                 -o {params.output_folder}/{wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}.{{name}}_1.fq.gz \
                 -p {params.output_folder}/{wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}.{{name}}_2.fq.gz \
                 --json {output.json} \
                 {input.fq1} {input.fq2} &> {log}
        """

# fastqc for junction classified reads
# -----------------------------------------------------
rule fastqc_junction_classification:
    input:
        PBL_r1=rules.junction_classification.output.PBL_r1,
        PBL_r2=rules.junction_classification.output.PBL_r2,
        PBR_r1=rules.junction_classification.output.PBR_r1,
        PBR_r2=rules.junction_classification.output.PBR_r2
    output:
        # Explicitly specify output files for better dependency tracking
        PBL_r1_html=f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBL_1_fastqc.html",
        PBL_r1_zip=f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBL_1_fastqc.zip",
        PBL_r2_html=f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBL_2_fastqc.html",
        PBL_r2_zip=f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBL_2_fastqc.zip",
        PBR_r1_html=f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBR_1_fastqc.html",
        PBR_r1_zip=f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBR_1_fastqc.zip",
        PBR_r2_html=f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBR_2_fastqc.html",
        PBR_r2_zip=f"reports/{project_name}/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBR_2_fastqc.zip"
    log:
        f"logs/{project_name}/preprocessing/fastqc/{{sample}}_{{timepoint}}_{{condition}}_fastqc_demultiplexed.log"
    conda:
        "../envs/fastqc.yml"
    params:
        output_dir = f"reports/{project_name}/fastqc"
    threads: 4  # Increased threads for parallel processing
    message:
        "*** Running FastQC for junction classified paired-end reads for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        echo "*** Running FastQC for PBL R1 reads..." > {log}
        fastqc --threads {threads} --noextract {input.PBL_r1} -o {params.output_dir} &>> {log}
        
        echo "*** Running FastQC for PBL R2 reads..." >> {log}
        fastqc --threads {threads} --noextract {input.PBL_r2} -o {params.output_dir} &>> {log}
        
        echo "*** Running FastQC for PBR R1 reads..." >> {log}
        fastqc --threads {threads} --noextract {input.PBR_r1} -o {params.output_dir} &>> {log}
        
        echo "*** Running FastQC for PBR R2 reads..." >> {log}
        fastqc --threads {threads} --noextract {input.PBR_r2} -o {params.output_dir} &>> {log}
        
        echo "*** FastQC analysis completed for all reads" >> {log}
        """

# BWA-MEM mapping, sorting and indexing
# -----------------------------------------------------
rule bwa_mem_mapping:
    input:
        ref=rules.download_pombase_data.output.fasta.format(release_version=config["Pombase_release_version"]),
        ref_index=expand(rules.bwa_index.output, release_version=config["Pombase_release_version"]), # ensure the reference index is built
        PBL_fq1=rules.junction_classification.output.PBL_r1,
        PBL_fq2=rules.junction_classification.output.PBL_r2,
        PBR_fq1=rules.junction_classification.output.PBR_r1,
        PBR_fq2=rules.junction_classification.output.PBR_r2,
        PBL_fastqc_r1=rules.fastqc_junction_classification.output.PBL_r1_html, # make sure the quality control before mapping
        PBL_fastqc_r2=rules.fastqc_junction_classification.output.PBL_r2_html, # make sure the quality control before mapping
        PBR_fastqc_r1=rules.fastqc_junction_classification.output.PBR_r1_html, # make sure the quality control before mapping
        PBR_fastqc_r2=rules.fastqc_junction_classification.output.PBR_r2_html # make sure the quality control before mapping
    output:
        PBL=f"results/{project_name}/3_mapped/{{sample}}_{{timepoint}}_{{condition}}.PBL.name_sorted.bam",
        PBR=f"results/{project_name}/3_mapped/{{sample}}_{{timepoint}}_{{condition}}.PBR.name_sorted.bam"
    log:
        f"logs/{project_name}/preprocessing/bwa_mem_mapping/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/bwa_mapping.yml"
    threads: 8
    message:
        "*** Mapping {input.PBL_fq1} and {input.PBL_fq2} to {input.ref}..."
    shell:
        """
        echo "*** Mapping PBL reads..." > {log}
        bwa-mem2 mem -t {threads} {input.ref} {input.PBL_fq1} {input.PBL_fq2} 2>> {log} | samtools sort -n -@ {threads} -O BAM -o {output.PBL} &>> {log}
        echo "*** Mapping PBR reads..." >> {log}
        bwa-mem2 mem -t {threads} {input.ref} {input.PBR_fq1} {input.PBR_fq2} 2>> {log} | samtools sort -n -@ {threads} -O BAM -o {output.PBR} &>> {log}
        """

# Samtools sorting and indexing for BAM files
# -----------------------------------------------------
rule samtools_sorting_and_indexing:
    input:
        PBL=rules.bwa_mem_mapping.output.PBL,
        PBR=rules.bwa_mem_mapping.output.PBR,
        ref_index=rules.samtools_faidx.output[0].format(release_version=config["Pombase_release_version"])
    output:
        PBL_sorted=f"results/{project_name}/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.PBL.sorted.bam",
        PBR_sorted=f"results/{project_name}/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.PBR.sorted.bam",
        PBL_index=f"results/{project_name}/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.PBL.sorted.bam.bai",
        PBR_index=f"results/{project_name}/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.PBR.sorted.bam.bai"
    log:
        f"logs/{project_name}/preprocessing/samtools_sorting_and_indexing/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/bwa_mapping.yml"
    threads: 2
    message:
        "*** Sorting and indexing {input.PBL} and {input.PBR}..."
    shell:
        """
        echo "*** Sorting PBL reads..." > {log}
        samtools sort -@ {threads} {input.PBL} -O BAM -o {output.PBL_sorted} &>> {log}
        echo "*** Indexing PBL reads..." >> {log}
        samtools index -@ {threads} {output.PBL_sorted} &>> {log}
        echo "*** Sorting PBR reads..." >> {log}
        samtools sort -@ {threads} {input.PBR} -O BAM -o {output.PBR_sorted} &>> {log}
        echo "*** Indexing PBR reads..." >> {log}
        samtools index -@ {threads} {output.PBR_sorted} &>> {log}
        """


# samtools for mapping statistics
# -----------------------------------------------------
rule samtools_mapping_statistics:
    input:
        PBL=rules.samtools_sorting_and_indexing.output.PBL_sorted,
        PBR=rules.samtools_sorting_and_indexing.output.PBR_sorted
    output:
        PBL_stats=f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBL.stats.txt",
        PBR_stats=f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBR.stats.txt",
        PBL_flagstat=f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBL.flagstat.txt",
        PBR_flagstat=f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBR.flagstat.txt",
        PBL_idxstats=f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBL.idxstats.txt",
        PBR_idxstats=f"reports/{project_name}/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBR.idxstats.txt"
    log:
        f"logs/{project_name}/preprocessing/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/bwa_mapping.yml"
    threads: 2
    message:
        "*** Running Samtools for mapping statistics for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        echo "*** Running Samtools for mapping statistics for PBL..." > {log}
        samtools stats {input.PBL} > {output.PBL_stats} 2>> {log}
        samtools flagstat {input.PBL} > {output.PBL_flagstat} 2>> {log}
        samtools idxstats {input.PBL} > {output.PBL_idxstats} 2>> {log}
        echo "*** Running Samtools for mapping statistics for PBR..." >> {log}
        samtools stats {input.PBR} > {output.PBR_stats} 2>> {log}
        samtools flagstat {input.PBR} > {output.PBR_flagstat} 2>> {log}
        samtools idxstats {input.PBR} > {output.PBR_idxstats} 2>> {log}
        """

# Picard collect insertion size metrics
# -----------------------------------------------------
rule insert_size:
    input:
        f"results/{project_name}/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.sorted.bam"
    output:
        txt=f"reports/{project_name}/picard_insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.txt",
        pdf=f"reports/{project_name}/picard_insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.pdf",
    log:
        f"logs/{project_name}/preprocessing/insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.log",
    params:
        # optional parameters (e.g. relax checks as below)
        extra="--VALIDATION_STRINGENCY LENIENT --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        f"{snakemake_wrapper_version}/bio/picard/collectinsertsizemetrics"

# transform bam to tsv for filtering and extracting insertions
# -----------------------------------------------------
rule bam_to_tsv:
    input:
        PBL=rules.bwa_mem_mapping.output.PBL,
        PBR=rules.bwa_mem_mapping.output.PBR
    output:
        PBL_tsv=f"results/{project_name}/5_tabulated/{{sample}}_{{timepoint}}_{{condition}}.PBL.tsv",
        PBR_tsv=f"results/{project_name}/5_tabulated/{{sample}}_{{timepoint}}_{{condition}}.PBR.tsv",
    log:
        f"logs/{project_name}/preprocessing/bam_to_tsv/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/pysam.yml"
    threads: 8
    resources:
        mem_mb=200000
    message:
        "*** Transforming BAM to TSV for filtering and extracting insertions for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        echo "*** Transforming PBL BAM to TSV..." > {log}
        python workflow/scripts/preprocessing/parse_bam_to_tsv.py -i {input.PBL} -o {output.PBL_tsv} -t {threads} &>> {log}
        echo "*** Transforming PBR BAM to TSV..." >> {log}
        python workflow/scripts/preprocessing/parse_bam_to_tsv.py -i {input.PBR} -o {output.PBR_tsv} -t {threads} &>> {log}
        """

# Filtering aligned read pairs for PBL and PBR
# -----------------------------------------------------
rule filter_aligned_reads:
    input:
        PBL_tsv=rules.bam_to_tsv.output.PBL_tsv,
        PBR_tsv=rules.bam_to_tsv.output.PBR_tsv
    output:
        PBL_filtered=f"results/{project_name}/6_filtered/{{sample}}_{{timepoint}}_{{condition}}.PBL.filtered.tsv",
        PBR_filtered=f"results/{project_name}/6_filtered/{{sample}}_{{timepoint}}_{{condition}}.PBR.filtered.tsv"
    log:
        f"logs/{project_name}/preprocessing/filter_aligned_reads/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Filtering aligned read pairs for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    params:
        snakemake_config_file=snakemake_config_file,
        chunk_size=config["chunk_size"]
    shell:
        """
        python workflow/scripts/preprocessing/filter_aligned_reads.py -i {input.PBL_tsv} \
                                                                      -o {output.PBL_filtered} \
                                                                      -c {params.chunk_size} \
                                                                      --config {params.snakemake_config_file} &> {log}
        python workflow/scripts/preprocessing/filter_aligned_reads.py -i {input.PBR_tsv} \
                                                                      -o {output.PBR_filtered} \
                                                                      -c {params.chunk_size} \
                                                                      --config {params.snakemake_config_file} &>> {log}
        """

# Extract insertion sites from filtered read pairs
# -----------------------------------------------------
rule extract_insertion_sites:
    input:
        PBL_filtered=rules.filter_aligned_reads.output.PBL_filtered,
        PBR_filtered=rules.filter_aligned_reads.output.PBR_filtered
    output:
        PBL_insertions=f"results/{project_name}/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.PBL.tsv",
        PBR_insertions=f"results/{project_name}/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.PBR.tsv"
    log:
        f"logs/{project_name}/preprocessing/extract_insertion_sites/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    params:
        chunk_size=config["chunk_size"]
    message:
        "*** Extracting insertion sites for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        echo "*** Extracting PBL insertion sites..." > {log}
        python workflow/scripts/preprocessing/extract_insertion_sites.py -i {input.PBL_filtered} -o {output.PBL_insertions} -c {params.chunk_size} &>> {log}
        echo "*** Extracting PBR insertion sites..." >> {log}
        python workflow/scripts/preprocessing/extract_insertion_sites.py -i {input.PBR_filtered} -o {output.PBR_insertions} -c {params.chunk_size} &>> {log}
        """

# Merge strand-specific insertions from PBL and PBR reads
# -----------------------------------------------------
rule merge_strand_insertions:
    input:
        PBL_insertions=rules.extract_insertion_sites.output.PBL_insertions,
        PBR_insertions=rules.extract_insertion_sites.output.PBR_insertions
    output:
        f"results/{project_name}/8_merged/{{sample}}_{{timepoint}}_{{condition}}.tsv"
    log:
        f"logs/{project_name}/preprocessing/merge_strand_insertions/{{sample}}_{{timepoint}}_{{condition}}.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Merging strand-specific insertions for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/preprocessing/merge_strand_insertions.py \
            -i {input.PBL_insertions} \
            -j {input.PBR_insertions} \
            -o {output} &> {log}
        """

# concat the insertions from the different timepoints
# -----------------------------------------------------
rule concat_timepoints:
    input:
        counts = lambda wildcards: expand(rules.merge_strand_insertions.output, sample=wildcards.sample, timepoint=timepoints, condition=wildcards.condition),
        ref = rules.download_pombase_data.output.fasta.format(release_version=config["Pombase_release_version"])
    output:
        PBL = f"results/{project_name}/9_concatenated/{{sample}}_{{condition}}.PBL.tsv",
        PBR = f"results/{project_name}/9_concatenated/{{sample}}_{{condition}}.PBR.tsv",
        Reads = f"results/{project_name}/9_concatenated/{{sample}}_{{condition}}.Reads.tsv"
    log:
        f"logs/{project_name}/preprocessing/concat_timepoints/{{sample}}_{{condition}}.log"
    conda:
        "../envs/biopython.yml"
    params:
        timepoints = " ".join(timepoints)
    shell:
        """
        python workflow/scripts/preprocessing/concatenate_timepoint_data.py \
                    -s {wildcards.sample}_{wildcards.condition} \
                    -i {input.counts} \
                    -tp {params.timepoints} \
                    -g {input.ref} \
                    -ol {output.PBL} \
                    -or {output.PBR} \
                    -o {output.Reads} &> {log}
        """        

# annotate the insertions with the genome region
# -----------------------------------------------------
rule annotate_insertions:
    input:
        insertions=rules.concat_timepoints.output.Reads,
        genome_region=rules.extract_genome_region.output.genome_intervals_bed.format(release_version=config["Pombase_release_version"])
    output:
        f"results/{project_name}/10_annotated/{{sample}}_{{condition}}.annotated.tsv"
    log:
        f"logs/{project_name}/preprocessing/annotate_insertions/{{sample}}_{{condition}}.log"
    conda:
        "../envs/pybedtools.yml"
    message:
        "*** Annotating insertions for {wildcards.sample}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/preprocessing/annotate_genomic_features.py -i {input.insertions} -g {input.genome_region} -o {output} &> {log}
        """

# Optional: merge the similar time points
# -----------------------------------------------------
if config["merge_similar_timepoints"]:
    rule merge_similar_timepoints:
        input:
            rules.concat_timepoints.output.Reads
        output:
            f"results/{project_name}/11_merged/{{sample}}_{{condition}}.merged.tsv"
        log:
            f"logs/{project_name}/preprocessing/merge_similar_timepoints/{{sample}}_{{condition}}.log"
        params:
            similar_timepoints = config["similar_timepoints"],
            merged_timepoint = config["merged_timepoint"],
            drop_columns = config["drop_columns"]
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** Merging similar time points for {wildcards.sample}_{wildcards.condition}..."
        run:
            import pandas as pd
            df = pd.read_csv(input[0], sep="\t", index_col=[0,1,2,3])
            print("*** Before merging similar time points")
            print(df.head())
            print(df.columns)
            df[params.merged_timepoint] = df[params.similar_timepoints].sum(axis=1)
            print("*** After merging similar time points")
            print("Merged timepoint: ", params.merged_timepoint)
            print("Drop columns: ", params.drop_columns)
            print(df.head())
            print(df.columns)
            df.drop(columns=params.drop_columns, inplace=True)
            df.sort_index(axis=1, inplace=True)
            print("*** After sorting the index")
            print(df.head())
            print(df.columns)
            df.to_csv(output[0], sep="\t", index=True)

# concat the sample counts and annotations
# -----------------------------------------------------
rule concat_counts_and_annotations:
    input:
        counts = branch(
            config["merge_similar_timepoints"],
            expand(f"results/{project_name}/11_merged/{{sample}}_{{condition}}.merged.tsv", sample=samples, condition=conditions),
            expand(rules.concat_timepoints.output.Reads, sample=samples, condition=conditions)
        ),
        annotations = expand(rules.annotate_insertions.output, sample=samples, condition=conditions)
    output:
        counts = f"results/{project_name}/12_concatenated/raw_reads.tsv",
        annotations = f"results/{project_name}/12_concatenated/annotations.tsv"
    log:
        f"logs/{project_name}/preprocessing/concat_counts_and_annotations.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Concatenating counts and annotations..."
    run:
        import pandas as pd

        all_counts = str(input.counts).split(" ")
        all_annotations = str(input.annotations).split(" ")
        print("*** All count files: ", all_counts)
        print("*** All annotation files: ", all_annotations)

        counts_df = {}
        for input_counts in all_counts:
            counts_df[Path(input_counts).name.split(".")[0]] = pd.read_csv(input_counts, sep="\t", index_col=[0,1,2,3])
            print("*** Counts file: ", input_counts)
            print("*** Counts file shape: ", counts_df[Path(input_counts).name.split(".")[0]].shape[0])
        counts_df = pd.concat(counts_df, axis=1).rename_axis(["Sample", "Timepoint"], axis=1)
        print("*** Counts dataframe shape: ", counts_df.shape[0])

        annotations_df = {}
        for input_annotations in all_annotations:
            annotations_df[Path(input_annotations).name.split(".")[0]] = pd.read_csv(input_annotations, index_col=[0,1,2,3], sep="\t")
            print("*** Annotations file: ", input_annotations)
            print("*** Annotations file shape: ", annotations_df[Path(input_annotations).name.split(".")[0]].shape[0])
        annotations_df = pd.concat(list(annotations_df.values()), axis=0).reset_index().drop_duplicates().set_index(["Chr", "Coordinate", "Strand", "Target"])
        print("*** Annotations dataframe shape: ", annotations_df.shape[0])

        counts_df.to_csv(output.counts, sep="\t", index=True)
        annotations_df.to_csv(output.annotations, sep="\t", index=True)

        print("*** Concatenated counts and annotations completed")
        print("*** Counts dataframe shape: ", counts_df.shape[0])
        print("*** Annotations dataframe shape: ", annotations_df.shape[0])


# hard filtering the insertions
# -----------------------------------------------------
rule hard_filtering:
    input:
        rules.concat_counts_and_annotations.output.counts
    output:
        f"results/{project_name}/13_filtered/raw_reads.filtered.tsv"
    log:
        f"logs/{project_name}/preprocessing/hard_filtering.log"
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    params:
        cutoff = config["hard_filtering_cutoff"],
        init_timepoint = config["initial_time_point"]
    message:
        "*** Hard filtering the insertions..."
    shell:
        """
        python workflow/scripts/preprocessing/reads_hard_filtering.py -i {input} -o {output} -c {params.cutoff} -itp {params.init_timepoint} &> {log}
        """



        
        