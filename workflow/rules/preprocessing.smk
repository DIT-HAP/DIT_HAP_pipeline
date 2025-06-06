# Fastp preprocessing for QC and adapter trimming
# -----------------------------------------------------
rule fastp_preprocessing:
    input:
        fq1=lambda wildcards: sample_sheet_dict[wildcards.sample][wildcards.timepoint][wildcards.condition]["fq1"],
        fq2=lambda wildcards: sample_sheet_dict[wildcards.sample][wildcards.timepoint][wildcards.condition]["fq2"]
    output:
        fq1=f"results/{project_name}/1_fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp_1.fq.gz",
        fq2=f"results/{project_name}/1_fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp_2.fq.gz",
        html=f"reports/{project_name}/fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp.html",
        json=f"reports/{project_name}/fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp.json"
    log:
        "logs/preprocessing/fastp/{sample}_{timepoint}_{condition}.log"
    conda:
        "../envs/fastp.yml"
    params:
        adapter_sequence=config["adapter_sequence"],
        adapter_sequence_r2=config["adapter_sequence_r2"]
    threads: 4
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

# Cutadapt preprocessing for demultiplexing of PBL and PBR
# -----------------------------------------------------
rule demultiplexing:
    input:
        fq1 = rules.fastp_preprocessing.output.fq1,
        fq2 = rules.fastp_preprocessing.output.fq2
    output:
        PBL_r1=f"results/{project_name}/2_demultiplexed/{{sample}}_{{timepoint}}_{{condition}}.PBL_1.fq.gz",
        PBL_r2=f"results/{project_name}/2_demultiplexed/{{sample}}_{{timepoint}}_{{condition}}.PBL_2.fq.gz",
        PBR_r1=f"results/{project_name}/2_demultiplexed/{{sample}}_{{timepoint}}_{{condition}}.PBR_1.fq.gz",
        PBR_r2=f"results/{project_name}/2_demultiplexed/{{sample}}_{{timepoint}}_{{condition}}.PBR_2.fq.gz",
        json=f"reports/{project_name}/demultiplexing/{{sample}}_{{timepoint}}_{{condition}}.json",
    log:
        "logs/preprocessing/demultiplexing/{sample}_{timepoint}_{condition}.log"
    conda:
        "../envs/cutadapt.yml"
    params:
        PBL_adapter=config["PBL_adapter"],
        PBR_adapter=config["PBR_adapter"],
        PBL_reverseComplement_adapter=config["PBL_reverseComplement_adapter"],
        PBR_reverseComplement_adapter=config["PBR_reverseComplement_adapter"],
        output_folder = f"results/{project_name}/2_demultiplexed"
    threads: 4
    message:
        "*** Demultiplexing {input.fq1} and {input.fq2}..."
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

# fastqc for demultiplexed reads
# -----------------------------------------------------
rule fastqc_demultiplexed:
    input:
        PBL_r1=rules.demultiplexing.output.PBL_r1,
        PBL_r2=rules.demultiplexing.output.PBL_r2,
        PBR_r1=rules.demultiplexing.output.PBR_r1,
        PBR_r2=rules.demultiplexing.output.PBR_r2
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
        "logs/preprocessing/fastqc/{sample}_{timepoint}_{condition}_fastqc_demultiplexed.log"
    conda:
        "../envs/fastqc.yml"
    params:
        output_dir = f"reports/{project_name}/fastqc"
    threads: 4  # Increased threads for parallel processing
    message:
        "*** Running FastQC for demultiplexed paired-end reads for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
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
        ref_index=expand(rules.bwa_index.output, release_version=config["Pombase_release_version"]),
        PBL_fq1=rules.demultiplexing.output.PBL_r1,
        PBL_fq2=rules.demultiplexing.output.PBL_r2,
        PBR_fq1=rules.demultiplexing.output.PBR_r1,
        PBR_fq2=rules.demultiplexing.output.PBR_r2,
        PBL_fastqc_r1=rules.fastqc_demultiplexed.output.PBL_r1_html, # make sure the quality control before mapping
        PBL_fastqc_r2=rules.fastqc_demultiplexed.output.PBL_r2_html, # make sure the quality control before mapping
        PBR_fastqc_r1=rules.fastqc_demultiplexed.output.PBR_r1_html, # make sure the quality control before mapping
        PBR_fastqc_r2=rules.fastqc_demultiplexed.output.PBR_r2_html # make sure the quality control before mapping
    output:
        PBL=f"results/{project_name}/3_mapped/{{sample}}_{{timepoint}}_{{condition}}.PBL.name_sorted.bam",
        PBR=f"results/{project_name}/3_mapped/{{sample}}_{{timepoint}}_{{condition}}.PBR.name_sorted.bam"
    log:
        "logs/preprocessing/bwa_mem_mapping/{sample}_{timepoint}_{condition}.log"
    conda:
        "../envs/bwa_mapping.yml"
    threads: 8
    message:
        "*** Mapping {input.PBL_fq1} and {input.PBL_fq2} to {input.ref}..."
    shell:
        """
        echo "*** Mapping PBL reads..." > {log}
        bwa mem -t {threads} {input.ref} {input.PBL_fq1} {input.PBL_fq2} 2>> {log} | samtools sort -n -@ {threads} -O BAM -o {output.PBL} &>> {log}
        echo "*** Mapping PBR reads..." >> {log}
        bwa mem -t {threads} {input.ref} {input.PBR_fq1} {input.PBR_fq2} 2>> {log} | samtools sort -n -@ {threads} -O BAM -o {output.PBR} &>> {log}
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
        "logs/preprocessing/samtools_sorting_and_indexing/{sample}_{timepoint}_{condition}.log"
    conda:
        "../envs/samtools.yml"
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
        "logs/preprocessing/samtools_mapping_statistics/{sample}_{timepoint}_{condition}.log"
    conda:
        "../envs/samtools.yml"
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
        "logs/preprocessing/bam_to_tsv/{sample}_{timepoint}_{condition}.log"
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

# Filtering read pairs for PBL and PBR
# -----------------------------------------------------
rule read_pair_filtering:
    input:
        PBL_tsv=rules.bam_to_tsv.output.PBL_tsv,
        PBR_tsv=rules.bam_to_tsv.output.PBR_tsv
    output:
        PBL_filtered=f"results/{project_name}/6_filtered/{{sample}}_{{timepoint}}_{{condition}}.PBL.filtered.tsv",
        PBR_filtered=f"results/{project_name}/6_filtered/{{sample}}_{{timepoint}}_{{condition}}.PBR.filtered.tsv"
    log:
        "logs/preprocessing/read_pair_filtering/{sample}_{timepoint}_{condition}.log"
    conda:
        "../envs/tabular_operations.yml"
    message:
        "*** Filtering read pairs for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    params:
        r1_mapq_threshold=config["r1_mapq_threshold"],
        r1_ncigar_value=config["r1_ncigar_value"],
        r1_nm_threshold=config["r1_nm_threshold"],
        r1_no_sa=config["r1_no_sa"],
        r1_no_xa=config["r1_no_xa"],
        r2_mapq_threshold=config["r2_mapq_threshold"],
        r2_disable_ncigar = config["r2_disable_ncigar"],
        r2_disable_nm = config["r2_disable_nm"],
        r2_no_sa=config["r2_no_sa"],
        r2_no_xa=config["r2_no_xa"],
        require_proper_pair=config["require_proper_pair"],
        chunk_size=config["chunk_size"]
    shell:
        """
        python workflow/scripts/preprocessing/tabular_read_pair_filtering.py -i {input.PBL_tsv} \
                                                                             -o {output.PBL_filtered} \
                                                                             -c {params.chunk_size} \
                                                                             --r1-mapq-threshold {params.r1_mapq_threshold} \
                                                                             --r1-ncigar-value {params.r1_ncigar_value} \
                                                                             --r1-nm-threshold {params.r1_nm_threshold} \
                                                                             {params.r1_no_sa} \
                                                                             {params.r1_no_xa} \
                                                                             --r2-mapq-threshold {params.r2_mapq_threshold} \
                                                                             {params.r2_disable_ncigar} \
                                                                             {params.r2_disable_nm} \
                                                                             {params.r2_no_sa} \
                                                                             {params.r2_no_xa} \
                                                                             {params.require_proper_pair} &> {log}
        python workflow/scripts/preprocessing/tabular_read_pair_filtering.py -i {input.PBR_tsv} \
                                                                             -o {output.PBR_filtered} \
                                                                             -c {params.chunk_size} \
                                                                             --r1-mapq-threshold {params.r1_mapq_threshold} \
                                                                             --r1-ncigar-value {params.r1_ncigar_value} \
                                                                             --r1-nm-threshold {params.r1_nm_threshold} \
                                                                             {params.r1_no_sa} \
                                                                             {params.r1_no_xa} \
                                                                             --r2-mapq-threshold {params.r2_mapq_threshold} \
                                                                             {params.r2_disable_ncigar} \
                                                                             {params.r2_disable_nm} \
                                                                             {params.r2_no_sa} \
                                                                             {params.r2_no_xa} \
                                                                             {params.require_proper_pair} &>> {log}
        """

# extract the insertions from the filtered read pairs
# -----------------------------------------------------
rule extract_insertions:
    input:
        PBL_filtered=rules.read_pair_filtering.output.PBL_filtered,
        PBR_filtered=rules.read_pair_filtering.output.PBR_filtered
    output:
        PBL_insertions=f"results/{project_name}/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.PBL.tsv",
        PBR_insertions=f"results/{project_name}/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.PBR.tsv"
    log:
        "logs/preprocessing/extract_insertions/{sample}_{timepoint}_{condition}.log"
    conda:
        "../envs/tabular_operations.yml"
    params:
        chunk_size=config["chunk_size"]
    message:
        "*** Extracting insertions for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        echo "*** Extracting PBL insertions..." > {log}
        python workflow/scripts/preprocessing/process_insertion_counts.py -i {input.PBL_filtered} -o {output.PBL_insertions} -c {params.chunk_size} &>> {log}
        echo "*** Extracting PBR insertions..." >> {log}
        python workflow/scripts/preprocessing/process_insertion_counts.py -i {input.PBR_filtered} -o {output.PBR_insertions} -c {params.chunk_size} &>> {log}
        """

# merge the insertions from the PBL and PBR reads
# -----------------------------------------------------
rule merge_insertions:
    input:
        PBL_insertions=rules.extract_insertions.output.PBL_insertions,
        PBR_insertions=rules.extract_insertions.output.PBR_insertions
    output:
        f"results/{project_name}/8_merged/{{sample}}_{{timepoint}}_{{condition}}.tsv"
    log:
        "logs/preprocessing/merge_insertions/{sample}_{timepoint}_{condition}.log"
    conda:
        "../envs/tabular_operations.yml"
    message:
        "*** Merging insertions for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/preprocessing/merge_insertions.py \
            -il {input.PBL_insertions} \
            -ir {input.PBR_insertions} \
            -o {output} &> {log}
        """

# concat the insertions from the different timepoints
# -----------------------------------------------------
rule concat_timepoints:
    input:
        counts = lambda wildcards: expand(rules.merge_insertions.output, sample=wildcards.sample, timepoint=timepoints, condition=wildcards.condition),
        ref = rules.download_pombase_data.output.fasta.format(release_version=config["Pombase_release_version"])
    output:
        PBL = f"results/{project_name}/9_concatenated/{{sample}}_{{condition}}.PBL.tsv",
        PBR = f"results/{project_name}/9_concatenated/{{sample}}_{{condition}}.PBR.tsv",
        Reads = f"results/{project_name}/9_concatenated/{{sample}}_{{condition}}.Reads.tsv"
    log:
        "logs/preprocessing/concat_timepoints/{sample}_{condition}.log"
    conda:
        "../envs/biopython.yml"
    params:
        timepoints = " ".join(timepoints)
    shell:
        """
        python workflow/scripts/preprocessing/concat_timepoints.py \
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
        "logs/preprocessing/annotate_insertions/{sample}_{condition}.log"
    conda:
        "../envs/pybedtools.yml"
    message:
        "*** Annotating insertions for {wildcards.sample}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/preprocessing/annotate_insertions.py -i {input.insertions} -g {input.genome_region} -o {output} &> {log}
        """

# Optional: merge the similar time points
# -----------------------------------------------------
rule merge_similar_timepoints:
    input:
        rules.concat_timepoints.output.Reads
    output:
        f"results/{project_name}/11_merged/{{sample}}_{{condition}}.merged.tsv"
    log:
        "logs/preprocessing/merge_similar_timepoints/{sample}_{condition}.log"
    params:
        similar_timepoints = config["similar_timepoints"],
        merged_timepoint = config["merged_timepoint"],
        drop_columns = config["drop_columns"]
    conda:
        "../envs/tabular_operations.yml"
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


# hard filtering the insertions
# -----------------------------------------------------
rule hard_filtering:
    input:
        rules.merge_similar_timepoints.output
    output:
        f"results/{project_name}/12_filtered/{{sample}}_{{condition}}.filtered.tsv"
    log:
        "logs/preprocessing/hard_filtering/{sample}_{condition}.log"
    conda:
        "../envs/tabular_operations.yml"
    params:
        cutoff = config["hard_filtering_cutoff"],
        init_timepoint = config["initial_time_point"]
    message:
        "*** Hard filtering the insertions for {wildcards.sample}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/preprocessing/reads_hard_filtering.py -i {input} -o {output} -c {params.cutoff} -itp {params.init_timepoint} &> {log}
        """

# concat the sample counts and annotations
# -----------------------------------------------------
rule concat_counts_and_annotations:
    input:
        counts = expand(rules.hard_filtering.output, sample=samples, condition=conditions),
        annotations = expand(rules.annotate_insertions.output, sample=samples, condition=conditions)
    output:
        counts = f"results/{project_name}/13_concatenated/raw_reads.tsv",
        annotations = f"results/{project_name}/13_concatenated/annotations.tsv"
    log:
        "logs/preprocessing/concat_counts_and_annotations.log"
    conda:
        "../envs/tabular_operations.yml"
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
        annotations_df = pd.concat(list(annotations_df.values()), axis=0).drop_duplicates()
        print("*** Annotations dataframe shape: ", annotations_df.shape[0])

        counts_df.to_csv(output.counts, sep="\t", index=True)
        annotations_df.to_csv(output.annotations, sep="\t", index=True)

        print("*** Concatenated counts and annotations completed")
        print("*** Counts dataframe shape: ", counts_df.shape[0])
        print("*** Annotations dataframe shape: ", annotations_df.shape[0])

        
        