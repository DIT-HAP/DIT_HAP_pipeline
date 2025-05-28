# Fastp preprocessing for QC and adapter trimming
# -----------------------------------------------------
rule fastp_preprocessing:
    input:
        fq1=lambda wildcards: file_to_rename_dict[wildcards.sample][0],
        fq2=lambda wildcards: file_to_rename_dict[wildcards.sample][1]
    output:
        fq1=f"results/{project_name}/1_fastp/{{sample}}.fastp_1.fq.gz",
        fq2=f"results/{project_name}/1_fastp/{{sample}}.fastp_2.fq.gz",
        html=f"reports/{project_name}/fastp/{{sample}}.fastp.html",
        json=f"reports/{project_name}/fastp/{{sample}}.fastp.json"
    log:
        f"logs/preprocessing/fastp/{{sample}}.log"
    conda:
        "../envs/fastp.yml"
    params:
        adapter_sequence=config["adapter_sequence"],
        adapter_sequence_r2=config["adapter_sequence_r2"]
    threads: 4
    message:
        f"*** Preprocessing fastp for {{input.fq1}} and {{input.fq2}}..."
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
        PBL_r1=f"results/{project_name}/2_demultiplexed/{{sample}}.PBL_1.fq.gz",
        PBL_r2=f"results/{project_name}/2_demultiplexed/{{sample}}.PBL_2.fq.gz",
        PBR_r1=f"results/{project_name}/2_demultiplexed/{{sample}}.PBR_1.fq.gz",
        PBR_r2=f"results/{project_name}/2_demultiplexed/{{sample}}.PBR_2.fq.gz",
        json=f"reports/{project_name}/demultiplexing/{{sample}}.json",
    log:
        "logs/preprocessing/demultiplexing/{sample}.log"
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
        f"*** Demultiplexing {{input.fq1}} and {{input.fq2}}..."
    shell:
        """
        cutadapt --cores {threads} \
                 -g PBL={params.PBL_adapter} \
                 -g PBR={params.PBR_adapter} \
                 -A PBL={params.PBL_reverseComplement_adapter} \
                 -A PBR={params.PBR_reverseComplement_adapter} \
                 -o {params.output_folder}/{wildcards.sample}.{{name}}_1.fq.gz \
                 -p {params.output_folder}/{wildcards.sample}.{{name}}_2.fq.gz \
                 --json {output.json} \
                 {input.fq1} {input.fq2} &> {log}
        """

# fastqc for demultiplexed reads
# -----------------------------------------------------
rule fastqc_demultiplexed:
    input:
        expand(f"results/{project_name}/2_demultiplexed/{{sample}}.{{fragment}}_{{read}}.fq.gz", sample=list(file_to_rename_dict.keys()), fragment=["PBL", "PBR"], read=["1", "2"]),
    output:
        directory(f"reports/{project_name}/fastqc")
    log:
        "logs/preprocessing/fastqc/fastqc_demultiplexed.log"
    conda:
        "../envs/fastqc.yml"
    threads: 16
    message:
        f"*** Running FastQC for demultiplexed reads..."
    shell:
        """
        mkdir -p {output}
        fastqc --threads {threads} {input} -o {output} &> {log}
        """

# BWA-MEM mapping, sorting and indexing
# -----------------------------------------------------
rule bwa_mem_mapping:
    input:
        ref=rules.download_pombase_data.output.fasta.format(release_version=config["Pombase_release_version"]),
        PBL_fq1=rules.demultiplexing.output.PBL_r1,
        PBL_fq2=rules.demultiplexing.output.PBL_r2,
        PBR_fq1=rules.demultiplexing.output.PBR_r1,
        PBR_fq2=rules.demultiplexing.output.PBR_r2
    output:
        PBL_sam=temp(f"results/{project_name}/3_mapped/{{sample}}.PBL.sam"),
        PBR_sam=temp(f"results/{project_name}/3_mapped/{{sample}}.PBR.sam"),
        PBL_sorted=f"results/{project_name}/3_mapped/{{sample}}.PBL.sorted.bam",
        PBR_sorted=f"results/{project_name}/3_mapped/{{sample}}.PBR.sorted.bam",
        PBL_index=f"results/{project_name}/3_mapped/{{sample}}.PBL.sorted.bam.bai",
        PBR_index=f"results/{project_name}/3_mapped/{{sample}}.PBR.sorted.bam.bai",
    log:
        f"logs/preprocessing/bwa_mem_mapping/{{sample}}.log"
    conda:
        "../envs/bwa_mapping.yml"
    threads: 8
    message:
        f"*** Mapping {{input.PBL_fq1}} and {{input.PBL_fq2}} to {{input.ref}}..."
    shell:
        """
        # Command for PBL reads
        # All stdout/stderr from this subshell (including bwa mem errors piped to samtools sort)
        # will be written to the log file, overwriting if it exists.
        bwa mem -t {threads} {input.ref} {input.PBL_fq1} {input.PBL_fq2} > {output.PBL_sam} 2> {log}
        samtools sort -@ {threads} {output.PBL_sam} -O BAM -o {output.PBL_sorted} &>> {log}
        samtools index -@ {threads} {output.PBL_sorted} &>> {log}

        # Add a separator in the log file for clarity
        echo "--- Finished PBL mapping. Starting PBR mapping ---" >> {log} 2>&1

        # Command for PBR reads
        # All stdout/stderr from this subshell will be appended to the log file.
        bwa mem -t {threads} {input.ref} {input.PBR_fq1} {input.PBR_fq2} > {output.PBR_sam} 2>> {log}
        samtools sort -@ {threads} {output.PBR_sam} -O BAM -o {output.PBR_sorted} &>> {log}
        samtools index -@ {threads} {output.PBR_sorted} &>> {log}
        """


# samtools for mapping statistics
# -----------------------------------------------------
rule samtools_mapping_statistics:
    input:
        PBL=rules.bwa_mem_mapping.output.PBL_sorted,
        PBR=rules.bwa_mem_mapping.output.PBR_sorted
    output:
        PBL_stats=f"reports/{project_name}/samtools_mapping_statistics_before_filtering/{{sample}}.PBL.stats.txt",
        PBR_stats=f"reports/{project_name}/samtools_mapping_statistics_before_filtering/{{sample}}.PBR.stats.txt",
        PBL_flagstat=f"reports/{project_name}/samtools_mapping_statistics_before_filtering/{{sample}}.PBL.flagstat.txt",
        PBR_flagstat=f"reports/{project_name}/samtools_mapping_statistics_before_filtering/{{sample}}.PBR.flagstat.txt",
        PBL_idxstats=f"reports/{project_name}/samtools_mapping_statistics_before_filtering/{{sample}}.PBL.idxstats.txt",
        PBR_idxstats=f"reports/{project_name}/samtools_mapping_statistics_before_filtering/{{sample}}.PBR.idxstats.txt"
    log:
        "logs/preprocessing/samtools_mapping_statistics_before_filtering/{sample}.log"
    conda:
        "../envs/samtools.yml"
    threads: 1
    message:
        f"*** Running Samtools for mapping statistics..."
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

# BAM TAG value extraction
# -----------------------------------------------------
rule bam_tag_value_extraction:
    input:
        PBL=rules.bwa_mem_mapping.output.PBL_sorted,
        PBR=rules.bwa_mem_mapping.output.PBR_sorted
    output:
        PBL_tag_value=f"results/{project_name}/4_filtered/{{sample}}.PBL.tag_value.txt",
        PBR_tag_value=f"results/{project_name}/4_filtered/{{sample}}.PBR.tag_value.txt"
    log:
        "logs/preprocessing/bam_tag_value_extraction/{sample}.log"
    conda:
        "../envs/samtools.yml"
    threads: 4
    message:
        "*** Extracting BAM tag values..."
    shell:
        """
        samtools view -H {input.PBL} | grep -E '^@PG' > {output.PBL_tag_value}
        samtools view -H {input.PBR} | grep -E '^@PG' > {output.PBR_tag_value}
        """

# Samtools filtering
# Filtering criteria for single-end reads:
# 1. Properly mapped reads (flag==0 or flag==16)
# 2. No more than 3 mismatches (NM<=3)
# 3. No secondary alignments (XA and SA are empty)
# 4. No supplementary alignments (ncigar==1)
# 5. Quality score >= 1
# Filtering criteria for paired-end reads:
# 1. Properly mapped reads -f 2
# 2. No more than 3 mismatches -e "[NM] <= 3"
# 3. No secondary alignments -e "!defined(XA) && !defined(SA)"
# 4. No supplementary alignments -e "ncigar == 1"
# 5. Quality score >= 1
# -----------------------------------------------------
rule samtools_filtering:
    input:
        PBL=rules.bwa_mem_mapping.output.PBL_sorted,
        PBR=rules.bwa_mem_mapping.output.PBR_sorted
    output:
        PBL=f"results/{project_name}/4_filtered/{{sample}}.PBL.filtered.bam",
        PBR=f"results/{project_name}/4_filtered/{{sample}}.PBR.filtered.bam"
    log:
        "logs/preprocessing/samtools_filtering/{sample}.log"
    conda:
        "../envs/samtools.yml"
    params:
        filter = "([NM]&&[NM]<=3)&&(!(([XA])||([SA])))&&(ncigar==1)"
    threads: 4
    message:
        "*** Filtering {input.PBL} and {input.PBR}..."
    shell:
        """
        echo "*** Filtering PBL reads..." > {log}
        samtools view -@ {threads} -h -q 1 -f 2 \
            -e '{params.filter}' \
            -bo {output.PBL} {input.PBL} &>> {log}

        echo "*** Filtering PBR reads..." >> {log}
        samtools view -@ {threads} -h -q 1 -f 2 \
            -e '{params.filter}' \
            -bo {output.PBR} {input.PBR} &>> {log}
        """

# Samtools sorting and indexing for BAM files
# -----------------------------------------------------
rule samtools_sorting_and_indexing:
    input:
        PBL=rules.samtools_filtering.output.PBL,
        PBR=rules.samtools_filtering.output.PBR
    output:
        PBL_sorted=f"results/{project_name}/5_sorted/{{sample}}.PBL.sorted.bam",
        PBR_sorted=f"results/{project_name}/5_sorted/{{sample}}.PBR.sorted.bam",
        PBL_index=f"results/{project_name}/5_sorted/{{sample}}.PBL.sorted.bam.bai",
        PBR_index=f"results/{project_name}/5_sorted/{{sample}}.PBR.sorted.bam.bai"
    log:
        "logs/preprocessing/samtools_sorting_and_indexing/{sample}.log"
    conda:
        "../envs/samtools.yml"
    threads: 4
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

use rule samtools_mapping_statistics as samtools_mapping_statistics_after_filtering with:
    input:
        PBL=rules.samtools_sorting_and_indexing.output.PBL_sorted,
        PBR=rules.samtools_sorting_and_indexing.output.PBR_sorted
    output:
        PBL_stats=f"reports/{project_name}/samtools_mapping_statistics_after_filtering/{{sample}}.PBL.stats.txt",
        PBR_stats=f"reports/{project_name}/samtools_mapping_statistics_after_filtering/{{sample}}.PBR.stats.txt",
        PBL_flagstat=f"reports/{project_name}/samtools_mapping_statistics_after_filtering/{{sample}}.PBL.flagstat.txt",
        PBR_flagstat=f"reports/{project_name}/samtools_mapping_statistics_after_filtering/{{sample}}.PBR.flagstat.txt",
        PBL_idxstats=f"reports/{project_name}/samtools_mapping_statistics_after_filtering/{{sample}}.PBL.idxstats.txt",
        PBR_idxstats=f"reports/{project_name}/samtools_mapping_statistics_after_filtering/{{sample}}.PBR.idxstats.txt"
    log:
        "logs/preprocessing/samtools_mapping_statistics_after_filtering/{sample}.log"

# MultiQC for preprocessing reports
# -----------------------------------------------------
rule multiqc_preprocessing:
    input:
        fastp_json=expand(rules.fastp_preprocessing.output.json, sample=list(file_to_rename_dict.keys())),
        demux_json=expand(rules.demultiplexing.output.json, sample=list(file_to_rename_dict.keys())),
        fastqc_dir=f"reports/{project_name}/fastqc",
    output:
        report(f"reports/{project_name}/multiqc/preprocessing_multiqc_report.html")
    log:
        f"logs/preprocessing/multiqc_preprocessing.log"
    conda:
        "../envs/multiqc.yml"
    params:
        outdir=f"reports/{project_name}/multiqc",
        title="Preprocessing Quality Control Report",
        cl_config=config["multiqc_config"]
    message:
        "*** Generating MultiQC report for preprocessing steps..."
    shell:
        """
        multiqc \
            --title "{params.title}" \
            --filename preprocessing_multiqc_report.html \
            --cl-config {params.cl_config} \
            --outdir {params.outdir} \
            --force \
            --verbose \
            --fn_as_s_name \
            {input.fastp_json} \
            {input.demux_json} \
            &> {log}
        """


# rule extract_insertions:
#     input:
#         rules.samtools_filtering.output
#     output:
#         project_dir/"5_insertions/{sample}/{sample}-{timepoint}.{PBtype}.bed"
#     log:
#         f"logs/{project_name}/extract_insertions/{{sample}}/{{sample}}-{{timepoint}}.{{PBtype}}.log"
#     threads: 8
#     shell:
#         """
#         python workflow/scripts/preprocessing/extract_insertions.py -i {input} -o {output} -n {threads} > {log}
#         """

# # rule cp_insertions:
# #     input:
# #         PBL=input_dir/"{sample}/{sample}-{timepoint}.PBL.bed",
# #         PBR=input_dir/"{sample}/{sample}-{timepoint}.PBR.bed"
# #     output:
# #         PBL=project_dir/"5_insertions/{sample}/{sample}-{timepoint}.PBL.bed",
# #         PBR=project_dir/"5_insertions/{sample}/{sample}-{timepoint}.PBR.bed"
# #     shell:
# #         """
# #         cp {input.PBL} {output.PBL}
# #         cp {input.PBR} {output.PBR}
# #         """

# rule merge_insertions:
#     input:
#         PBL_insertions=project_dir/"5_insertions/{sample}/{sample}-{timepoint}.PBL.bed",
#         PBR_insertions=project_dir/"5_insertions/{sample}/{sample}-{timepoint}.PBR.bed"
#     output:
#         project_dir/"6_merged/{sample}/{sample}-{timepoint}.bed"
#     log:
#         f"logs/{project_name}/merge_insertions/{{sample}}/{{sample}}-{{timepoint}}.log"
#     shell:
#         """
#         python workflow/scripts/preprocessing/merge_insertions.py \
#             -il {input.PBL_insertions} \
#             -ir {input.PBR_insertions} \
#             -o {output} > {log}
#         """

# rule concat_timepoints:
#     input:
#         tp_files = lambda wildcards: expand(project_dir/"6_merged/{sample}/{sample}-{timepoint}.bed", sample=wildcards.sample, timepoint=timepoints),
#         ref = genome_reference
#     output:
#         PBL=project_dir/"7_concated/{sample}/{sample}.PBL",
#         PBR=project_dir/"7_concated/{sample}/{sample}.PBR",
#         Reads=project_dir/"7_concated/{sample}/{sample}.Reads"
#     log:
#         f"logs/{project_name}/concat_timepoints/{{sample}}.log"
#     params:
#         timepoints = timepoints
#     shell:
#         """
#         python workflow/scripts/preprocessing/concat_timepoints.py \
#             -s {wildcards.sample} \
#             -i {input.tp_files} \
#             -tp {params.timepoints} \
#             -g {input.ref} \
#             -ol {output.PBL} \
#             -or {output.PBR} \
#             -o {output.Reads} > {log}
#         """

# rule compare_PBL_PBR:
#     input:
#         PBL=rules.concat_timepoints.output.PBL,
#         PBR=rules.concat_timepoints.output.PBR,
#         Reads=rules.concat_timepoints.output.Reads
#     output:
#         report("reports/compare_PBL_PBR/{sample}.pdf")
#     log:
#         "logs/compare_PBL_PBR/{sample}.log"
#     shell:
#         """
#         python src/compare_PBL_PBR_reads.py -pbl {input.PBL} -pbr {input.PBR} -reads {input.Reads} -pdf {output} > {log}
#         """