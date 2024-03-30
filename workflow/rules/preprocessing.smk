
# rule fastp_preprocessing:
#     input:
#         input_dir+"/{sample}-{timepoint}_1.fq.gz"
#     output:
#         fq="data/1_trimmed/{sample}/{sample}-{timepoint}.fq.gz",
#         html="reports/fastp/{sample}/{sample}-{timepoint}.fastp.html",
#         json="reports/fastp/{sample}/{sample}-{timepoint}.fastp.json"
#     log:
#         "logs/fastp/{sample}/{sample}-{timepoint}.log"
#     threads: 8
#     shell:
#         """
#         fastp --adapter_sequence CTGTCTCTTATACACATCT \
#               --disable_quality_filtering \
#               --disable_length_filtering \
#               --overrepresentation_analysis \
#               -j {output.json} \
#               -h {output.html} \
#               --thread {threads} \
#               --in1 {input} \
#               -o {output.fq} > {log}
#         """

# rule dissect_PBL_PBR:
#     input:
#         rules.fastp_preprocessing.output.fq
#     output:
#         outPBL="data/2_dissected/{sample}/{sample}-{timepoint}.PBL.fq.gz",
#         outPBR="data/2_dissected/{sample}/{sample}-{timepoint}.PBR.fq.gz",
#         summary="reports/dissect_PBL_PBR/{sample}/{sample}-{timepoint}_reads_summary.csv"
#     log:
#         "logs/dissect_PBL_PBR/{sample}/{sample}-{timepoint}.log"
#     shell:
#         "python src/dissect_PBL_and_PBR.py -i {input} -ol {output.outPBL} -or {output.outPBR} -s {output.summary} > {log}"

# rule merge_read_summary:
#     # input:
#         # expand("reports/dissect_PBL_PBR/{sample}/{sample}-{timepoint}_reads_summary.csv", sample=config["samples"], timepoint=config["timepoints"])
#         # expand(rules.dissect_PBL_PBR.output.summary, sample=config["samples"], timepoint=config["timepoints"])
#     params:
#         input_dir = "reports/dissect_PBL_PBR"
#     output:
#         report("reports/dissect_PBL_PBR/reads_summary.csv")
#     run:
#         import pandas as pd
#         from pathlib import Path
#         input_dir = Path(params.input_dir)
#         input_files = input_dir.glob("*/*_reads_summary.csv")
#         df = pd.concat(
#             [pd.read_csv(f, sep="\t") for f in input_files], 
#             axis=0
#         )
#         df.sort_values(by=df.columns[0], inplace=True)
#         df.to_csv(output[0], index=False)

# rule bwa_mem_mapping:
#     input:
#         ref=genome_reference,
#         fqIn="data/2_dissected/{sample}/{sample}-{timepoint}.{PBtype}.fq.gz"
#     output:
#         "data/3_mapped/{sample}/{sample}-{timepoint}.{PBtype}.sam"
#     log:
#         "logs/bwa_mem_mapping/{sample}/{sample}-{timepoint}.{PBtype}.log"
#     threads: 8
#     shell:
#         """
#         bwa mem -t {threads} {input.ref} {input.fqIn} \
#         | samtools sort -@ {threads} -O SAM -o {output} - > {log}
#         """

# rule samtools_filtering:
#     input:
#         rules.bwa_mem_mapping.output
#     output:
#         "data/4_filtered/{sample}/{sample}-{timepoint}.{PBtype}.unique.sam"
#     log:
#         "logs/samtools_filtering/{sample}/{sample}-{timepoint}.{PBtype}.log"
#     params:
#         filter="(flag==0||flag==16)&&([NM]&&[NM]<=3)&&(!(([XA])||([SA])))&&(ncigar==1)"
#     threads: 8
#     shell:
#         """
#         samtools view -@ {threads} -h -q 1 -e "{params.filter}" {input} > {output} 2> {log}
#         """

# rule extract_insertions:
#     input:
#         rules.samtools_filtering.output
#     output:
#         "data/5_insertions/{sample}/{sample}-{timepoint}.{PBtype}.bed"

#     log:
#         "logs/extract_insertions/{sample}/{sample}-{timepoint}.{PBtype}.log"
#     threads: 8
#     shell:
#         """
#         python src/extract_insertions.py -i {input} -o {output} -n {threads} > {log}
#         """

rule cp_insertions:
    input:
        PBL=input_dir/"{sample}/{sample}-{timepoint}.PBL.bed",
        PBR=input_dir/"{sample}/{sample}-{timepoint}.PBR.bed"
    output:
        PBL=project_dir/"5_insertions/{sample}/{sample}-{timepoint}.PBL.bed",
        PBR=project_dir/"5_insertions/{sample}/{sample}-{timepoint}.PBR.bed"
    shell:
        """
        cp {input.PBL} {output.PBL}
        cp {input.PBR} {output.PBR}
        """

rule merge_insertions:
    input:
        PBL_insertions=rules.cp_insertions.output.PBL,
        PBR_insertions=rules.cp_insertions.output.PBR
    output:
        project_dir/"6_merged/{sample}/{sample}-{timepoint}.bed"
    log:
        "workflow/logs/merge_insertions/{sample}/{sample}-{timepoint}.log"
    shell:
        """
        python workflow/scripts/preprocessing/merge_insertions.py \
            -il {input.PBL_insertions} \
            -ir {input.PBR_insertions} \
            -o {output} > {log}
        """

rule concat_timepoints:
    input:
        tp_files = lambda wildcards: expand(project_dir/"6_merged/{sample}/{sample}-{timepoint}.bed", sample=wildcards.sample, timepoint=timepoints),
        ref = genome_reference
    output:
        PBL=project_dir/"7_concated/{sample}/{sample}.PBL",
        PBR=project_dir/"7_concated/{sample}/{sample}.PBR",
        Reads=project_dir/"7_concated/{sample}/{sample}.Reads"
    log:
        "workflow/logs/concat_timepoints/{sample}.log"
    params:
        timepoints = timepoints
    shell:
        """
        python workflow/scripts/preprocessing/concat_timepoints.py \
            -s {wildcards.sample} \
            -i {input.tp_files} \
            -tp {params.timepoints} \
            -g {input.ref} \
            -ol {output.PBL} \
            -or {output.PBR} \
            -o {output.Reads} > {log}
        """

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