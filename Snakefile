import os
import pandas as pd

PIPELINE_DIR = workflow.basedir

def get_fastq_input(wildcards):
    sample = wildcards.sample
    for ext in ['fastq.gz', 'fq.gz']:
        r1 = os.path.join(INPUT_DIR, f"{sample}_R1.{ext}")
        r2 = os.path.join(INPUT_DIR, f"{sample}_R2.{ext}")
        if os.path.exists(r1) and os.path.exists(r2):
            return {'r1': r1, 'r2': r2}
    raise FileNotFoundError(f"未找到样本 {sample} 的 R1/R2 文件，需要 .fastq.gz 或 .fq.gz")

INPUT_DIR = config["input_dir"]
OUT_DIR   = config["out_dir"]
REF_FASTA = config["reference"]
META_FILE = config["metadata"]

samples_df = pd.read_csv(META_FILE, sep='\t')
SAMPLE_NAMES = samples_df['SampleID'].astype(str).tolist()

ENV_UPSTREAM = os.path.join(PIPELINE_DIR, "envs/upstream.yaml")
ENV_DOWNSTREAM = os.path.join(PIPELINE_DIR, "envs/downstream.yaml")
ENV_EGGNOG = os.path.join(PIPELINE_DIR, "envs/eggnog.yaml")

TARGETS =[
    expand(os.path.join(OUT_DIR, "salmon_quant", "{sample}", "quant.sf"), sample=SAMPLE_NAMES),
    os.path.join(OUT_DIR, "MultiQC", "multiqc_report.html"),
    os.path.join(OUT_DIR, "DESeq2_Analysis", "Master_DEG_Significant.csv")
]

# 如果用户开启了富集分析、WGCNA，增加输出目标
if config.get("run_enrichment"):
    TARGETS.append(os.path.join(OUT_DIR, "Enrichment", "Final_Master_Table_Annotated.csv"))
if config.get("run_wgcna"):
    TARGETS.append(os.path.join(OUT_DIR, "WGCNA", "Module_Trait_Heatmap.pdf"))


rule all:
    input:TARGETS

# 1. 建立 Salmon 转录组索引
rule salmon_index:
    input:
        fasta = REF_FASTA
    output:
        directory(os.path.join(OUT_DIR, "reference", "salmon_index"))
    conda: 
        ENV_UPSTREAM
    threads: 
        config.get("threads", 8)
    shell:
        """
        salmon index -t {input.fasta} -i {output} -p {threads}
        """

# 2. 质控与接头过滤 (Fastp) 
rule fastp:
    input:
        unpack(get_fastq_input)   # 动态获取文件
    output:
        r1_clean = os.path.join(OUT_DIR, "trimmed", "{sample}_R1_clean.fq.gz"),
        r2_clean = os.path.join(OUT_DIR, "trimmed", "{sample}_R2_clean.fq.gz"),
        html = os.path.join(OUT_DIR, "trimmed", "{sample}_fastp.html"),
        json = os.path.join(OUT_DIR, "trimmed", "{sample}_fastp.json")
    conda: ENV_UPSTREAM
    threads: 4
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1_clean} -O {output.r2_clean} \
              -w {threads} \
              -h {output.html} -j {output.json}
        """

# 3. 表达量定量 (Salmon Quant) 
rule salmon_quant:
    input:
        index = os.path.join(OUT_DIR, "reference", "salmon_index"),
        # 衔接 fastp 产生的 clean data
        r1 = os.path.join(OUT_DIR, "trimmed", "{sample}_R1_clean.fq.gz"),
        r2 = os.path.join(OUT_DIR, "trimmed", "{sample}_R2_clean.fq.gz")
    output:
        quant = os.path.join(OUT_DIR, "salmon_quant", "{sample}", "quant.sf")
    params:
        outdir = os.path.join(OUT_DIR, "salmon_quant", "{sample}")
    conda: 
        ENV_UPSTREAM
    threads: 
        config.get("threads", 8)
    shell:
        """
        salmon quant -i {input.index} -l A \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} --validateMappings -o {params.outdir}
        """

# 4. 汇总质控报告 (MultiQC) 
rule multiqc:
    input:
        quant = expand(os.path.join(OUT_DIR, "salmon_quant", "{sample}", "quant.sf"), sample=SAMPLE_NAMES),
        fastp_json = expand(os.path.join(OUT_DIR, "trimmed", "{sample}_fastp.json"), sample=SAMPLE_NAMES)
    output:
        html = os.path.join(OUT_DIR, "MultiQC", "multiqc_report.html")
    params:
        search_dir = OUT_DIR,
        out_dir = os.path.join(OUT_DIR, "MultiQC")
    conda: 
        ENV_UPSTREAM
    shell:
        """
        multiqc {params.search_dir} -f -o {params.out_dir}
        """

# 5. 下游差异表达分析 (DESeq2)
rule deseq2:
    input:
        quant = expand(os.path.join(OUT_DIR, "salmon_quant", "{sample}", "quant.sf"), sample=SAMPLE_NAMES),
        meta = config["metadata"],
        gtf = config["gtf"]
    output:
        all_csv = os.path.join(OUT_DIR, "DESeq2_Analysis", "Master_DEG_All_Genes.csv"),
        sig_csv = os.path.join(OUT_DIR, "DESeq2_Analysis", "Master_DEG_Significant.csv"),
        vsd_csv = os.path.join(OUT_DIR, "DESeq2_Analysis", "DESeq2_vsd_matrix.csv")
    params:
        salmon_dir = os.path.join(OUT_DIR, "salmon_quant"),
        out_dir = os.path.join(OUT_DIR, "DESeq2_Analysis"),
        ref_group = config["ref_group"],
        r_script = os.path.join(PIPELINE_DIR, "scripts/run_DESeq2.R")
    conda:
        ENV_DOWNSTREAM
    shell:
        """
        Rscript {params.r_script} \
            --salmon_dir {params.salmon_dir} \
            --meta_file {input.meta} \
            --gtf_file {input.gtf} \
            --out_dir {params.out_dir} \
            --ref_group {params.ref_group}
        """

# 6. eggNOG 注释 (仅在 run_enrichment=True 时执行)
rule eggnog_anno:
    input:
        pep = config.get("pep_file")
    output:
        anno = os.path.join(OUT_DIR, "Annotation", "eggnog.emapper.annotations")
    params:
        db = config.get("eggnog_db"),
        out_prefix = os.path.join(OUT_DIR, "Annotation", "eggnog")
    conda: ENV_EGGNOG
    threads: config.get("threads", 8)
    shell:
        """
        emapper.py -i {input.pep} -o {params.out_prefix} \
            --data_dir {params.db} \
            --cpu {threads} -m diamond
        """

# 7. 功能富集与 GSEA
rule enrichment:
    input:
        deg_res = os.path.join(OUT_DIR, "DESeq2_Analysis", "Master_DEG_All_Genes.csv"),
        anno = os.path.join(OUT_DIR, "Annotation", "eggnog.emapper.annotations"),
        gtf = config["gtf"]  
    output:
        final_table = os.path.join(OUT_DIR, "Enrichment", "Final_Master_Table_Annotated.csv"),
        gsea_plot = os.path.join(OUT_DIR, "Enrichment", "GSEA_Plot.pdf")
    params:
        out_dir = os.path.join(OUT_DIR, "Enrichment"),
        r_script = os.path.join(PIPELINE_DIR, "scripts/run_enrichment.R")
    conda: ENV_DOWNSTREAM
    shell:
        """
        Rscript {params.r_script} \
            --deg_file {input.deg_res} \
            --anno_file {input.anno} \
            --gtf_file {input.gtf} \
            --out_dir {params.out_dir}
        """

# 8. WGCNA 分析 
rule wgcna:
    input:
        vsd_mat = os.path.join(OUT_DIR, "DESeq2_Analysis", "DESeq2_vsd_matrix.csv"),
        meta = config["metadata"]
    output:
        os.path.join(OUT_DIR, "WGCNA", "Module_Trait_Heatmap.pdf")
    params:
        out_dir = os.path.join(OUT_DIR, "WGCNA"),
        r_script = os.path.join(PIPELINE_DIR, "scripts/run_WGCNA.R")
    conda: ENV_DOWNSTREAM
    shell:
        """
        Rscript {params.r_script} \
            --vsd_file {input.vsd_mat} \
            --meta_file {input.meta} \
            --out_dir {params.out_dir}
        """
