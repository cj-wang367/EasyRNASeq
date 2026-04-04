# EasyRNASeq

![alt text](https://img.shields.io/badge/License-MIT-green.svg)
![alt text](https://img.shields.io/badge/Snakemake-≥7.0-blue.svg)
![alt text](https://img.shields.io/badge/Python-3.8+-blue.svg)
![alt text](https://img.shields.io/badge/R-4.1+-blue.svg)

EasyRNASeq 是一个基于 Snakemake 开发的 RNA-Seq 一键式自动化分析工具。输入Raw Data 自动获得差异基因分析、GO/KEGG/GSEA 富集分析、WGCNA 共表达网络分析结果。

## :arrow_double_down:安装
```
# 0.安装Snakemake环境
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
conda activate snakemake

# 1. 克隆本仓库到本地服务器
git clone https://github.com/cj-wang367/EasyRNASeq.git

# 2. 进入项目目录
cd EasyRNASeq

# 3. 为主程序和测试脚本赋予执行权限
chmod +x EasyRNASeq run_test.sh

# 4.快速测试
./run_test.sh
```

## :play_or_pause_button:用法
```
# 以测试数据为例：
# 1. 解压压缩的参考文件
cd test_data
gzip -d -k genomic.gtf.gz rna.fna.gz
cd ..

# 2. 运行基础差异分析流程
./EasyRNASeq \
  -input_dir test_data \
  -reference test_data/rna.fna \
  -gtf test_data/genomic.gtf \
  -meta test_data/metadata.txt \
  -ref_group Control \
  -out_dir test_result \
  -threads 4
```
### 完整参数说明

**必需参数:**

- `-input_dir` : 存放双端测序数据 (命名规则：*R1.fastq.gz, R2.fastq.gz* 或 *R1.fq.gz, R2.fq.gz*) 的目录路径。
    
- `-reference` : 转录组参考序列文件 (.fasta / .fna)。
    
- `-gtf` : 参考基因组注释文件 (.gtf)。
    
- `-meta` : 样本元信息表（制表符 \t 分隔，必须包含 SampleID 和 condition 两列）。
    
- `-ref_group` : condition 列中作为对照组的名称（例如 Control）。
    

**基础参数:**

- `-out_dir` : 结果输出主目录（默认: ./EasyRNASeq_Output）。
    
- `-threads` : 最大并行线程数（默认: 8）。
    

**可选分析模块 :**

- `-run_wgcna` : 开启 WGCNA 模块，自动构建共表达网络与性状热图。
    
- `-run_enrichment` : 开启功能富集模块（含 GO/KEGG 及 GSEA 分析）。
    
- `-pep` : 蛋白质序列文件 (.faa)，富集模块**必填**（用于 eggNOG 同源注释）。
    
- `-eggnog_db` : 本地 eggNOG 数据库的绝对路径（[eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)）。

### 输出结构
```
Output_Directory/
├── trimmed/                  # Fastp 质控过滤后的 Clean data
├── salmon_quant/             # Salmon 样本定量结果
├── MultiQC/                  # 全局质控汇总报告 (HTML)
├── DESeq2_Analysis/          # 差异表达分析结果 (火山图, PCA图, DEG总表)
├── Annotation/               # [可选] eggNOG 功能注释表
├── Enrichment/               # [可选] GO/KEGG/GSEA 富集分析图表
└── WGCNA/                    # [可选] WGCNA 共表达网络结果及关联热图
```
