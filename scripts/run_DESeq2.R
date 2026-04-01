#!/usr/bin/env Rscript
library(optparse)

# 1. 命令行参数解析 
option_list <- list(
  make_option(c("-s", "--salmon_dir"), type="character", default=NULL, help="Salmon定量目录"),
  make_option(c("-m", "--meta_file"), type="character", default=NULL, help="样本元信息表"),
  make_option(c("-g", "--gtf_file"), type="character", default=NULL, help="GTF注释文件路径"),
  make_option(c("-o", "--out_dir"), type="character", default="DESeq2_output", help="输出目录"),
  make_option(c("-r", "--ref_group"), type="character", default=NULL, help="对照组的名称 (如 Control)")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$salmon_dir) || is.null(opt$meta_file) || is.null(opt$gtf_file) || is.null(opt$ref_group)){
  print_help(opt_parser)
  stop("缺少必须参数！", call.=FALSE)
}

# 2. 载入核心包
suppressMessages({
  library(tximport)
  library(rtracklayer)
  library(DESeq2)
  library(apeglm)
  library(vsn)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
})

if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)
set.seed(1)

# 3. 读取元数据与解析 GTF
message("=> 正在读取元信息与 GTF 文件...")
coldata <- read_tsv(opt$meta_file, show_col_types = FALSE) |> as.data.frame()
rownames(coldata) <- coldata$SampleID

coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = opt$ref_group)

gtf <- import(opt$gtf_file)
tx2gene <- as.data.frame(mcols(gtf)[, c("transcript_id", "gene_id")])
tx2gene <- na.omit(unique(tx2gene))

# 4. tximport 导入数据
message("=> 正在导入 Salmon 定量结果...")
files <- file.path(opt$salmon_dir, coldata$SampleID, "quant.sf")
names(files) <- coldata$SampleID
if(any(!file.exists(files))) stop("某些样本的 quant.sf 文件不存在！")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

# 5. DESeq2 核心计算与差异表达分析
message("=> 正在执行 DESeq2 差异分析...")
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# 初步过滤 (保留总 count >= 10 的基因)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
message(sprintf("=> 过滤低表达基因后，保留 %d 个基因", sum(keep)))

# 运行拟合
dds <- DESeq(dds, parallel = FALSE)

# 获取收缩后的 log2FC (apeglm)
res_names <- resultsNames(dds)
target_coef <- res_names[grep("condition_", res_names)][1]
message("=> 正在计算贝叶斯收缩系数: ", target_coef)
res_shrunk <- lfcShrink(dds, coef = target_coef, type = "apeglm")

# 6. 构建总表 (Master Table) 供下游模块使用
message("=> 正在构建基础大总表 (包含 Count, TPM, Means, Stats)...")

# (A) 提取 Raw Counts 并改列名
raw_counts <- as.data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste0(colnames(raw_counts), "_count")

# (B) 提取 TPM 并改列名 (只保留在 dds 中过滤后剩下的基因)
tpm_matrix <- as.data.frame(txi$abundance)[rownames(dds), ]
colnames(tpm_matrix) <- paste0(colnames(tpm_matrix), "_tpm")

# (C) 提取差异表达统计量
stats_df <- as.data.frame(res_shrunk) |>
  select(log2FoldChange, pvalue, padj) |>
  rename(log2FC = log2FoldChange, PValue = pvalue, FDR = padj)

# (D) 计算各组 TPM 均值
groups <- levels(coldata$condition)
mean_tpm_list <- list()
for (grp in groups) {
  samples_in_grp <- coldata$SampleID[coldata$condition == grp]
  tpm_cols <- paste0(samples_in_grp, "_tpm")
  mean_tpm_list[[paste0(grp, "_mean")]] <- rowMeans(tpm_matrix[, tpm_cols, drop=FALSE])
}
mean_tpm_df <- as.data.frame(mean_tpm_list)

# (E) 组装并标记显著性
master_table <- cbind(raw_counts, tpm_matrix, mean_tpm_df, stats_df) |>
  rownames_to_column(var = "id") |>
  mutate(Significance = case_when(
    FDR < 0.05 & log2FC > 1  ~ "Up",
    FDR < 0.05 & log2FC < -1 ~ "Down",
    TRUE ~ "Not_Sig"
  ))

# 导出供下游 Enrichment 使用
write_csv(master_table, file.path(opt$out_dir, "Master_DEG_All_Genes.csv"))

sig_table <- master_table |> filter(Significance != "Not_Sig")
write_csv(sig_table, file.path(opt$out_dir, "Master_DEG_Significant.csv"))

# 7. 数据转换绘图 & 导出 WGCNA 接口矩阵
message("=> 正在绘制质控图表并导出 WGCNA 矩阵...")
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
assay(vsd)[assay(vsd) < 0] <- 0 

# 导出给 WGCNA 用的标准化矩阵 
vsd_matrix <- as.data.frame(assay(vsd)) |> rownames_to_column(var = "id")
write_csv(vsd_matrix, file.path(opt$out_dir, "DESeq2_vsd_matrix.csv"))

# 火山图
p_vol <- ggplot(master_table, aes(x = log2FC, y = -log10(PValue), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not_Sig" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_bw() + labs(title = "Volcano plot (shrunken LFC)")
ggsave(file.path(opt$out_dir, "Volcano_plot.pdf"), p_vol, width=6, height=5)

# PCA 图
p_pca <- plotPCA(vsd, intgroup = "condition") + 
         ggtitle("PCA on VST-transformed counts") + theme_bw()
ggsave(file.path(opt$out_dir, "PCA_plot.pdf"), p_pca, width=6, height=5)

# 样本距离热图
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pdf(file.path(opt$out_dir, "Sample_Distances_Heatmap.pdf"), width=6, height=5)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distances")
dev.off()

message("核心差异分析模块已全部完成！输出目录: ", opt$out_dir)
