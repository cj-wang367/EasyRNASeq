#!/usr/bin/env Rscript
library(optparse)

# 1. 命令行参数解析
option_list <- list(
  make_option(c("-v", "--vsd_file"), type="character", default=NULL, help="DESeq2输出的标准化矩阵 (DESeq2_vsd_matrix.csv)"),
  make_option(c("-m", "--meta_file"), type="character", default=NULL, help="样本元信息表"),
  make_option(c("-o", "--out_dir"), type="character", default="WGCNA_Output", help="输出目录")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$vsd_file) || is.null(opt$meta_file)){
  print_help(opt_parser)
  stop("缺少必须参数！", call.=FALSE)
}

# 2. 载入核心包
suppressMessages({
  library(WGCNA)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

enableWGCNAThreads()
options(stringsAsFactors = FALSE)

# 3. 数据预处理与高变基因过滤
message("=> [1/6] 正在读取 VST 标准化矩阵与样本信息...")
vsd <- read_csv(opt$vsd_file, show_col_types = FALSE) %>% as.data.frame()
rownames(vsd) <- vsd$id
vsd$id <- NULL

meta <- read_tsv(opt$meta_file, show_col_types = FALSE) %>% as.data.frame()

# WGCNA 要求矩阵的行是样本，列是基因
datExpr <- as.data.frame(t(vsd))

# 匹配样本顺序
meta <- meta[match(rownames(datExpr), meta$SampleID), ]
if(!all(rownames(datExpr) == meta$SampleID)) stop("元信息表与矩阵样本名不匹配！")

message("=> [2/6] 正在通过 MAD (中位数绝对偏差) 过滤高变基因...")
# WGCNA 官方建议过滤低表达/无变化的基因以降低噪音。这里提取 top 10000 基因。
m.mad <- apply(datExpr, 2, mad)
datExpr <- datExpr[, which(m.mad > 0)] # 剔除绝对无变化的基因
if(ncol(datExpr) > 10000) {
  datExpr <- datExpr[, order(m.mad, decreasing=TRUE)[1:10000]]
}
message(sprintf("=> 保留了 %d 个高变基因用于网络构建。", ncol(datExpr)))

# 检查是否存在缺失值或离群样本
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 4. 自动确定软阈值 (Soft Threshold)
message("=> [3/6] 正在自动扫描最优 Soft Threshold Power...")
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

# 自动选取 R^2 达到 0.85 的最小 power
power <- sft$powerEstimate
if(is.na(power)) {
  message("  -> 警告: 自动选取失败(可能由于样本太少)，智能回退至默认阈值 Power=6")
  power <- 6
} else {
  message(sprintf("  -> 成功选取最佳软阈值: Power = %d", power))
}

# 绘制软阈值图
pdf(file.path(opt$out_dir, "Soft_Threshold_Plot.pdf"), width=9, height=5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()

# 5. 一步法构建共表达网络与模块划分
message("=> [4/6] 正在构建共表达网络与模块聚类 (Blockwise Modules)...")
net <- blockwiseModules(datExpr, power = power,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamStage = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = FALSE, verbose = 3)

# 提取模块颜色
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# 绘制基因树与模块颜色分配图
pdf(file.path(opt$out_dir, "Gene_Dendrogram_and_Colors.pdf"), width=12, height=9)
plotDendroAndColors(geneTree, moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 6. 计算模块与性状的相关性 (Module-Trait Relationships)
message("=>[5/6] 正在计算模块与性状(分组)的相关性...")

trait_tmp <- meta %>% select(-SampleID)
traitData <- as.data.frame(model.matrix(~ . - 1, data = trait_tmp))
colnames(traitData) <- gsub("condition|batch", "", colnames(traitData))

# 重新排序特征向量 MEs
MEs <- orderMEs(MEs)

# 计算相关系数矩阵和 P 值矩阵
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 绘制模块-性状关联热图
message("=> [6/6] 正在生成模块-性状热图...")
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file.path(opt$out_dir, "Module_Trait_Heatmap.pdf"), width=8, height=8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# 7. 导出基因-模块映射表
gene_module_df <- data.frame(
  id = names(datExpr),
  ModuleColor = moduleColors
) %>% arrange(ModuleColor)

write_csv(gene_module_df, file.path(opt$out_dir, "WGCNA_Gene_Modules.csv"))
message("🎉 WGCNA 共表达网络分析全部完成！结果保存在 ", opt$out_dir)
