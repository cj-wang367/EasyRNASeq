#!/usr/bin/env Rscript
library(optparse)

# 1. 命令行参数解析
option_list <- list(
  make_option(c("-d", "--deg_file"), type="character", default=NULL, help="DESeq2 输出的大总表 (Master_DEG_All_Genes.csv)"),
  make_option(c("-a", "--anno_file"), type="character", default=NULL, help="eggNOG 输出的注释表 (*.emapper.annotations)"),
  make_option(c("-g", "--gtf_file"), type="character", default=NULL, help="GTF 文件 (用于桥接 protein_id 和 gene_id)"),
  make_option(c("-o", "--out_dir"), type="character", default="Enrichment_Output", help="输出目录")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$deg_file) || is.null(opt$anno_file) || is.null(opt$gtf_file)){
  print_help(opt_parser)
  stop("缺少必须参数！", call.=FALSE)
}

# 2. 载入核心包
suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rtracklayer)
  library(clusterProfiler)
  library(enrichplot)
})

if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

# 3. 解析 GTF 搭建 ID 桥梁 (核心逻辑)
message("=> [1/5] 正在解析 GTF 提取 gene_id 与 protein_id 的对应关系...")
gtf <- import(opt$gtf_file)
gtf_df <- as.data.frame(mcols(gtf))

if(! "protein_id" %in% colnames(gtf_df)) {
  warning("GTF 文件中未找到 protein_id 属性，将尝试直接使用 gene_id 进行匹配...")
  id_map <- data.frame(gene_id = unique(gtf_df$gene_id), protein_id = unique(gtf_df$gene_id))
} else {
  id_map <- gtf_df %>%
    filter(!is.na(protein_id) & !is.na(gene_id)) %>%
    select(protein_id, gene_id) %>%
    distinct()
}

# 4. 读取数据并合并注释
message("=> [2/5] 正在读取并合并 eggNOG 注释与 DESeq2 结果...")

eggnog <- read_delim(opt$anno_file, delim="\t", comment="##", show_col_types=FALSE)
colnames(eggnog)[1] <- "protein_id"

eggnog_mapped <- eggnog %>% inner_join(id_map, by="protein_id")

# 读取 DESeq2 结果并合并
deg_data <- read_csv(opt$deg_file, show_col_types = FALSE)
final_table <- deg_data %>% left_join(eggnog_mapped, by=c("id" = "gene_id"))

write_csv(final_table, file.path(opt$out_dir, "Final_Master_Table_Annotated.csv"))

# 提取显著差异基因 (FDR < 0.05 且 |log2FC| > 1)
sig_genes <- final_table %>% 
  filter(FDR < 0.05 & abs(log2FC) > 1) %>% 
  pull(id) %>% 
  unique()

message(sprintf("=> 提取到 %d 个显著差异基因，准备进行富集分析...", length(sig_genes)))

# 5. 提取 clusterProfiler 需要的 TERM2GENE 背景
message("=> [3/5] 正在从 eggNOG 构建 GO 和 KEGG 的背景词典...")

# 构建 GO 背景
go_terms <- final_table %>% 
  select(id, GOs) %>% 
  filter(!is.na(GOs) & GOs != "-") %>%
  separate_rows(GOs, sep=",") %>% 
  select(GOs, id)

# 使用 clusterProfiler 内置函数获取 GO 名字
go_names <- go2term(unique(go_terms$GOs))
colnames(go_names) <- c("GOs", "Name")

# 构建 KEGG背景
kegg_terms <- final_table %>% 
  select(id, KEGG_Pathway) %>% 
  filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep=",") %>% 
  select(KEGG_Pathway, id) %>%
  # 清理带有映射前缀的 ID (如 map00010, ko00010 -> 00010)
  mutate(KEGG_Pathway = gsub("^[a-zA-Z]+", "", KEGG_Pathway))

# 6. 运行 GO 与 KEGG 气泡图富集
message("=> [4/5] 正在运行超几何分布富集分析...")

if(length(sig_genes) > 0) {
  # --- GO 富集 ---
  tryCatch({
    go_res <- enricher(sig_genes, TERM2GENE = go_terms, TERM2NAME = go_names, pvalueCutoff = 0.05)
    if(!is.null(go_res) && nrow(go_res@result %>% filter(p.adjust < 0.05)) > 0) {
      p_go <- dotplot(go_res, showCategory=20, title="GO Enrichment")
      ggsave(file.path(opt$out_dir, "GO_Enrichment_Dotplot.pdf"), p_go, width=8, height=6)
      write_csv(as.data.frame(go_res), file.path(opt$out_dir, "GO_Enrichment_Results.csv"))
    } else {
      message("  -> 提示: 没有显著富集的 GO Term。")
    }
  }, error = function(e) message("  -> GO 富集报错: ", e$message))
  
  # --- KEGG 富集 ---
  tryCatch({
    kegg_res <- enricher(sig_genes, TERM2GENE = kegg_terms, pvalueCutoff = 0.05)
    if(!is.null(kegg_res) && nrow(kegg_res@result %>% filter(p.adjust < 0.05)) > 0) {
      p_kegg <- dotplot(kegg_res, showCategory=20, title="KEGG Pathway Enrichment")
      ggsave(file.path(opt$out_dir, "KEGG_Enrichment_Dotplot.pdf"), p_kegg, width=8, height=6)
      write_csv(as.data.frame(kegg_res), file.path(opt$out_dir, "KEGG_Enrichment_Results.csv"))
    } else {
      message("  -> 提示: 没有显著富集的 KEGG Pathway。")
    }
  }, error = function(e) message("  -> KEGG 富集报错: ", e$message))
} else {
  message("  -> 提示: 没有差异基因，跳过气泡图富集分析。")
}

# 7. 运行 GSEA (Gene Set Enrichment Analysis)
message("=>[5/5] 正在运行全局基因的 GSEA (GO/KEGG) 分析...")

# 构造排序好的 named vector (GSEA 必须格式)
gene_list <- final_table$log2FC
names(gene_list) <- final_table$id
# 剔除 NA 并按 log2FC 降序排列
gene_list <- sort(na.omit(gene_list), decreasing = TRUE)

if(length(gene_list) > 100) {
  tryCatch({
    # 这里以 GO 为例跑全局 GSEA
    gsea_res <- GSEA(gene_list, TERM2GENE = go_terms, TERM2NAME = go_names, pvalueCutoff = 0.05, pAdjustMethod = "BH")
    
    if(!is.null(gsea_res) && nrow(gsea_res@result %>% filter(p.adjust < 0.05)) > 0) {
      # 获取最显著的一条通路画山峰图
      top_pathway <- gsea_res@result$ID[1]
      p_gsea <- gseaplot2(gsea_res, geneSetID = top_pathway, title = gsea_res@result$Description[1])
      ggsave(file.path(opt$out_dir, "GSEA_Plot.pdf"), p_gsea, width=8, height=6)
      write_csv(as.data.frame(gsea_res), file.path(opt$out_dir, "GSEA_GO_Results.csv"))
    } else {
      # 即使没有显著结果，创建一个空文件告诉 Snakemake 任务跑完了
      pdf(file.path(opt$out_dir, "GSEA_Plot.pdf"))
      plot(1, type="n", main="No significant GSEA result", xlab="", ylab="")
      dev.off()
      message("  -> 提示: GSEA 没有显著富集的通路。")
    }
  }, error = function(e) {
    message("  -> GSEA 分析报错: ", e$message)
    pdf(file.path(opt$out_dir, "GSEA_Plot.pdf")); plot(1, type="n", main="GSEA Error"); dev.off()
  })
}

message("富集分析全部完成！结果已保存在 ", opt$out_dir)
