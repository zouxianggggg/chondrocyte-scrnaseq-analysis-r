# =============================================================================
# 软骨细胞和髓核细胞共同保护通路分析及HDAC9验证
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")
load_packages()
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(VennDiagram)

print_progress("软骨细胞和髓核细胞共同保护通路分析及HDAC9验证")

# 参数设置 -------------------------------------------------------------------
TARGET_GENE <- "HDAC9"
PVALUE_CUTOFF <- 0.05
QVALUE_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0.25

# =============================================================================
# 一、加载数据并提取差异基因
# =============================================================================
print_progress("步骤1: 加载数据并识别差异基因")

# 1.1 加载软骨细胞数据
message("\n加载软骨细胞数据...")
load(file.path(DATA_DIR, "chondrocyte_clustered.Rdata"))
chondrocyte_data <- chondrocyte_subset
rm(chondrocyte_subset)

message(paste("软骨细胞数:", ncol(chondrocyte_data)))
message(paste("基因数:", nrow(chondrocyte_data)))
message("分组统计:")
print(table(chondrocyte_data$group))

# 1.2 加载髓核细胞数据
message("\n加载髓核细胞数据...")
load(file.path(DATA_DIR, "np_clustered.Rdata"))
np_data <- np_subset
rm(np_subset)

message(paste("髓核细胞数:", ncol(np_data)))
message(paste("基因数:", nrow(np_data)))
message("分组统计:")
print(table(np_data$group))

# 1.3 软骨细胞差异分析
message("\n进行软骨细胞差异分析 (SDD vs MDD)...")
Idents(chondrocyte_data) <- chondrocyte_data$group

chondrocyte_deg <- FindMarkers(
  chondrocyte_data,
  ident.1 = "SDD",
  ident.2 = "MDD",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1,
  verbose = TRUE
)

chondrocyte_deg$gene <- rownames(chondrocyte_deg)
chondrocyte_deg$cell_type <- "Chondrocyte"

# 筛选显著差异基因
chondrocyte_deg_sig <- chondrocyte_deg %>%
  filter(p_val_adj < PVALUE_CUTOFF, abs(avg_log2FC) > LOG2FC_CUTOFF)

message(paste("软骨细胞显著差异基因数:", nrow(chondrocyte_deg_sig)))
message(paste("  上调基因:", sum(chondrocyte_deg_sig$avg_log2FC > 0)))
message(paste("  下调基因:", sum(chondrocyte_deg_sig$avg_log2FC < 0)))

# 1.4 髓核细胞差异分析
message("\n进行髓核细胞差异分析 (SDD vs MDD)...")
Idents(np_data) <- np_data$group

np_deg <- FindMarkers(
  np_data,
  ident.1 = "SDD",
  ident.2 = "MDD",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1,
  verbose = TRUE
)

np_deg$gene <- rownames(np_deg)
np_deg$cell_type <- "Nucleus_Pulposus"

# 筛选显著差异基因
np_deg_sig <- np_deg %>%
  filter(p_val_adj < PVALUE_CUTOFF, abs(avg_log2FC) > LOG2FC_CUTOFF)

message(paste("髓核细胞显著差异基因数:", nrow(np_deg_sig)))
message(paste("  上调基因:", sum(np_deg_sig$avg_log2FC > 0)))
message(paste("  下调基因:", sum(np_deg_sig$avg_log2FC < 0)))

# 1.5 保存差异基因结果
write.csv(chondrocyte_deg, 
          file = file.path(RESULT_DIR, "10_软骨细胞_差异基因_全部.csv"),
          row.names = FALSE)
write.csv(chondrocyte_deg_sig, 
          file = file.path(RESULT_DIR, "10_软骨细胞_差异基因_显著.csv"),
          row.names = FALSE)
write.csv(np_deg, 
          file = file.path(RESULT_DIR, "10_髓核细胞_差异基因_全部.csv"),
          row.names = FALSE)
write.csv(np_deg_sig, 
          file = file.path(RESULT_DIR, "10_髓核细胞_差异基因_显著.csv"),
          row.names = FALSE)

# 1.6 检查HDAC9
message("\n检查HDAC9在两种细胞中的差异表达:")
if(TARGET_GENE %in% chondrocyte_deg$gene) {
  hdac9_chon <- chondrocyte_deg[chondrocyte_deg$gene == TARGET_GENE, ]
  message(paste("软骨细胞 - HDAC9 log2FC:", round(hdac9_chon$avg_log2FC, 3),
                ", p_adj:", format(hdac9_chon$p_val_adj, scientific = TRUE)))
} else {
  message("软骨细胞 - HDAC9: 未检测到")
  hdac9_chon <- NULL
}

if(TARGET_GENE %in% np_deg$gene) {
  hdac9_np <- np_deg[np_deg$gene == TARGET_GENE, ]
  message(paste("髓核细胞 - HDAC9 log2FC:", round(hdac9_np$avg_log2FC, 3),
                ", p_adj:", format(hdac9_np$p_val_adj, scientific = TRUE)))
} else {
  message("髓核细胞 - HDAC9: 未检测到")
  hdac9_np <- NULL
}

# =============================================================================
# 二、识别共同差异基因（下调基因，保护性基因）
# =============================================================================
print_progress("步骤2: 识别共同下调的保护性基因")

# 2.1 提取下调基因（在病变中表达降低 = 保护性基因）
chondrocyte_down <- chondrocyte_deg_sig %>%
  filter(avg_log2FC < 0) %>%
  pull(gene)

np_down <- np_deg_sig %>%
  filter(avg_log2FC < 0) %>%
  pull(gene)

message(paste("\n软骨细胞下调基因数:", length(chondrocyte_down)))
message(paste("髓核细胞下调基因数:", length(np_down)))

# 2.2 共同下调基因
common_down_genes <- intersect(chondrocyte_down, np_down)
message(paste("共同下调基因数:", length(common_down_genes)))

if(length(common_down_genes) > 0) {
  message("\n共同下调基因列表:")
  print(head(common_down_genes, 20))
  
  # 检查HDAC9是否在其中
  if(TARGET_GENE %in% common_down_genes) {
    message(paste("\n*** ", TARGET_GENE, " 在共同下调基因中 ***"))
  } else {
    message(paste("\n警告: ", TARGET_GENE, " 不在共同下调基因中"))
  }
}

# 2.3 保存共同基因
write.csv(
  data.frame(Gene = common_down_genes),
  file = file.path(RESULT_DIR, "10_共同下调基因.csv"),
  row.names = FALSE
)

# 2.4 Venn图展示
venn_list <- list(
  Chondrocyte = chondrocyte_down,
  NucleusPulposus = np_down
)

pdf(file.path(FIGURE_DIR, "10_共同下调基因_Venn图.pdf"), width = 8, height = 8)
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("Chondrocyte\nDownregulated", "Nucleus Pulposus\nDownregulated"),
  filename = NULL,
  fill = c("#00A087", "#DC0000"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.1
)
grid.draw(venn.plot)
dev.off()

# =============================================================================
# 三、共同下调基因的通路富集分析（保护性通路）
# =============================================================================
print_progress("步骤3: 共同保护性通路富集分析")

if(length(common_down_genes) < 10) {
  message("警告: 共同下调基因数量太少(<10)，富集分析可能不可靠")
  message("将使用更宽松的标准重新筛选...")
  
  # 使用更宽松的标准
  chondrocyte_down_relaxed <- chondrocyte_deg %>%
    filter(p_val_adj < 0.1, avg_log2FC < -0.1) %>%
    pull(gene)
  
  np_down_relaxed <- np_deg %>%
    filter(p_val_adj < 0.1, avg_log2FC < -0.1) %>%
    pull(gene)
  
  common_down_genes <- intersect(chondrocyte_down_relaxed, np_down_relaxed)
  message(paste("放宽标准后共同下调基因数:", length(common_down_genes)))
}

# 3.1 基因ID转换
message("\n转换基因ID...")
gene_id_conversion <- bitr(
  common_down_genes,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "ENSEMBL"),
  OrgDb = org.Hs.eg.db
)

message(paste("成功转换:", nrow(gene_id_conversion), "/", length(common_down_genes)))
common_entrez <- gene_id_conversion$ENTREZID

# 3.2 GO富集分析
message("\n进行GO富集分析...")
go_result <- enrichGO(
  gene = common_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = PVALUE_CUTOFF,
  qvalueCutoff = QVALUE_CUTOFF,
  readable = TRUE
)

if(!is.null(go_result) && nrow(go_result@result) > 0) {
  message(paste("GO BP 富集到", nrow(go_result@result), "个通路"))
  
  # 保存结果
  write.csv(
    go_result@result,
    file = file.path(RESULT_DIR, "10_共同下调基因_GO富集.csv"),
    row.names = FALSE
  )
  
  # 可视化
  if(nrow(go_result@result) >= 10) {
    p_go_bar <- barplot(go_result, showCategory = 20) +
      labs(title = "")
    save_plot(p_go_bar, "10_共同保护通路_GO_barplot.pdf", width = 12, height = 10)
    
    p_go_dot <- dotplot(go_result, showCategory = 20) +
      labs(title = "")
    save_plot(p_go_dot, "10_共同保护通路_GO_dotplot.pdf", width = 12, height = 10)
  }
} else {
  message("GO富集分析未找到显著通路")
}

# 3.3 KEGG富集分析
message("\n进行KEGG富集分析...")
kegg_result <- enrichKEGG(
  gene = common_entrez,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = PVALUE_CUTOFF,
  qvalueCutoff = QVALUE_CUTOFF
)

if(!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  # 转换为可读基因名
  kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  message(paste("KEGG 富集到", nrow(kegg_result@result), "个通路"))
  
  # 保存结果
  write.csv(
    kegg_result@result,
    file = file.path(RESULT_DIR, "10_共同下调基因_KEGG富集.csv"),
    row.names = FALSE
  )
  
  # 可视化
  if(nrow(kegg_result@result) >= 5) {
    p_kegg_bar <- barplot(kegg_result, showCategory = 20) +
      labs(title = "")
    save_plot(p_kegg_bar, "10_共同保护通路_KEGG_barplot.pdf", width = 12, height = 10)
    
    p_kegg_dot <- dotplot(kegg_result, showCategory = 20) +
      labs(title = "")
    save_plot(p_kegg_dot, "10_共同保护通路_KEGG_dotplot.pdf", width = 12, height = 10)
  }
} else {
  message("KEGG富集分析未找到显著通路")
}

# =============================================================================
# 四、HDAC9在富集通路中的验证
# =============================================================================
print_progress("步骤4: 验证HDAC9在富集通路中的作用")

# 4.1 检查HDAC9是否在富集通路中
hdac9_pathways <- list()

# 检查GO通路
if(!is.null(go_result) && nrow(go_result@result) > 0) {
  message("\n检查HDAC9在GO通路中...")
  
  go_with_hdac9 <- go_result@result %>%
    filter(grepl(TARGET_GENE, geneID)) %>%
    select(ID, Description, GeneRatio, pvalue, p.adjust, geneID)
  
  if(nrow(go_with_hdac9) > 0) {
    message(paste("*** HDAC9 出现在", nrow(go_with_hdac9), "个GO通路中 ***"))
    hdac9_pathways$GO <- go_with_hdac9
    
    write.csv(
      go_with_hdac9,
      file = file.path(RESULT_DIR, "10_HDAC9相关GO通路.csv"),
      row.names = FALSE
    )
    
    print(go_with_hdac9[, c("Description", "pvalue", "p.adjust")])
  } else {
    message("HDAC9 未出现在显著的GO通路中")
  }
}

# 检查KEGG通路
if(!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  message("\n检查HDAC9在KEGG通路中...")
  
  kegg_with_hdac9 <- kegg_result@result %>%
    filter(grepl(TARGET_GENE, geneID)) %>%
    select(ID, Description, GeneRatio, pvalue, p.adjust, geneID)
  
  if(nrow(kegg_with_hdac9) > 0) {
    message(paste("*** HDAC9 出现在", nrow(kegg_with_hdac9), "个KEGG通路中 ***"))
    hdac9_pathways$KEGG <- kegg_with_hdac9
    
    write.csv(
      kegg_with_hdac9,
      file = file.path(RESULT_DIR, "10_HDAC9相关KEGG通路.csv"),
      row.names = FALSE
    )
    
    print(kegg_with_hdac9[, c("Description", "pvalue", "p.adjust")])
  } else {
    message("HDAC9 未出现在显著的KEGG通路中")
  }
}

# 4.2 即使HDAC9不在显著富集的通路中，也可以分析与它功能相关的通路
message("\n识别与HDAC9功能相关的通路（表观遗传、转录调控等）...")

# 定义HDAC9相关的功能关键词
hdac9_keywords <- c(
  "histone", "acetylation", "deacetylase", "chromatin", "transcription",
  "gene expression", "epigenetic", "regulation of transcription",
  "DNA-templated transcription", "RNA polymerase"
)

# 在GO结果中搜索
if(!is.null(go_result) && nrow(go_result@result) > 0) {
  go_hdac9_related <- go_result@result %>%
    filter(grepl(paste(hdac9_keywords, collapse = "|"), Description, ignore.case = TRUE))
  
  if(nrow(go_hdac9_related) > 0) {
    message(paste("找到", nrow(go_hdac9_related), "个与HDAC9功能相关的GO通路"))
    
    write.csv(
      go_hdac9_related,
      file = file.path(RESULT_DIR, "10_HDAC9功能相关GO通路.csv"),
      row.names = FALSE
    )
    
    # 可视化
    if(nrow(go_hdac9_related) >= 5) {
      go_subset <- go_result
      go_subset@result <- go_hdac9_related[1:min(20, nrow(go_hdac9_related)), ]
      
      p_go_hdac9 <- dotplot(go_subset, showCategory = 20) +
        labs(title = "")
      save_plot(p_go_hdac9, "10_HDAC9功能相关通路_dotplot.pdf", width = 12, height = 10)
    }
  }
}

# 在KEGG结果中搜索
if(!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  kegg_hdac9_related <- kegg_result@result %>%
    filter(grepl(paste(hdac9_keywords, collapse = "|"), Description, ignore.case = TRUE))
  
  if(nrow(kegg_hdac9_related) > 0) {
    message(paste("找到", nrow(kegg_hdac9_related), "个与HDAC9功能相关的KEGG通路"))
    
    write.csv(
      kegg_hdac9_related,
      file = file.path(RESULT_DIR, "10_HDAC9功能相关KEGG通路.csv"),
      row.names = FALSE
    )
  }
}

# =============================================================================
# 五、HDAC9表达验证
# =============================================================================
print_progress("步骤5: HDAC9表达模式验证")

# 5.1 提取HDAC9表达数据
if(TARGET_GENE %in% rownames(chondrocyte_data)) {
  chondrocyte_hdac9 <- FetchData(
    chondrocyte_data,
    vars = c(TARGET_GENE, "group")
  )
  chondrocyte_hdac9$cell_type <- "Chondrocyte"
} else {
  chondrocyte_hdac9 <- NULL
  message("警告: HDAC9不在软骨细胞数据中")
}

if(TARGET_GENE %in% rownames(np_data)) {
  np_hdac9 <- FetchData(
    np_data,
    vars = c(TARGET_GENE, "group")
  )
  np_hdac9$cell_type <- "Nucleus_Pulposus"
} else {
  np_hdac9 <- NULL
  message("警告: HDAC9不在髓核细胞数据中")
}

# 5.2 合并数据
if(!is.null(chondrocyte_hdac9) && !is.null(np_hdac9)) {
  hdac9_expr <- rbind(chondrocyte_hdac9, np_hdac9)
  colnames(hdac9_expr)[1] <- "HDAC9"
  
  # 统计
  hdac9_stats <- hdac9_expr %>%
    group_by(cell_type, group) %>%
    summarise(
      n = n(),
      mean = mean(HDAC9),
      median = median(HDAC9),
      sd = sd(HDAC9),
      .groups = "drop"
    )
  
  message("\nHDAC9表达统计:")
  print(hdac9_stats)
  
  write.csv(
    hdac9_stats,
    file = file.path(RESULT_DIR, "10_HDAC9表达统计.csv"),
    row.names = FALSE
  )
  
  # 5.3 可视化HDAC9表达
  # 小提琴图
  p_hdac9_violin <- ggplot(hdac9_expr, aes(x = cell_type, y = HDAC9, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), 
                 fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000")) +
    theme_classic() +
    labs(
      x = "Cell Type",
      y = "HDAC9 Expression",
      fill = "Group"
    ) +
    stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test") +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot(p_hdac9_violin, "10_HDAC9表达_violin图.pdf", width = 10, height = 8)
  
  # 分面小提琴图
  p_hdac9_facet <- ggplot(hdac9_expr, aes(x = group, y = HDAC9, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    facet_wrap(~cell_type, scales = "free_y") +
    scale_fill_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000")) +
    theme_bw() +
    labs(
      x = "Group",
      y = "HDAC9 Expression"
    ) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
    theme(
      text = element_text(size = 14),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  save_plot(p_hdac9_facet, "10_HDAC9表达_分面图.pdf", width = 12, height = 6)
  
  # 柱状图 + 误差棒
  p_hdac9_bar <- ggplot(hdac9_stats, aes(x = cell_type, y = mean, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  position = position_dodge(0.9), width = 0.2) +
    scale_fill_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000")) +
    theme_classic() +
    labs(
      x = "Cell Type",
      y = "Mean HDAC9 Expression",
      fill = "Group"
    ) +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot(p_hdac9_bar, "10_HDAC9表达_柱状图.pdf", width = 10, height = 8)
}

# 5.4 UMAP可视化HDAC9表达
if(TARGET_GENE %in% rownames(chondrocyte_data)) {
  p_chon_umap <- FeaturePlot(
    chondrocyte_data,
    features = TARGET_GENE,
    split.by = "group",
    cols = c("lightgrey", "red")
  ) +
    plot_annotation(title = "")
  
  save_plot(p_chon_umap, "10_HDAC9表达_软骨细胞UMAP.pdf", width = 14, height = 6)
}

if(TARGET_GENE %in% rownames(np_data)) {
  p_np_umap <- FeaturePlot(
    np_data,
    features = TARGET_GENE,
    split.by = "group",
    cols = c("lightgrey", "red")
  ) +
    plot_annotation(title = "")
  
  save_plot(p_np_umap, "10_HDAC9表达_髓核细胞UMAP.pdf", width = 14, height = 6)
}

# =============================================================================
# 六、关键保护性通路的深入分析
# =============================================================================
print_progress("步骤6: 关键保护性通路深入分析")

# 6.1 识别top保护性通路
message("\n识别最显著的保护性通路...")

top_pathways <- list()

if(!is.null(go_result) && nrow(go_result@result) > 0) {
  top_go <- go_result@result %>%
    arrange(p.adjust) %>%
    head(10)
  
  top_pathways$GO <- top_go
  
  message("\nTop 10 GO通路:")
  print(top_go[, c("Description", "Count", "pvalue", "p.adjust")])
}

if(!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  top_kegg <- kegg_result@result %>%
    arrange(p.adjust) %>%
    head(10)
  
  top_pathways$KEGG <- top_kegg
  
  message("\nTop 10 KEGG通路:")
  print(top_kegg[, c("Description", "Count", "pvalue", "p.adjust")])
}

# 6.2 通路网络图
if(!is.null(go_result) && nrow(go_result@result) >= 5) {
  message("\n生成GO通路网络图...")
  
  # 使用前30个通路
  go_subset <- go_result
  go_subset@result <- go_result@result[1:min(30, nrow(go_result@result)), ]
  
  pdf(file.path(FIGURE_DIR, "10_保护通路_GO网络图.pdf"), width = 14, height = 12)
  print(emapplot(pairwise_termsim(go_subset), showCategory = 30, cex_label_category = 0.6))
  dev.off()
  
  # Cnetplot - 显示基因和通路的关系
  pdf(file.path(FIGURE_DIR, "10_保护通路_GO基因网络.pdf"), width = 16, height = 14)
  print(cnetplot(go_subset, showCategory = 10, 
                 cex_label_gene = 0.7, cex_label_category = 1))
  dev.off()
}

if(!is.null(kegg_result) && nrow(kegg_result@result) >= 5) {
  message("\n生成KEGG通路网络图...")
  
  kegg_subset <- kegg_result
  kegg_subset@result <- kegg_result@result[1:min(20, nrow(kegg_result@result)), ]
  
  pdf(file.path(FIGURE_DIR, "10_保护通路_KEGG网络图.pdf"), width = 14, height = 12)
  print(emapplot(pairwise_termsim(kegg_subset), showCategory = 20, cex_label_category = 0.7))
  dev.off()
  
  pdf(file.path(FIGURE_DIR, "10_保护通路_KEGG基因网络.pdf"), width = 16, height = 14)
  print(cnetplot(kegg_subset, showCategory = 10,
                 cex_label_gene = 0.7, cex_label_category = 1))
  dev.off()
}

# =============================================================================
# 七、生成综合报告
# =============================================================================
print_progress("生成综合分析报告")

summary_text <- c(
  "=============================================================================",
  "软骨细胞和髓核细胞共同保护通路分析及HDAC9验证报告",
  "=============================================================================",
  "",
  paste("分析日期:", Sys.time()),
  "",
  "--- 数据概况 ---",
  paste("软骨细胞数:", ncol(chondrocyte_data)),
  paste("  MDD:", sum(chondrocyte_data$group == "MDD")),
  paste("  SDD:", sum(chondrocyte_data$group == "SDD")),
  paste("髓核细胞数:", ncol(np_data)),
  paste("  MDD:", sum(np_data$group == "MDD")),
  paste("  SDD:", sum(np_data$group == "SDD")),
  "",
  "--- 差异基因分析 ---",
  paste("软骨细胞显著差异基因:", nrow(chondrocyte_deg_sig)),
  paste("  上调:", sum(chondrocyte_deg_sig$avg_log2FC > 0)),
  paste("  下调:", sum(chondrocyte_deg_sig$avg_log2FC < 0)),
  paste("髓核细胞显著差异基因:", nrow(np_deg_sig)),
  paste("  上调:", sum(np_deg_sig$avg_log2FC > 0)),
  paste("  下调:", sum(np_deg_sig$avg_log2FC < 0)),
  "",
  "--- 共同保护性基因 ---",
  paste("共同下调基因数:", length(common_down_genes)),
  ifelse(TARGET_GENE %in% common_down_genes,
         paste("*** HDAC9在共同下调基因中 ***"),
         paste("HDAC9不在共同下调基因中")),
  ""
)

# HDAC9差异表达
if(!is.null(hdac9_chon)) {
  summary_text <- c(
    summary_text,
    "--- HDAC9差异表达 ---",
    paste("软骨细胞:"),
    paste("  log2FC:", round(hdac9_chon$avg_log2FC, 3)),
    paste("  p_adj:", format(hdac9_chon$p_val_adj, scientific = TRUE)),
    paste("  趋势:", ifelse(hdac9_chon$avg_log2FC < 0, "病变中降低 (保护性)", "病变中升高"))
  )
}

if(!is.null(hdac9_np)) {
  summary_text <- c(
    summary_text,
    paste("髓核细胞:"),
    paste("  log2FC:", round(hdac9_np$avg_log2FC, 3)),
    paste("  p_adj:", format(hdac9_np$p_val_adj, scientific = TRUE)),
    paste("  趋势:", ifelse(hdac9_np$avg_log2FC < 0, "病变中降低 (保护性)", "病变中升高")),
    ""
  )
}

# 富集分析结果
if(!is.null(go_result) && nrow(go_result@result) > 0) {
  summary_text <- c(
    summary_text,
    "--- GO富集分析 ---",
    paste("富集到显著GO通路:", nrow(go_result@result)),
    paste("Top 5 GO通路:"),
    paste0("  ", 1:min(5, nrow(go_result@result)), ". ", 
           go_result@result$Description[1:min(5, nrow(go_result@result))]),
    ""
  )
}

if(!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
  summary_text <- c(
    summary_text,
    "--- KEGG富集分析 ---",
    paste("富集到显著KEGG通路:", nrow(kegg_result@result)),
    paste("Top 5 KEGG通路:"),
    paste0("  ", 1:min(5, nrow(kegg_result@result)), ". ", 
           kegg_result@result$Description[1:min(5, nrow(kegg_result@result))]),
    ""
  )
}

# HDAC9在通路中的验证
if(length(hdac9_pathways) > 0) {
  summary_text <- c(
    summary_text,
    "--- HDAC9通路验证 ---"
  )
  
  if(!is.null(hdac9_pathways$GO)) {
    summary_text <- c(
      summary_text,
      paste("*** HDAC9出现在", nrow(hdac9_pathways$GO), "个显著GO通路中 ***"),
      paste("相关通路:"),
      paste0("  - ", hdac9_pathways$GO$Description)
    )
  }
  
  if(!is.null(hdac9_pathways$KEGG)) {
    summary_text <- c(
      summary_text,
      paste("*** HDAC9出现在", nrow(hdac9_pathways$KEGG), "个显著KEGG通路中 ***"),
      paste("相关通路:"),
      paste0("  - ", hdac9_pathways$KEGG$Description)
    )
  }
  summary_text <- c(summary_text, "")
}

summary_text <- c(
  summary_text,
  "=============================================================================",
  "结论:",
  "1. 识别出两种细胞类型的共同下调基因（保护性基因）",
  ifelse(TARGET_GENE %in% common_down_genes,
         "2. HDAC9在两种细胞类型中均下调，符合保护性基因特征",
         "2. HDAC9表达变化需进一步验证"),
  "3. 共同下调基因富集到的通路代表共同保护机制",
  ifelse(length(hdac9_pathways) > 0,
         "4. *** HDAC9在显著富集通路中得到验证 ***",
         "4. HDAC9在富集通路中的作用需要进一步实验验证"),
  "5. 详细结果请查看生成的图表和数据文件",
  "============================================================================="
)

writeLines(summary_text, file.path(RESULT_DIR, "10_共同保护通路_HDAC9验证_总结.txt"))

print_progress("分析完成！")

message("\n生成的结果文件:")
message(paste("- 数据表格:", RESULT_DIR))
message(paste("- 可视化图表:", FIGURE_DIR))
message(paste("- 综合报告:", file.path(RESULT_DIR, "10_共同保护通路_HDAC9验证_总结.txt")))

message("\n关键发现:")
if(TARGET_GENE %in% common_down_genes) {
  message(paste("✓", TARGET_GENE, "在两种细胞类型中均下调"))
}
if(length(hdac9_pathways) > 0) {
  message(paste("✓", TARGET_GENE, "在显著富集通路中得到验证"))
}
message("✓ 识别出共同保护性通路")
