# =============================================================================
# 软骨细胞亚群15差异基因分析
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")
load_packages()
library(pheatmap)
library(EnhancedVolcano)

print_progress("软骨细胞亚群15差异基因分析")

# 1. 加载数据 ----------------------------------------------------------------
message("加载软骨细胞数据...")
load(file.path(DATA_DIR, "chondrocyte_final.Rdata"))

# 确定使用的分辨率
if("RNA_snn_res.0.5" %in% colnames(chondrocyte_subset@meta.data)) {
  resolution <- "RNA_snn_res.0.5"
} else {
  resolution <- "RNA_snn_res.0.7"
}

Idents(chondrocyte_subset) <- resolution

# 查看聚类分布
message("聚类分布:")
print(table(chondrocyte_subset@meta.data[[resolution]]))

# 检查亚群15是否存在
available_clusters <- unique(chondrocyte_subset@meta.data[[resolution]])
target_cluster <- "15"

if(!target_cluster %in% available_clusters) {
  stop(paste0("错误: 亚群 '", target_cluster, "' 不存在！\n",
              "可用的亚群有: ", paste(available_clusters, collapse = ", ")))
}

message(paste("\n目标亚群:", target_cluster))
message(paste("亚群15细胞数:", sum(chondrocyte_subset@meta.data[[resolution]] == target_cluster)))

# 2. 差异基因分析 ------------------------------------------------------------
print_progress("步骤1: 差异基因分析")

message("分析亚群15 vs 其他所有亚群...")

# 使用Wilcoxon秩和检验
deg_wilcox <- FindMarkers(
  object = chondrocyte_subset,
  ident.1 = target_cluster,
  ident.2 = NULL,  # NULL表示与所有其他亚群比较
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = FALSE
)

deg_wilcox$gene <- rownames(deg_wilcox)
deg_wilcox$cluster <- target_cluster
deg_wilcox <- deg_wilcox %>% arrange(p_val_adj)

message(paste("发现", nrow(deg_wilcox), "个差异基因"))

# 保存完整结果
write.csv(
  deg_wilcox,
  file = file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_vs_others_all.csv")),
  row.names = FALSE
)

# 3. 筛选显著差异基因 --------------------------------------------------------
print_progress("步骤2: 筛选显著差异基因")

# 筛选显著差异基因
deg_sig <- deg_wilcox %>%
  filter(p_val_adj < 0.05) %>%
  mutate(
    change = case_when(
      avg_log2FC > 0.25 ~ "Up",
      avg_log2FC < -0.25 ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

n_up <- sum(deg_sig$change == "Up")
n_down <- sum(deg_sig$change == "Down")

message(paste("显著上调基因:", n_up))
message(paste("显著下调基因:", n_down))
message(paste("总显著差异基因:", n_up + n_down))

# 保存显著差异基因
write.csv(
  deg_sig,
  file = file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_vs_others_significant.csv")),
  row.names = FALSE
)

# 提取top基因
top_up_genes <- deg_sig %>%
  filter(change == "Up") %>%
  arrange(desc(avg_log2FC)) %>%
  head(50)

top_down_genes <- deg_sig %>%
  filter(change == "Down") %>%
  arrange(avg_log2FC) %>%
  head(50)

write.csv(
  top_up_genes,
  file = file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_top50_upregulated.csv")),
  row.names = FALSE
)

write.csv(
  top_down_genes,
  file = file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_top50_downregulated.csv")),
  row.names = FALSE
)

# 4. 火山图 ------------------------------------------------------------------
print_progress("步骤3: 绘制火山图")

# 准备火山图数据
volcano_data <- deg_wilcox %>%
  mutate(
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

# 基础火山图
p_volcano <- ggplot(volcano_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "Not Sig" = "grey80")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  labs(
    title = paste0("Chondrocyte Cluster ", target_cluster, " vs Others"),
    subtitle = paste0("Up: ", n_up, " genes | Down: ", n_down, " genes"),
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-value)",
    color = "Expression"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  )

save_plot(p_volcano, paste0("chondrocyte_cluster", target_cluster, "_volcano_plot.pdf"))

# EnhancedVolcano（如果包可用）
if(requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  # 标注top基因
  top_genes_to_label <- c(
    head(top_up_genes$gene, 10),
    head(top_down_genes$gene, 10)
  )
  
  p_volcano_enhanced <- EnhancedVolcano(
    volcano_data,
    lab = volcano_data$gene,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = paste0('Chondrocyte Cluster ', target_cluster, ' vs Others'),
    subtitle = paste0("Up: ", n_up, " | Down: ", n_down),
    pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 2.0,
    labSize = 4.0,
    selectLab = top_genes_to_label,
    col = c('grey80', 'grey80', '#4DBBD5', '#E64B35'),
    colAlpha = 0.5,
    legendPosition = 'right',
    drawConnectors = TRUE,
    widthConnectors = 0.5
  )
  
  save_plot(p_volcano_enhanced, 
            paste0("chondrocyte_cluster", target_cluster, "_volcano_enhanced.pdf"),
            width = 12, height = 10)
}

# 5. 热图 --------------------------------------------------------------------
print_progress("步骤4: 绘制热图")

# 提取top50差异基因
top_genes_heatmap <- c(
  head(top_up_genes$gene, 25),
  head(top_down_genes$gene, 25)
)

if(length(top_genes_heatmap) > 0) {
  # 提取表达矩阵
  expr_matrix <- GetAssayData(chondrocyte_subset, slot = "data")[top_genes_heatmap, ]
  
  # 标准化
  expr_scaled <- t(scale(t(as.matrix(expr_matrix))))
  
  # 准备列注释
  col_anno <- data.frame(
    Cluster = chondrocyte_subset@meta.data[[resolution]],
    Group = chondrocyte_subset$group,
    row.names = colnames(chondrocyte_subset)
  )
  col_anno$IsCluster15 <- ifelse(col_anno$Cluster == target_cluster, "Cluster15", "Others")
  
  # 绘制热图
  pdf(file.path(FIGURE_DIR, paste0("chondrocyte_cluster", target_cluster, "_heatmap_top50.pdf")),
      width = 12, height = 10)
  
  pheatmap(
    expr_scaled,
    color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = col_anno[, c("IsCluster15", "Group"), drop = FALSE],
    annotation_colors = list(
      IsCluster15 = c("Cluster15" = "#E64B35", "Others" = "grey70"),
      Group = c("MDD" = "#00A087", "SDD" = "#DC0000")
    ),
    fontsize_row = 6,
    main = paste0("Top 50 DEGs in Cluster ", target_cluster),
    border_color = NA
  )
  
  dev.off()
}

# 6. 小提琴图 ----------------------------------------------------------------
print_progress("步骤5: 绘制小提琴图")

# 选择top12个基因
top_genes_violin <- c(
  head(top_up_genes$gene, 6),
  head(top_down_genes$gene, 6)
)

if(length(top_genes_violin) > 0) {
  # 设置ident为目标cluster vs others
  chondrocyte_subset$cluster_group <- ifelse(
    chondrocyte_subset@meta.data[[resolution]] == target_cluster,
    paste0("Cluster_", target_cluster),
    "Others"
  )
  
  Idents(chondrocyte_subset) <- "cluster_group"
  
  p_violin <- VlnPlot(
    chondrocyte_subset,
    features = top_genes_violin,
    pt.size = 0,
    ncol = 3
  )
  
  for(i in 1:length(p_violin)) {
    p_violin[[i]] <- p_violin[[i]] + scale_fill_custom()
  }
  
  p_violin <- p_violin &
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 10)
    )
  
  save_plot(p_violin, 
            paste0("chondrocyte_cluster", target_cluster, "_violin_top12.pdf"),
            width = 15, height = 12)
}

# 7. 特征图 ------------------------------------------------------------------
print_progress("步骤6: 绘制特征图")

# 选择top9上调基因
top9_up <- head(top_up_genes$gene, 9)

if(length(top9_up) > 0) {
  # 恢复原始ident
  Idents(chondrocyte_subset) <- resolution
  
  p_feature <- FeaturePlot(
    chondrocyte_subset,
    features = top9_up,
    reduction = "umap",
    ncol = 3,
    pt.size = 0.5
  ) &
    scale_color_gradientn(colors = c("grey90", "yellow", "orange", "red", "darkred"))
  
  save_plot(p_feature, 
            paste0("chondrocyte_cluster", target_cluster, "_feature_top9.pdf"),
            width = 15, height = 15)
}

# 8. DotPlot -----------------------------------------------------------------
print_progress("步骤7: 绘制DotPlot")

# 选择top20基因
top20_genes <- c(
  head(top_up_genes$gene, 10),
  head(top_down_genes$gene, 10)
)

if(length(top20_genes) > 0) {
  Idents(chondrocyte_subset) <- resolution
  
  p_dot <- DotPlot(
    chondrocyte_subset,
    features = top20_genes,
    cols = c("grey90", "#E64B35"),
    dot.scale = 8
  ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Top 20 DEGs - Cluster ", target_cluster))
  
  save_plot(p_dot, 
            paste0("chondrocyte_cluster", target_cluster, "_dotplot_top20.pdf"),
            width = 12, height = 8)
}

# 9. 按分组比较（可选）------------------------------------------------------
print_progress("步骤8: 亚群15在SDD vs MDD中的差异")

# 提取亚群15的细胞
cluster15_cells <- WhichCells(chondrocyte_subset, 
                               expression = chondrocyte_subset@meta.data[[resolution]] == target_cluster)
cluster15_subset <- subset(chondrocyte_subset, cells = cluster15_cells)

message("亚群15细胞分组分布:")
print(table(cluster15_subset$group))

# 如果两组都有足够的细胞，进行组间比较
group_table <- table(cluster15_subset$group)
if(all(group_table >= 10)) {
  message("\n在亚群15内进行SDD vs MDD差异分析...")
  
  Idents(cluster15_subset) <- "group"
  
  deg_cluster15_group <- FindMarkers(
    cluster15_subset,
    ident.1 = "SDD",
    ident.2 = "MDD",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = FALSE
  )
  
  deg_cluster15_group$gene <- rownames(deg_cluster15_group)
  deg_cluster15_group <- deg_cluster15_group %>% arrange(p_val_adj)
  
  # 保存结果
  write.csv(
    deg_cluster15_group,
    file = file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_SDD_vs_MDD.csv")),
    row.names = FALSE
  )
  
  message(paste("在亚群15中发现", sum(deg_cluster15_group$p_val_adj < 0.05), "个SDD vs MDD差异基因"))
  
  # 火山图
  p_volcano_group <- plot_volcano(
    deg_cluster15_group,
    title = paste0("Cluster ", target_cluster, ": SDD vs MDD")
  )
  save_plot(p_volcano_group, 
            paste0("chondrocyte_cluster", target_cluster, "_SDD_vs_MDD_volcano.pdf"))
}

# 10. 统计摘要 ---------------------------------------------------------------
print_progress("步骤9: 生成统计摘要")

summary_report <- data.frame(
  Metric = c(
    "Target Cluster",
    "Cluster15 Cell Count",
    "Total Genes Tested",
    "Significant DEGs (FDR < 0.05)",
    "Upregulated Genes",
    "Downregulated Genes",
    "Max log2FC (Up)",
    "Min log2FC (Down)",
    "Most Significant Gene (Up)",
    "Most Significant Gene (Down)"
  ),
  Value = c(
    target_cluster,
    sum(chondrocyte_subset@meta.data[[resolution]] == target_cluster),
    nrow(deg_wilcox),
    nrow(deg_sig),
    n_up,
    n_down,
    round(max(deg_sig$avg_log2FC[deg_sig$change == "Up"], na.rm = TRUE), 3),
    round(min(deg_sig$avg_log2FC[deg_sig$change == "Down"], na.rm = TRUE), 3),
    ifelse(n_up > 0, deg_sig$gene[deg_sig$change == "Up"][1], "NA"),
    ifelse(n_down > 0, deg_sig$gene[deg_sig$change == "Down"][1], "NA")
  )
)

write.csv(
  summary_report,
  file = file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_summary.csv")),
  row.names = FALSE
)

print(summary_report)

# 11. 保存工作空间 -----------------------------------------------------------
save.image(file = file.path(DATA_DIR, paste0("cluster", target_cluster, "_DEG_analysis.RData")))

print_progress("分析完成")

message("\n结果文件:")
message(paste("- 完整差异基因:", file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_vs_others_all.csv"))))
message(paste("- 显著差异基因:", file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_vs_others_significant.csv"))))
message(paste("- Top50上调基因:", file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_top50_upregulated.csv"))))
message(paste("- Top50下调基因:", file.path(RESULT_DIR, paste0("chondrocyte_cluster", target_cluster, "_top50_downregulated.csv"))))
message(paste("- 图片保存于:", FIGURE_DIR))

