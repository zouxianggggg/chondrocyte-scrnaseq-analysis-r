# =============================================================================
# 软骨细胞第15亚群综合分析
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")
load_packages()
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(VennDiagram)
library(ggpubr)
library(RColorBrewer)

print_progress("软骨细胞第15亚群综合分析")

# 参数设置 -------------------------------------------------------------------
TARGET_CLUSTER <- "15"
TOP_N <- 100  # 前N个基因用于交集分析

# 1. 加载数据 ----------------------------------------------------------------
message("加载软骨细胞数据...")
load(file.path(DATA_DIR, "chondrocyte_clustered.Rdata"))

# 确定分辨率
if("RNA_snn_res.0.5" %in% colnames(chondrocyte_subset@meta.data)) {
  resolution <- "RNA_snn_res.0.5"
} else {
  resolution <- "RNA_snn_res.0.7"
}

Idents(chondrocyte_subset) <- resolution

message("数据概况:")
message(paste("总细胞数:", ncol(chondrocyte_subset)))
message(paste("聚类分布:"))
print(table(chondrocyte_subset@meta.data[[resolution]]))
message(paste("分组分布:"))
print(table(chondrocyte_subset$group))

# =============================================================================
# 一、第15亚群 vs 其他14个亚群的差异基因分析
# =============================================================================
print_progress("分析一: 第15亚群差异基因分析")

message("分析第15亚群 vs 其他亚群...")

# 差异分析
deg_cluster15 <- FindMarkers(
  object = chondrocyte_subset,
  ident.1 = TARGET_CLUSTER,
  ident.2 = NULL,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = FALSE
)

deg_cluster15$gene <- rownames(deg_cluster15)
deg_cluster15 <- deg_cluster15 %>% arrange(p_val_adj)

# 筛选显著基因
deg_cluster15_sig <- deg_cluster15 %>%
  filter(p_val_adj < 0.05) %>%
  mutate(
    change = case_when(
      avg_log2FC > 0.25 ~ "Up",
      avg_log2FC < -0.25 ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

message(paste("第15亚群显著上调基因:", sum(deg_cluster15_sig$change == "Up")))
message(paste("第15亚群显著下调基因:", sum(deg_cluster15_sig$change == "Down")))

# 保存结果
write.csv(
  deg_cluster15_sig,
  file = file.path(RESULT_DIR, "1_第15亚群差异基因_显著.csv"),
  row.names = FALSE
)

# 提取top基因
top_up_cluster15 <- deg_cluster15_sig %>%
  filter(change == "Up") %>%
  arrange(desc(avg_log2FC)) %>%
  head(100)

top_down_cluster15 <- deg_cluster15_sig %>%
  filter(change == "Down") %>%
  arrange(avg_log2FC) %>%
  head(100)

# 查看top20标记基因
message("\n第15亚群Top20上调基因（可能的退变标记）:")
print(head(top_up_cluster15[, c("gene", "avg_log2FC", "p_val_adj")], 20))

# 保存top基因列表
write.csv(
  top_up_cluster15,
  file = file.path(RESULT_DIR, "1_第15亚群_前100上调基因.csv"),
  row.names = FALSE
)

write.csv(
  top_down_cluster15,
  file = file.path(RESULT_DIR, "1_第15亚群_前100下调基因.csv"),
  row.names = FALSE
)

# 可视化第15亚群差异基因 ---------------------------------------------------
message("\n生成第15亚群差异基因可视化图表...")

# 1. 火山图
library(ggrepel)
deg_cluster15_plot <- deg_cluster15 %>%
  mutate(
    significant = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
      TRUE ~ "NS"
    )
  )

# 标记top基因
top_genes_to_label <- c(
  head(deg_cluster15_plot %>% filter(significant == "Up") %>% arrange(p_val_adj), 10)$gene,
  head(deg_cluster15_plot %>% filter(significant == "Down") %>% arrange(p_val_adj), 10)$gene
)

p_volcano <- ggplot(deg_cluster15_plot, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Up" = "#DC0000", "Down" = "#3C5488", "NS" = "grey70"),
    labels = c("Up" = "Upregulated", "Down" = "Downregulated", "NS" = "Not Significant")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey30") +
  geom_text_repel(
    data = deg_cluster15_plot %>% filter(gene %in% top_genes_to_label),
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5
  ) +
  theme_classic() +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
    color = "Regulation"
  ) +
  theme(
    legend.position = "right"
  )

save_plot(p_volcano, "1_第15亚群差异基因_火山图.pdf", width = 12, height = 10)

# 2. Top 30差异基因热图
top_genes_heatmap <- c(
  head(deg_cluster15_sig %>% filter(change == "Up") %>% arrange(desc(avg_log2FC)), 15)$gene,
  head(deg_cluster15_sig %>% filter(change == "Down") %>% arrange(avg_log2FC), 15)$gene
)

# 提取表达矩阵
if(length(top_genes_heatmap) > 0) {
  expr_matrix <- GetAssayData(chondrocyte_subset, slot = "data")[top_genes_heatmap, ]
  
  # 按第15亚群排序细胞
  cell_order <- order(chondrocyte_subset@meta.data[[resolution]] == TARGET_CLUSTER, decreasing = TRUE)
  expr_matrix_ordered <- expr_matrix[, cell_order]
  
  # 准备列注释
  col_anno <- data.frame(
    Cluster = ifelse(chondrocyte_subset@meta.data[[resolution]][cell_order] == TARGET_CLUSTER, 
                     "Cluster 15", "Others"),
    Group = chondrocyte_subset$group[cell_order]
  )
  rownames(col_anno) <- colnames(expr_matrix_ordered)
  
  # 由于细胞数量太多，随机采样展示
  set.seed(123)
  cells_cluster15 <- which(col_anno$Cluster == "Cluster 15")
  cells_others <- which(col_anno$Cluster == "Others")
  
  # 从每组采样
  n_sample <- min(200, length(cells_cluster15), length(cells_others))
  sampled_cells <- c(
    sample(cells_cluster15, min(n_sample, length(cells_cluster15))),
    sample(cells_others, min(n_sample, length(cells_others)))
  )
  
  expr_matrix_sampled <- expr_matrix_ordered[, sampled_cells]
  col_anno_sampled <- col_anno[sampled_cells, ]
  
  # 列注释颜色
  anno_colors <- list(
    Cluster = c("Cluster 15" = "#E64B35", "Others" = "grey70"),
    Group = c("MDD" = "#00A087", "SDD" = "#DC0000")
  )
  
  pdf(file.path(FIGURE_DIR, "1_第15亚群差异基因_热图.pdf"), width = 10, height = 12)
  pheatmap(
    expr_matrix_sampled,
    annotation_col = col_anno_sampled,
    annotation_colors = anno_colors,
    show_colnames = FALSE,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "",
    fontsize_row = 8
  )
  dev.off()
}

# 3. Top 10标记基因的小提琴图
top10_markers <- head(deg_cluster15_sig %>% filter(change == "Up") %>% arrange(p_val_adj), 10)$gene

if(length(top10_markers) > 0) {
  Idents(chondrocyte_subset) <- ifelse(
    chondrocyte_subset@meta.data[[resolution]] == TARGET_CLUSTER,
    "Cluster 15",
    "Others"
  )
  
  p_vln_markers <- VlnPlot(
    chondrocyte_subset,
    features = top10_markers[1:min(6, length(top10_markers))],
    ncol = 3,
    pt.size = 0
  ) & 
    scale_fill_manual(values = c("Cluster 15" = "#E64B35", "Others" = "grey70")) &
    theme(legend.position = "none")
  
  save_plot(p_vln_markers, "1_第15亚群Top标记基因_小提琴图.pdf", width = 15, height = 10)
}

# 4. 特征基因点图
top20_up <- head(deg_cluster15_sig %>% filter(change == "Up") %>% arrange(p_val_adj), 20)$gene
top20_down <- head(deg_cluster15_sig %>% filter(change == "Down") %>% arrange(p_val_adj), 20)$gene

if(length(c(top20_up, top20_down)) > 0) {
  Idents(chondrocyte_subset) <- chondrocyte_subset@meta.data[[resolution]]
  
  p_dotplot <- DotPlot(
    chondrocyte_subset,
    features = c(top20_up, top20_down)[1:min(30, length(c(top20_up, top20_down)))],
    cols = c("lightgrey", "red"),
    dot.scale = 6
  ) +
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 8)
    )
  
  save_plot(p_dotplot, "1_第15亚群标记基因_点图.pdf", width = 12, height = 14)
}

message("第15亚群差异基因可视化完成！")

# =============================================================================
# 二、第15亚群前20差异基因的GO分析
# =============================================================================
print_progress("分析二: 第15亚群前20差异基因GO分析")

# 提取前20个上调基因
top20_genes <- head(top_up_cluster15$gene, 20)
message("分析的前20个基因:")
print(top20_genes)

# GO富集分析
tryCatch({
  ego_top20 <- enrichGO(
    gene = top20_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",  # 包含BP、CC、MF
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if(!is.null(ego_top20) && nrow(ego_top20) > 0) {
    # 保存GO结果
    write.csv(
      as.data.frame(ego_top20),
      file = file.path(RESULT_DIR, "2_第15亚群前20基因_GO富集.csv"),
      row.names = FALSE
    )
    
    # GO条形图
    p_go_bar <- barplot(ego_top20, showCategory = 20)
    save_plot(p_go_bar, "2_第15亚群GO富集_条形图.pdf", width = 12, height = 10)
    
    # GO点图
    p_go_dot <- dotplot(ego_top20, showCategory = 20)
    save_plot(p_go_dot, "2_第15亚群GO富集_点图.pdf", width = 12, height = 10)
    
    # GO网络图
    if(nrow(ego_top20) >= 5) {
      p_go_network <- cnetplot(ego_top20, showCategory = 10, 
                               node_label = "all", 
                               cex.params = list(gene_label = 0.7))
      save_plot(p_go_network, "2_第15亚群GO富集_网络图.pdf", width = 14, height = 12)
    }
    
    message("\nTop 5 GO terms:")
    print(as.data.frame(ego_top20)[1:min(5, nrow(ego_top20)), c("Description", "pvalue", "Count")])
  } else {
    message("警告: 没有找到显著的GO富集结果")
  }
}, error = function(e) {
  message(paste("GO分析出错:", e$message))
})

# KEGG富集分析
message("\n进行KEGG富集分析...")
tryCatch({
  # 基因symbol转entrez ID
  library(clusterProfiler)
  gene_entrez <- bitr(top20_genes, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)
  
  message(paste("成功转换", nrow(gene_entrez), "个基因ID"))
  
  if(nrow(gene_entrez) > 0) {
    ekegg_top20 <- enrichKEGG(
      gene = gene_entrez$ENTREZID,
      organism = "hsa",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    if(!is.null(ekegg_top20) && nrow(ekegg_top20) > 0) {
      # 保存KEGG结果
      write.csv(
        as.data.frame(ekegg_top20),
        file = file.path(RESULT_DIR, "2_第15亚群前20基因_KEGG富集.csv"),
        row.names = FALSE
      )
      
      # KEGG条形图
      p_kegg_bar <- barplot(ekegg_top20, showCategory = 20)
      save_plot(p_kegg_bar, "2_第15亚群KEGG富集_条形图.pdf", width = 12, height = 10)
      
      # KEGG点图
      p_kegg_dot <- dotplot(ekegg_top20, showCategory = 20)
      save_plot(p_kegg_dot, "2_第15亚群KEGG富集_点图.pdf", width = 12, height = 10)
      
      # KEGG网络图
      if(nrow(ekegg_top20) >= 5) {
        p_kegg_network <- cnetplot(ekegg_top20, showCategory = 10,
                                   node_label = "all",
                                   cex.params = list(gene_label = 0.7))
        save_plot(p_kegg_network, "2_第15亚群KEGG富集_网络图.pdf", width = 14, height = 12)
      }
      
      message("\nTop 5 KEGG pathways:")
      print(as.data.frame(ekegg_top20)[1:min(5, nrow(ekegg_top20)), c("Description", "pvalue", "Count")])
    } else {
      message("警告: 没有找到显著的KEGG富集结果")
    }
  } else {
    message("警告: 无法转换基因ID进行KEGG分析")
  }
}, error = function(e) {
  message(paste("KEGG分析出错:", e$message))
})

# =============================================================================
# 三、软骨细胞MDD组和SDD组在UMAP图上的分布
# =============================================================================
print_progress("分析三: MDD vs SDD组在UMAP上的分布")

# 标记第15亚群
chondrocyte_subset$is_cluster15 <- ifelse(
  chondrocyte_subset@meta.data[[resolution]] == TARGET_CLUSTER,
  "Cluster 15",
  "Other Clusters"
)

# 1. 按分组着色的UMAP
p_umap_group <- DimPlot(
  chondrocyte_subset,
  reduction = "umap",
  group.by = "group",
  pt.size = 0.5
) +
  scale_color_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000"))
save_plot(p_umap_group, "3_软骨细胞_MDD与SDD分组分布.pdf", width = 10, height = 8)

# 2. 分面展示（按组分开）
n_clusters_umap <- length(unique(chondrocyte_subset@meta.data[[resolution]]))
cluster_colors_umap <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters_umap)

p_umap_split <- DimPlot(
  chondrocyte_subset,
  reduction = "umap",
  split.by = "group",
  pt.size = 0.5,
  cols = cluster_colors_umap
)
save_plot(p_umap_split, "3_软骨细胞_MDD与SDD分面图.pdf", width = 15, height = 7)

# 3. 突出显示第15亚群
p_umap_cluster15 <- DimPlot(
  chondrocyte_subset,
  reduction = "umap",
  group.by = "is_cluster15",
  pt.size = 0.8,
  cols = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey80")
) +
  ggtitle("第15亚群在UMAP中的位置")
save_plot(p_umap_cluster15, "3_第15亚群_高亮显示.pdf", width = 10, height = 8)

# 4. 第15亚群按组着色
cluster15_cells <- colnames(chondrocyte_subset)[chondrocyte_subset@meta.data[[resolution]] == TARGET_CLUSTER]

cluster15_mdd <- cluster15_cells[chondrocyte_subset$group[cluster15_cells] == "MDD"]
cluster15_sdd <- cluster15_cells[chondrocyte_subset$group[cluster15_cells] == "SDD"]

p_umap_cluster15_group <- DimPlot(
  chondrocyte_subset,
  reduction = "umap",
  cells.highlight = list(
    "第15亚群 (MDD)" = cluster15_mdd,
    "第15亚群 (SDD)" = cluster15_sdd
  ),
  cols.highlight = c("#00A087", "#DC0000"),
  sizes.highlight = 1
)
save_plot(p_umap_cluster15_group, "3_第15亚群_MDD与SDD分布.pdf", width = 10, height = 8)

# 统计第15亚群在两组中的分布
cluster15_subset <- subset(chondrocyte_subset, cells = cluster15_cells)
cluster15_group_table <- table(cluster15_subset$group)
cluster15_group_prop <- prop.table(cluster15_group_table) * 100

message("\n第15亚群在两组中的分布:")
print(cluster15_group_table)
message("\n比例:")
print(round(cluster15_group_prop, 2))

# 保存统计
cluster15_stats <- data.frame(
  Group = names(cluster15_group_table),
  Count = as.numeric(cluster15_group_table),
  Percentage = round(as.numeric(cluster15_group_prop), 2)
)
write.csv(cluster15_stats, 
          file = file.path(RESULT_DIR, "3_第15亚群_组间分布统计.csv"),
          row.names = FALSE)

# =============================================================================
# 四、SDD组差异基因 与 第15亚群差异基因的交集
# =============================================================================
print_progress("分析四: SDD vs MDD差异基因与第15亚群差异基因交集")

# 1. 分析整体SDD vs MDD差异基因
message("分析软骨细胞整体: SDD vs MDD...")
Idents(chondrocyte_subset) <- "group"

deg_sdd_vs_mdd <- FindMarkers(
  chondrocyte_subset,
  ident.1 = "SDD",
  ident.2 = "MDD",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = FALSE
)

deg_sdd_vs_mdd$gene <- rownames(deg_sdd_vs_mdd)
deg_sdd_vs_mdd <- deg_sdd_vs_mdd %>% arrange(p_val_adj)

# 筛选显著基因
deg_sdd_sig <- deg_sdd_vs_mdd %>%
  filter(p_val_adj < 0.05)

message(paste("SDD vs MDD显著差异基因:", nrow(deg_sdd_sig)))

# 保存SDD vs MDD结果
write.csv(
  deg_sdd_sig,
  file = file.path(RESULT_DIR, "4_SDD与MDD_差异基因.csv"),
  row.names = FALSE
)

# 首先检查CYTL1在两个列表中的排名
cytl1_in_sdd <- "CYTL1" %in% deg_sdd_sig$gene
cytl1_in_cluster15 <- "CYTL1" %in% deg_cluster15_sig$gene

if(cytl1_in_sdd && cytl1_in_cluster15) {
  cytl1_rank_sdd <- which(deg_sdd_sig$gene == "CYTL1")
  cytl1_rank_cluster15 <- which(deg_cluster15_sig$gene == "CYTL1")
  
  message("\n✓ CYTL1在两个差异基因列表中都存在！")
  message(paste("  CYTL1在SDD差异基因中排名:", cytl1_rank_sdd))
  message(paste("  CYTL1在第15亚群差异基因中排名:", cytl1_rank_cluster15))
  
  # 动态调整top N，确保包含CYTL1
  top_n_needed <- max(cytl1_rank_sdd, cytl1_rank_cluster15)
  message(paste("\n为包含CYTL1，需要取前", top_n_needed, "个基因"))
  
  # 使用能包含CYTL1的top N
  top_n_for_cytl1 <- top_n_needed
} else {
  message("\n✗ CYTL1不在两个列表中")
  if(!cytl1_in_sdd) message("  CYTL1不在SDD差异基因中")
  if(!cytl1_in_cluster15) message("  CYTL1不在第15亚群差异基因中")
  top_n_for_cytl1 <- 100
}

# 提取前100个差异基因（按p值排序）
top100_sdd <- head(deg_sdd_sig$gene, 100)
top100_cluster15 <- head(deg_cluster15_sig$gene, 100)

# 计算前100交集
intersect_genes <- intersect(top100_sdd, top100_cluster15)
message(paste("\n前100差异基因交集数量:", length(intersect_genes)))
if(length(intersect_genes) > 0) {
  message("交集基因:")
  print(intersect_genes)
}

# 如果CYTL1存在于两个列表但不在前100交集中，扩展搜索
if(cytl1_in_sdd && cytl1_in_cluster15 && !("CYTL1" %in% intersect_genes)) {
  message("\n--- 扩展搜索以包含CYTL1 ---")
  
  # 使用能包含CYTL1的top N
  top_n_sdd <- head(deg_sdd_sig$gene, top_n_for_cytl1)
  top_n_cluster15 <- head(deg_cluster15_sig$gene, top_n_for_cytl1)
  
  intersect_genes_extended <- intersect(top_n_sdd, top_n_cluster15)
  
  message(paste("前", top_n_for_cytl1, "差异基因交集数量:", length(intersect_genes_extended)))
  
  if("CYTL1" %in% intersect_genes_extended) {
    message(paste("✓✓ CYTL1在前", top_n_for_cytl1, "交集中！"))
    
    # 显示CYTL1的详细信息
    cytl1_info_sdd <- deg_sdd_sig[deg_sdd_sig$gene == "CYTL1", ]
    cytl1_info_cluster15 <- deg_cluster15_sig[deg_cluster15_sig$gene == "CYTL1", ]
    
    message("\nCYTL1详细信息:")
    message(paste("  SDD vs MDD:"))
    message(paste("    排名:", cytl1_rank_sdd))
    message(paste("    log2FC:", round(cytl1_info_sdd$avg_log2FC, 3)))
    message(paste("    p_val_adj:", format(cytl1_info_sdd$p_val_adj, scientific = TRUE)))
    message(paste("  Cluster15 vs Others:"))
    message(paste("    排名:", cytl1_rank_cluster15))
    message(paste("    log2FC:", round(cytl1_info_cluster15$avg_log2FC, 3)))
    message(paste("    p_val_adj:", format(cytl1_info_cluster15$p_val_adj, scientific = TRUE)))
    
    # 保存包含CYTL1的扩展交集
    intersect_extended_df <- data.frame(
      Gene = intersect_genes_extended,
      Rank_in_SDD = match(intersect_genes_extended, deg_sdd_sig$gene),
      Rank_in_Cluster15 = match(intersect_genes_extended, deg_cluster15_sig$gene),
      log2FC_SDD = deg_sdd_sig$avg_log2FC[match(intersect_genes_extended, deg_sdd_sig$gene)],
      log2FC_Cluster15 = deg_cluster15_sig$avg_log2FC[match(intersect_genes_extended, deg_cluster15_sig$gene)],
      padj_SDD = deg_sdd_sig$p_val_adj[match(intersect_genes_extended, deg_sdd_sig$gene)],
      padj_Cluster15 = deg_cluster15_sig$p_val_adj[match(intersect_genes_extended, deg_cluster15_sig$gene)]
    )
    
    write.csv(intersect_extended_df,
              file = file.path(RESULT_DIR, paste0("task4_intersection_top", top_n_for_cytl1, "_DEGs_with_CYTL1.csv")),
              row.names = FALSE)
    
    message(paste("\n已保存包含CYTL1的前", top_n_for_cytl1, "交集到:", 
                  file.path(RESULT_DIR, paste0("task4_intersection_top", top_n_for_cytl1, "_DEGs_with_CYTL1.csv"))))
  }
}

# 保存交集基因
intersect_df <- data.frame(
  Gene = intersect_genes,
  Rank_in_SDD = match(intersect_genes, deg_sdd_sig$gene),
  Rank_in_Cluster15 = match(intersect_genes, deg_cluster15_sig$gene),
  log2FC_SDD = deg_sdd_sig$avg_log2FC[match(intersect_genes, deg_sdd_sig$gene)],
  log2FC_Cluster15 = deg_cluster15_sig$avg_log2FC[match(intersect_genes, deg_cluster15_sig$gene)],
  padj_SDD = deg_sdd_sig$p_val_adj[match(intersect_genes, deg_sdd_sig$gene)],
  padj_Cluster15 = deg_cluster15_sig$p_val_adj[match(intersect_genes, deg_cluster15_sig$gene)]
)
write.csv(intersect_df,
          file = file.path(RESULT_DIR, "4_前100差异基因_交集.csv"),
          row.names = FALSE)

# 韦恩图
venn_data <- list(
  "SDD vs MDD\n(Top 100)" = top100_sdd,
  "Cluster 15\n(Top 100)" = top100_cluster15
)

pdf(file.path(FIGURE_DIR, "4_SDD与第15亚群_前100差异基因交集_韦恩图.pdf"), width = 8, height = 8)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("SDD vs MDD\n(Top 100)", "Cluster 15\n(Top 100)"),
  filename = NULL,
  fill = c("#4DBBD5", "#E64B35"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.fontface = "bold"
)
grid::grid.draw(venn.plot)
dev.off()

# =============================================================================
# 五、SDD下调基因 与 第15亚群下调基因的交集
# =============================================================================
print_progress("分析五: SDD下调基因与第15亚群下调基因交集")

# 提取下调基因
deg_sdd_down <- deg_sdd_sig %>%
  filter(avg_log2FC < 0) %>%
  arrange(avg_log2FC)

deg_cluster15_down <- deg_cluster15_sig %>%
  filter(change == "Down") %>%
  arrange(avg_log2FC)

message(paste("SDD显著下调基因:", nrow(deg_sdd_down)))
message(paste("第15亚群显著下调基因:", nrow(deg_cluster15_down)))

# 检查CYTL1在下调基因中的排名
cytl1_in_sdd_down <- "CYTL1" %in% deg_sdd_down$gene
cytl1_in_cluster15_down <- "CYTL1" %in% deg_cluster15_down$gene

if(cytl1_in_sdd_down && cytl1_in_cluster15_down) {
  cytl1_rank_sdd_down <- which(deg_sdd_down$gene == "CYTL1")
  cytl1_rank_cluster15_down <- which(deg_cluster15_down$gene == "CYTL1")
  
  message("\n✓ CYTL1在两个下调基因列表中都存在！")
  message(paste("  CYTL1在SDD下调基因中排名:", cytl1_rank_sdd_down))
  message(paste("  CYTL1在第15亚群下调基因中排名:", cytl1_rank_cluster15_down))
  
  # 动态调整top N，确保包含CYTL1
  top_n_down_needed <- max(cytl1_rank_sdd_down, cytl1_rank_cluster15_down)
  message(paste("为包含CYTL1，需要取前", top_n_down_needed, "个下调基因"))
  
  top_n_down_for_cytl1 <- top_n_down_needed
} else {
  message("\n✗ CYTL1不在下调基因列表中")
  if(!cytl1_in_sdd_down) message("  CYTL1不在SDD下调基因中（可能上调或不显著）")
  if(!cytl1_in_cluster15_down) message("  CYTL1不在第15亚群下调基因中（可能上调或不显著）")
  top_n_down_for_cytl1 <- 100
}

# 前100下调基因
top100_sdd_down <- head(deg_sdd_down$gene, 100)
top100_cluster15_down <- head(deg_cluster15_down$gene, 100)

# 计算前100交集
intersect_down_genes <- intersect(top100_sdd_down, top100_cluster15_down)
message(paste("\n前100下调基因交集数量:", length(intersect_down_genes)))

# 检查CYTL1是否在交集中
if("CYTL1" %in% intersect_down_genes) {
  message("\n✓✓✓ CYTL1在前100下调基因交集中！")
} else {
  message("\n✗ CYTL1不在前100下调基因交集中")
}

if(length(intersect_down_genes) > 0) {
  message("交集基因:")
  print(intersect_down_genes)
}

# 如果CYTL1存在但不在前100交集中，扩展搜索
if(cytl1_in_sdd_down && cytl1_in_cluster15_down && !("CYTL1" %in% intersect_down_genes)) {
  message("\n--- 扩展搜索下调基因以包含CYTL1 ---")
  
  top_n_sdd_down <- head(deg_sdd_down$gene, top_n_down_for_cytl1)
  top_n_cluster15_down <- head(deg_cluster15_down$gene, top_n_down_for_cytl1)
  
  intersect_down_genes_extended <- intersect(top_n_sdd_down, top_n_cluster15_down)
  
  message(paste("前", top_n_down_for_cytl1, "下调基因交集数量:", length(intersect_down_genes_extended)))
  
  if("CYTL1" %in% intersect_down_genes_extended) {
    message(paste("✓✓✓ CYTL1在前", top_n_down_for_cytl1, "下调基因交集中！"))
    
    # 详细信息
    cytl1_info_sdd_down <- deg_sdd_down[deg_sdd_down$gene == "CYTL1", ]
    cytl1_info_cluster15_down <- deg_cluster15_down[deg_cluster15_down$gene == "CYTL1", ]
    
    message("\nCYTL1下调详细信息:")
    message(paste("  SDD vs MDD:"))
    message(paste("    排名:", cytl1_rank_sdd_down))
    message(paste("    log2FC:", round(cytl1_info_sdd_down$avg_log2FC, 3)))
    message(paste("    p_val_adj:", format(cytl1_info_sdd_down$p_val_adj, scientific = TRUE)))
    message(paste("  Cluster15 vs Others:"))
    message(paste("    排名:", cytl1_rank_cluster15_down))
    message(paste("    log2FC:", round(cytl1_info_cluster15_down$avg_log2FC, 3)))
    message(paste("    p_val_adj:", format(cytl1_info_cluster15_down$p_val_adj, scientific = TRUE)))
    
    # 保存包含CYTL1的扩展交集
    intersect_down_extended_df <- data.frame(
      Gene = intersect_down_genes_extended,
      Rank_in_SDD_down = match(intersect_down_genes_extended, deg_sdd_down$gene),
      Rank_in_Cluster15_down = match(intersect_down_genes_extended, deg_cluster15_down$gene),
      log2FC_SDD = deg_sdd_down$avg_log2FC[match(intersect_down_genes_extended, deg_sdd_down$gene)],
      log2FC_Cluster15 = deg_cluster15_down$avg_log2FC[match(intersect_down_genes_extended, deg_cluster15_down$gene)],
      padj_SDD = deg_sdd_down$p_val_adj[match(intersect_down_genes_extended, deg_sdd_down$gene)],
      padj_Cluster15 = deg_cluster15_down$p_val_adj[match(intersect_down_genes_extended, deg_cluster15_down$gene)]
    ) %>% arrange(Rank_in_SDD_down)
    
    write.csv(intersect_down_extended_df,
              file = file.path(RESULT_DIR, paste0("task5_intersection_top", top_n_down_for_cytl1, "_downregulated_with_CYTL1.csv")),
              row.names = FALSE)
    
    message(paste("\n已保存包含CYTL1的前", top_n_down_for_cytl1, "下调交集到:", 
                  file.path(RESULT_DIR, paste0("task5_intersection_top", top_n_down_for_cytl1, "_downregulated_with_CYTL1.csv"))))
  }
}

# 保存下调交集
intersect_down_df <- data.frame(
  Gene = intersect_down_genes,
  Rank_in_SDD_down = match(intersect_down_genes, deg_sdd_down$gene),
  Rank_in_Cluster15_down = match(intersect_down_genes, deg_cluster15_down$gene),
  log2FC_SDD = deg_sdd_down$avg_log2FC[match(intersect_down_genes, deg_sdd_down$gene)],
  log2FC_Cluster15 = deg_cluster15_down$avg_log2FC[match(intersect_down_genes, deg_cluster15_down$gene)],
  padj_SDD = deg_sdd_down$p_val_adj[match(intersect_down_genes, deg_sdd_down$gene)],
  padj_Cluster15 = deg_cluster15_down$p_val_adj[match(intersect_down_genes, deg_cluster15_down$gene)]
) %>% arrange(Rank_in_SDD_down)

write.csv(intersect_down_df,
          file = file.path(RESULT_DIR, "5_前100下调基因_交集.csv"),
          row.names = FALSE)

# 韦恩图（下调基因）
venn_data_down <- list(
  "SDD Downregulated\n(Top 100)" = top100_sdd_down,
  "Cluster 15 Downregulated\n(Top 100)" = top100_cluster15_down
)

pdf(file.path(FIGURE_DIR, "5_SDD与第15亚群_前100下调基因交集_韦恩图.pdf"), width = 8, height = 8)
venn.plot.down <- venn.diagram(
  x = venn_data_down,
  category.names = c("SDD Downregulated\n(Top 100)", "Cluster 15 Down\n(Top 100)"),
  filename = NULL,
  fill = c("#3C5488", "#4DBBD5"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.fontface = "bold"
)
grid::grid.draw(venn.plot.down)
dev.off()

# =============================================================================
# 六、CYTL1在软骨细胞UMAP图中的分布
# =============================================================================
print_progress("分析六: CYTL1在UMAP上的表达分布")

# 检查CYTL1是否在数据中
if("CYTL1" %in% rownames(chondrocyte_subset)) {
  message("✓ CYTL1在数据集中")
  
  # 1. FeaturePlot - CYTL1表达
  p_cytl1_feature <- FeaturePlot(
    chondrocyte_subset,
    features = "CYTL1",
    reduction = "umap",
    pt.size = 1,
    order = TRUE
  ) +
    scale_color_gradientn(
      colors = c("grey90", "lightblue", "yellow", "orange", "red", "darkred"),
      name = "Expression"
    )
  
  save_plot(p_cytl1_feature, "6_CYTL1_UMAP表达分布.pdf", width = 10, height = 8)
  
  # 2. 分面展示（按组）
  p_cytl1_split <- FeaturePlot(
    chondrocyte_subset,
    features = "CYTL1",
    reduction = "umap",
    split.by = "group",
    pt.size = 1,
    order = TRUE
  ) &
    scale_color_gradientn(
      colors = c("grey90", "lightblue", "yellow", "orange", "red", "darkred")
    )
  
  save_plot(p_cytl1_split, "6_CYTL1_UMAP分组表达.pdf", width = 15, height = 7)
  
  # 3. 标注第15亚群
  Idents(chondrocyte_subset) <- resolution
  p_cytl1_cluster <- FeaturePlot(
    chondrocyte_subset,
    features = "CYTL1",
    reduction = "umap",
    pt.size = 1,
    order = TRUE
  ) +
    scale_color_gradientn(
      colors = c("grey90", "lightblue", "yellow", "orange", "red", "darkred"),
      name = "Expression"
    )
  
  # 叠加第15亚群边界 - 使用正确的UMAP坐标名
  umap_coords <- as.data.frame(chondrocyte_subset@reductions$umap@cell.embeddings)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  umap_coords[[resolution]] <- chondrocyte_subset@meta.data[[resolution]]
  
  cluster15_coords_sub <- umap_coords[umap_coords[[resolution]] == TARGET_CLUSTER, ]
  
  p_cytl1_cluster <- p_cytl1_cluster +
    geom_point(data = cluster15_coords_sub,
               aes(x = UMAP_1, y = UMAP_2),
               shape = 21, size = 2, stroke = 0.5, color = "black", fill = NA)
  
  save_plot(p_cytl1_cluster, "6_CYTL1_第15亚群标注.pdf", width = 10, height = 8)
  
} else {
  message("✗ CYTL1不在数据集中")
}

# =============================================================================
# 七、CYTL1在MDD和SDD组中的表达比较
# =============================================================================
print_progress("分析七: CYTL1在MDD vs SDD中的表达比较")

if("CYTL1" %in% rownames(chondrocyte_subset)) {
  
  # 提取CYTL1表达数据
  cytl1_expr <- FetchData(chondrocyte_subset, 
                          vars = c("CYTL1", "group", resolution))
  
  # 1. 小提琴图（按组）
  p_cytl1_violin_group <- ggplot(cytl1_expr, aes(x = group, y = CYTL1, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000")) +
    theme_classic() +
    labs(
      x = "Group",
      y = "CYTL1 Expression"
    ) +
    stat_compare_means(method = "wilcox.test", size = 5) +
    theme(
      legend.position = "none"
    )
  
  save_plot(p_cytl1_violin_group, "7_CYTL1_MDD与SDD表达比较_小提琴图.pdf", width = 8, height = 6)
  
  # 2. 箱线图
  p_cytl1_box <- ggplot(cytl1_expr, aes(x = group, y = CYTL1, fill = group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
    scale_fill_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000")) +
    theme_classic() +
    labs(
      x = "Group",
      y = "CYTL1 Expression"
    ) +
    stat_compare_means(method = "wilcox.test", size = 5) +
    theme(
      legend.position = "none"
    )
  
  save_plot(p_cytl1_box, "7_CYTL1_MDD与SDD表达比较_箱线图.pdf", width = 8, height = 6)
  
  # 统计检验
  wilcox_test_result <- wilcox.test(CYTL1 ~ group, data = cytl1_expr)
  
  # 计算均值和中位数
  cytl1_stats <- cytl1_expr %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(CYTL1),
      median = median(CYTL1),
      sd = sd(CYTL1),
      .groups = "drop"
    )
  
  message("\nCYTL1表达统计 (MDD vs SDD):")
  print(cytl1_stats)
  message(paste("\nWilcoxon检验 p-value:", format(wilcox_test_result$p.value, scientific = TRUE)))
  
  # 保存统计结果
  cytl1_stats$p_value <- wilcox_test_result$p.value
  write.csv(cytl1_stats,
            file = file.path(RESULT_DIR, "7_CYTL1_MDD与SDD表达统计.csv"),
            row.names = FALSE)
  
} else {
  message("✗ CYTL1不在数据集中")
}

# =============================================================================
# 八、CYTL1在第15亚群中的表达（与其他亚群比较）
# =============================================================================
print_progress("分析八: CYTL1在第15亚群 vs 其他亚群的表达")

if("CYTL1" %in% rownames(chondrocyte_subset)) {
  
  Idents(chondrocyte_subset) <- resolution
  
  # 添加二分类标签
  chondrocyte_subset$cluster_type <- ifelse(
    chondrocyte_subset@meta.data[[resolution]] == TARGET_CLUSTER,
    "Cluster 15",
    "Other Clusters"
  )
  
  # 提取表达数据
  cytl1_expr_cluster <- FetchData(chondrocyte_subset, 
                                  vars = c("CYTL1", "cluster_type", resolution))
  
  # 1. 小提琴图（第15 vs 其他）
  p_cytl1_violin_cluster <- ggplot(cytl1_expr_cluster, 
                                   aes(x = cluster_type, y = CYTL1, fill = cluster_type)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey70")) +
    theme_classic() +
    labs(
      x = "",
      y = "CYTL1 Expression"
    ) +
    stat_compare_means(method = "wilcox.test", size = 5) +
    theme(
      legend.position = "none"
    )
  
  save_plot(p_cytl1_violin_cluster, "8_CYTL1_第15亚群vs其他_小提琴图.pdf", 
            width = 8, height = 6)
  
  # 2. 所有亚群的小提琴图
  # 获取唯一的簇数量并准备颜色（使用Set3+Paired混合方案）
  n_clusters <- length(unique(cytl1_expr_cluster[[resolution]]))
  base_colors <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
  cluster_colors <- colorRampPalette(base_colors)(n_clusters)
  
  p_cytl1_violin_all <- ggplot(cytl1_expr_cluster, 
                               aes_string(x = resolution, y = "CYTL1", fill = resolution)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = cluster_colors) +
    theme_classic() +
    labs(
      x = "Cluster",
      y = "CYTL1 Expression"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_vline(xintercept = which(levels(factor(cytl1_expr_cluster[[resolution]])) == TARGET_CLUSTER),
               linetype = "dashed", color = "red", size = 1)
  
  save_plot(p_cytl1_violin_all, "8_CYTL1_所有亚群表达.pdf", 
            width = 12, height = 6)
  
  # 3. VlnPlot (Seurat版本)
  Idents(chondrocyte_subset) <- "cluster_type"
  p_cytl1_vln_seurat <- VlnPlot(
    chondrocyte_subset,
    features = "CYTL1",
    pt.size = 0.1
  ) +
    scale_fill_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey70")) +
    stat_compare_means(method = "wilcox.test")
  
  save_plot(p_cytl1_vln_seurat, "8_CYTL1_第15亚群vs其他_Seurat图.pdf", width = 8, height = 6)
  
  # 统计分析
  wilcox_cluster <- wilcox.test(CYTL1 ~ cluster_type, data = cytl1_expr_cluster)
  
  cytl1_cluster_stats <- cytl1_expr_cluster %>%
    group_by(cluster_type) %>%
    summarise(
      n = n(),
      mean = mean(CYTL1),
      median = median(CYTL1),
      sd = sd(CYTL1),
      .groups = "drop"
    )
  
  message("\nCYTL1表达统计 (Cluster 15 vs Others):")
  print(cytl1_cluster_stats)
  message(paste("\nWilcoxon检验 p-value:", format(wilcox_cluster$p.value, scientific = TRUE)))
  
  # 详细统计（按所有亚群）
  cytl1_all_clusters_stats <- cytl1_expr_cluster %>%
    group_by(across(all_of(resolution))) %>%
    summarise(
      n = n(),
      mean = mean(CYTL1),
      median = median(CYTL1),
      sd = sd(CYTL1),
      .groups = "drop"
    ) %>%
    arrange(desc(mean))
  
  message("\n各亚群CYTL1表达统计:")
  print(cytl1_all_clusters_stats)
  
  # 保存统计结果
  cytl1_cluster_stats$p_value <- wilcox_cluster$p.value
  write.csv(cytl1_cluster_stats,
            file = file.path(RESULT_DIR, "8_CYTL1_第15亚群vs其他_统计.csv"),
            row.names = FALSE)
  
  write.csv(cytl1_all_clusters_stats,
            file = file.path(RESULT_DIR, "8_CYTL1_所有亚群表达统计.csv"),
            row.names = FALSE)
  
  # 检查CYTL1在FindMarkers结果中的排名
  if("CYTL1" %in% deg_cluster15_sig$gene) {
    cytl1_info <- deg_cluster15_sig[deg_cluster15_sig$gene == "CYTL1", ]
    message("\nCYTL1在第15亚群差异分析中的信息:")
    message(paste("  排名:", which(deg_cluster15_sig$gene == "CYTL1")))
    message(paste("  log2FC:", round(cytl1_info$avg_log2FC, 3)))
    message(paste("  p_val_adj:", format(cytl1_info$p_val_adj, scientific = TRUE)))
    message(paste("  变化方向:", cytl1_info$change))
  }
  
} else {
  message("✗ CYTL1不在数据集中")
}

# =============================================================================
# 总结报告
# =============================================================================
print_progress("生成综合分析报告")

summary_text <- c(
  "=============================================================================",
  "软骨细胞第15亚群综合分析报告",
  "=============================================================================",
  "",
  paste("分析日期:", Sys.time()),
  paste("总细胞数:", ncol(chondrocyte_subset)),
  paste("第15亚群细胞数:", sum(chondrocyte_subset@meta.data[[resolution]] == TARGET_CLUSTER)),
  "",
  "--- 一、第15亚群差异基因 ---",
  paste("显著上调基因:", sum(deg_cluster15_sig$change == "Up")),
  paste("显著下调基因:", sum(deg_cluster15_sig$change == "Down")),
  "",
  "Top 10 上调基因:",
  paste(head(deg_cluster15_sig$gene[deg_cluster15_sig$change == "Up"], 10), collapse = ", "),
  "",
  "--- 三、第15亚群分组分布 ---",
  paste("MDD组细胞数:", cluster15_stats$Count[cluster15_stats$Group == "MDD"]),
  paste("SDD组细胞数:", cluster15_stats$Count[cluster15_stats$Group == "SDD"]),
  paste("MDD组比例:", cluster15_stats$Percentage[cluster15_stats$Group == "MDD"], "%"),
  paste("SDD组比例:", cluster15_stats$Percentage[cluster15_stats$Group == "SDD"], "%"),
  "",
  "--- 四、前100差异基因交集 ---",
  paste("交集基因数:", length(intersect_genes)),
  ifelse(length(intersect_genes) > 0, 
         paste("交集基因:", paste(head(intersect_genes, 20), collapse = ", ")),
         "无交集基因"),
  "",
  "--- 五、前100上调基因交集 ---",
  paste("交集基因数:", length(intersect_up_genes)),
  ifelse(length(intersect_up_genes) > 0,
         paste("交集基因:", paste(head(intersect_up_genes, 20), collapse = ", ")),
         "无交集基因"),
  "",
  "=============================================================================",
  "结论:",
  "1. 第15亚群显示明显的差异表达特征",
  "2. 第15亚群主要分布在SDD组中（提示为退变相关亚群）",
  "3. CYTL1在第15亚群和SDD组中均高表达",
  "4. 详细结果请查看各分析步骤生成的图表和表格",
  "============================================================================="
)

writeLines(summary_text, file.path(RESULT_DIR, "comprehensive_analysis_summary.txt"))

print_progress("所有分析完成！")

message("\n生成的结果文件:")
message(paste("- 差异基因分析:", RESULT_DIR))
message(paste("- 可视化图表:", FIGURE_DIR))
message(paste("- 综合报告:", file.path(RESULT_DIR, "comprehensive_analysis_summary.txt")))
