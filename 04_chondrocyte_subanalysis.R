# =============================================================================
# 软骨细胞亚群深度分析 - 优化版
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")
load_packages()
library(clustree)

print_progress("软骨细胞亚群深度分析")

# 1. 加载软骨细胞数据 --------------------------------------------------------
message("加载软骨细胞数据...")
load(file.path(DATA_DIR, "chondrocyte_cells.Rdata"))

# 检查数据
message(paste("软骨细胞数量:", ncol(chondrocyte_subset)))
message("分组分布:")
print(table(chondrocyte_subset$group))
message("样本分布:")
print(table(chondrocyte_subset$orig.ident))

# 2. 数据预处理 --------------------------------------------------------------
print_progress("步骤1: 数据预处理")

# 标准化
message("标准化数据...")
chondrocyte_subset <- NormalizeData(chondrocyte_subset)

# 查找高变基因
message("查找高变基因...")
chondrocyte_subset <- FindVariableFeatures(
  chondrocyte_subset,
  selection.method = "vst",
  nfeatures = 2000
)

# 数据缩放（只缩放高变基因，加快速度）
message("缩放数据...")
# 只缩放高变基因可以大幅提速，对后续分析影响很小
chondrocyte_subset <- ScaleData(
  chondrocyte_subset,
  features = VariableFeatures(chondrocyte_subset),
  vars.to.regress = "percent.mt"  # 只回归线粒体比例
)

# 3. PCA降维 -----------------------------------------------------------------
print_progress("步骤2: PCA降维")

message("运行PCA...")
chondrocyte_subset <- RunPCA(
  chondrocyte_subset,
  features = VariableFeatures(chondrocyte_subset),
  npcs = 50,
  verbose = FALSE
)

# 可视化PCA结果
p_pca_group <- DimPlot(chondrocyte_subset, reduction = "pca", group.by = "group") +
  scale_color_custom() +
  ggtitle("Chondrocyte PCA - by Group")
save_plot(p_pca_group, "chondrocyte_pca_by_group.pdf")

p_pca_sample <- DimPlot(chondrocyte_subset, reduction = "pca", group.by = "orig.ident") +
  scale_color_custom() +
  ggtitle("Chondrocyte PCA - by Sample")
save_plot(p_pca_sample, "chondrocyte_pca_by_sample.pdf")

# ElbowPlot确定维度
p_elbow <- ElbowPlot(chondrocyte_subset, ndims = 50, reduction = "pca")
save_plot(p_elbow, "chondrocyte_elbow_plot.pdf")

# 根据ElbowPlot选择合适的维度数
n_pcs <- 20  # 可以根据ElbowPlot调整

# 4. 聚类分析 ----------------------------------------------------------------
print_progress("步骤3: 聚类分析")

message("查找邻居...")
chondrocyte_subset <- FindNeighbors(
  chondrocyte_subset,
  dims = 1:n_pcs,
  verbose = FALSE
)

message("多分辨率聚类...")
chondrocyte_subset <- FindClusters(
  chondrocyte_subset,
  resolution = seq(from = 0.1, to = 2.0, by = 0.1),
  verbose = FALSE
)

# 5. UMAP和tSNE降维 ----------------------------------------------------------
print_progress("步骤4: UMAP和tSNE降维")

message("运行UMAP...")
chondrocyte_subset <- RunUMAP(
  chondrocyte_subset,
  dims = 1:n_pcs,
  verbose = FALSE
)

message("运行tSNE...")
chondrocyte_subset <- RunTSNE(
  chondrocyte_subset,
  dims = 1:n_pcs,
  verbose = FALSE
)

# 保存中间结果
save(chondrocyte_subset, 
     file = file.path(DATA_DIR, "chondrocyte_clustered.Rdata"))

# 6. 聚类可视化 --------------------------------------------------------------
print_progress("步骤5: 聚类可视化")

# 聚类树
message("生成聚类树...")
p_clustree <- clustree(chondrocyte_subset, prefix = "RNA_snn_res.")
save_plot(p_clustree, "chondrocyte_clustree.pdf", width = 12, height = 12)

# 选择最佳分辨率（根据聚类树和生物学意义）
best_res <- "RNA_snn_res.0.5"  # 可以根据clustree调整
Idents(chondrocyte_subset) <- best_res

message(paste("使用分辨率:", best_res))
message("聚类分布:")
print(table(chondrocyte_subset@meta.data[[best_res]]))

# UMAP可视化
p_umap_cluster <- DimPlot(
  chondrocyte_subset,
  reduction = "umap",
  label = TRUE,
  label.size = 5
) +
  ggtitle(paste("Chondrocyte Subclusters -", best_res))
save_plot(p_umap_cluster, "chondrocyte_umap_clusters.pdf")

p_umap_group <- DimPlot(
  chondrocyte_subset,
  reduction = "umap",
  group.by = "group"
) +
  scale_color_custom() +
  ggtitle("Chondrocyte UMAP - by Group")
save_plot(p_umap_group, "chondrocyte_umap_by_group.pdf")

p_umap_sample <- DimPlot(
  chondrocyte_subset,
  reduction = "umap",
  group.by = "orig.ident"
) +
  scale_color_custom() +
  ggtitle("Chondrocyte UMAP - by Sample")
save_plot(p_umap_sample, "chondrocyte_umap_by_sample.pdf")

# tSNE可视化
p_tsne_cluster <- DimPlot(
  chondrocyte_subset,
  reduction = "tsne",
  label = TRUE,
  label.size = 5
) +
  ggtitle(paste("Chondrocyte Subclusters (tSNE) -", best_res))
save_plot(p_tsne_cluster, "chondrocyte_tsne_clusters.pdf")

# 组合图展示
p_combined <- (p_umap_sample + p_umap_group) / (p_tsne_cluster + p_umap_cluster)
save_plot(p_combined, "chondrocyte_clustering_overview.pdf", width = 15, height = 15)

# 按组分面
p_split <- DimPlot(
  chondrocyte_subset,
  reduction = "umap",
  split.by = "group",
  label = TRUE
) +
  ggtitle("Chondrocyte Subclusters: SDD vs MDD")
save_plot(p_split, "chondrocyte_umap_split_by_group.pdf", width = 15, height = 8)

# 7. 亚群标记基因分析 --------------------------------------------------------
print_progress("步骤6: 亚群标记基因分析")

message("查找各亚群的标记基因...")
Idents(chondrocyte_subset) <- best_res

# 查找所有亚群的标记基因
all_markers <- FindAllMarkers(
  chondrocyte_subset,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

# 保存标记基因
write.csv(
  all_markers,
  file = file.path(RESULT_DIR, "chondrocyte_subcluster_markers.csv"),
  row.names = FALSE
)

# 提取每个亚群的top10标记基因
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(
  top10_markers,
  file = file.path(RESULT_DIR, "chondrocyte_subcluster_top10_markers.csv"),
  row.names = FALSE
)

# DotPlot展示top标记基因
top5_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene)

if(length(top5_per_cluster) > 0) {
  p_dotplot <- DotPlot(
    chondrocyte_subset,
    features = unique(top5_per_cluster),
    cols = c("lightgrey", "red")
  ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Top 5 Markers per Subcluster")
  
  save_plot(p_dotplot, "chondrocyte_subcluster_markers_dotplot.pdf", 
            width = 20, height = 8)
}

# 热图展示top标记基因
top3_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  pull(gene)

if(length(top3_per_cluster) > 0) {
  p_heatmap <- DoHeatmap(
    chondrocyte_subset,
    features = unique(top3_per_cluster),
    size = 3
  ) +
    scale_fill_gradientn(colors = c("#4DBBD5", "white", "#E64B35"))
  
  save_plot(p_heatmap, "chondrocyte_subcluster_markers_heatmap.pdf", 
            width = 12, height = 15)
}

# 8. 亚群比例分析 ------------------------------------------------------------
print_progress("步骤7: 亚群比例分析")

# 统计各亚群在两组中的分布
subcluster_group_table <- table(
  chondrocyte_subset@meta.data[[best_res]],
  chondrocyte_subset$group
)

write.csv(
  subcluster_group_table,
  file = file.path(RESULT_DIR, "chondrocyte_subcluster_by_group.csv")
)

# 比例分析
subcluster_prop <- prop.table(subcluster_group_table, margin = 2) * 100

write.csv(
  subcluster_prop,
  file = file.path(RESULT_DIR, "chondrocyte_subcluster_proportion.csv")
)

# 可视化比较
library(reshape2)
prop_df <- melt(subcluster_prop)
colnames(prop_df) <- c("Subcluster", "Group", "Percentage")

p_prop <- ggplot(prop_df, aes(x = Subcluster, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_custom() +
  theme_classic() +
  labs(
    x = "Chondrocyte Subcluster",
    y = "Percentage (%)",
    title = "Chondrocyte Subcluster Composition: SDD vs MDD"
  ) +
  geom_text(
    aes(label = round(Percentage, 1)),
    position = position_dodge(0.9),
    vjust = -0.5,
    size = 3
  )

save_plot(p_prop, "chondrocyte_subcluster_proportion.pdf", width = 10, height = 6)

# 9. 软骨相关基因表达 --------------------------------------------------------
print_progress("步骤8: 软骨相关基因表达")

# 定义软骨相关基因
chondrocyte_genes <- c(
  "SOX9", "SOX5", "SOX6",           # SOX转录因子
  "COL2A1", "COL9A1", "COL11A1",    # 胶原蛋白
  "ACAN", "COMP", "MIA",             # 软骨基质
  "PRG4", "CILP", "FMOD"             # 其他软骨标记
)

# 检查哪些基因存在
available_genes <- chondrocyte_genes[chondrocyte_genes %in% rownames(chondrocyte_subset)]

if(length(available_genes) > 0) {
  # FeaturePlot
  p_feature <- FeaturePlot(
    chondrocyte_subset,
    features = available_genes[1:min(9, length(available_genes))],
    reduction = "umap",
    ncol = 3
  ) &
    scale_color_gradientn(colors = c("lightgrey", "yellow", "red", "darkred"))
  
  save_plot(p_feature, "chondrocyte_marker_genes_featureplot.pdf", 
            width = 15, height = 15)
  
  # VlnPlot
  p_vln <- VlnPlot(
    chondrocyte_subset,
    features = available_genes[1:min(6, length(available_genes))],
    ncol = 3,
    pt.size = 0
  )
  
  save_plot(p_vln, "chondrocyte_marker_genes_violinplot.pdf", 
            width = 15, height = 10)
}

# 10. 保存最终结果 -----------------------------------------------------------
save(chondrocyte_subset, 
     file = file.path(DATA_DIR, "chondrocyte_final.Rdata"))

print_progress("软骨细胞亚群分析完成")

# 生成报告
generate_analysis_report(
  chondrocyte_subset,
  analysis_name = "软骨细胞亚群分析",
  output_file = file.path(RESULT_DIR, "chondrocyte_subanalysis_report.txt")
)

message("\n结果保存位置:")
message(paste("- 数据:", file.path(DATA_DIR, "chondrocyte_final.Rdata")))
message(paste("- 图片:", FIGURE_DIR))
message(paste("- 表格:", RESULT_DIR))

