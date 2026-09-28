# =============================================================================
# 单细胞测序细胞类型注释 - 优化版
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")
load_packages()

print_progress("细胞类型注释分析")

# 1. 加载数据 ----------------------------------------------------------------
load(file.path(DATA_DIR, "scRNA_harmony_final.Rdata"))
sce <- combined_harmony

# 设置分辨率
resolution <- paste0("RNA_snn_res.", CLUSTER_PARAMS$resolution)
Idents(sce) <- resolution

# 查看聚类分布
message("聚类分布:")
print(table(sce@meta.data[[resolution]]))

# 2. 标记基因可视化 ----------------------------------------------------------
print_progress("步骤1: 标记基因表达分析")

# DotPlot展示所有标记基因
p_dotplot <- DotPlot(
  sce, 
  features = CELL_TYPE_MARKERS,
  assay = "RNA",
  group.by = resolution
) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Type Marker Genes Expression")

save_plot(p_dotplot, "marker_genes_dotplot.pdf", width = 30, height = 10)

# 3. 细胞类型注释 ------------------------------------------------------------
print_progress("步骤2: 细胞类型注释")

# 创建注释表格（根据marker基因表达和聚类结果）
celltype_annotation <- data.frame(
  ClusterID = 0:15,
  celltype = 'unknown',
  stringsAsFactors = FALSE
)

# 基于你的原始脚本的注释结果
celltype_annotation$celltype[celltype_annotation$ClusterID %in% c(0,1,2,3,4,5,6,11,12)] <- "Chondrocyte"
celltype_annotation$celltype[celltype_annotation$ClusterID %in% c(7,15)] <- "AF cell"
celltype_annotation$celltype[celltype_annotation$ClusterID == 8] <- "Endothelia"
celltype_annotation$celltype[celltype_annotation$ClusterID == 9] <- "Erythrocyte"
celltype_annotation$celltype[celltype_annotation$ClusterID == 10] <- "smooth muscle cells"
celltype_annotation$celltype[celltype_annotation$ClusterID == 13] <- "T cell"
celltype_annotation$celltype[celltype_annotation$ClusterID == 14] <- "NP Cell"

# 应用注释
sce <- annotate_cell_types(sce, resolution, celltype_annotation)

# 保存注释表格
write.csv(
  celltype_annotation,
  file = file.path(RESULT_DIR, "celltype_annotation.csv"),
  row.names = FALSE
)

# 4. 注释结果可视化 ----------------------------------------------------------
print_progress("步骤3: 注释结果可视化")

# UMAP图
p_umap_celltype <- DimPlot(
  sce, 
  reduction = "umap", 
  group.by = "celltype", 
  label = TRUE,
  label.size = 5
) + 
  scale_color_npg() +
  ggtitle("Cell Type Annotation (UMAP)") +
  theme(
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.key = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.spacing.y = unit(0.8, "cm"),
    legend.margin = margin(l = 5, r = 5, b = 5, t = 5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold")
  )

ggsave(
  filename = file.path(FIGURE_DIR, "celltype_annotation_umap.pdf"),
  plot = p_umap_celltype,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "transparent"
)

# tSNE图
p_tsne_celltype <- DimPlot(
  sce, 
  reduction = "tsne", 
  group.by = "celltype", 
  label = TRUE,
  label.size = 5
) + 
  scale_color_npg() +
  ggtitle("Cell Type Annotation (tSNE)")

save_plot(p_tsne_celltype, "celltype_annotation_tsne.pdf")

# 按组分面
p_celltype_split <- DimPlot(
  sce, 
  reduction = "umap", 
  group.by = "celltype",
  split.by = "group",
  label = TRUE
) + 
  scale_color_npg()

save_plot(p_celltype_split, "celltype_by_group.pdf", width = 15, height = 8)

# 5. 统计分析 ----------------------------------------------------------------
print_progress("步骤4: 细胞类型统计")

# 细胞类型分布
celltype_count <- table(sce$celltype)
celltype_prop <- prop.table(celltype_count) * 100

stat_df <- data.frame(
  CellType = names(celltype_count),
  Count = as.numeric(celltype_count),
  Percentage = round(as.numeric(celltype_prop), 2)
)

write.csv(stat_df, file = file.path(RESULT_DIR, "celltype_statistics.csv"), 
          row.names = FALSE)

# 按组统计
celltype_by_group <- table(sce$celltype, sce$group)
write.csv(celltype_by_group, 
          file = file.path(RESULT_DIR, "celltype_by_group.csv"))

# 可视化细胞类型分布
p_bar <- ggplot(stat_df, aes(x = reorder(CellType, -Count), y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type", y = "Cell Count", title = "Cell Type Distribution") +
  geom_text(aes(label = Count), vjust = -0.5)

save_plot(p_bar, "celltype_barplot.pdf", width = 10, height = 8)

# 6. 各细胞类型的特征基因 ----------------------------------------------------
print_progress("步骤5: 查找各细胞类型的标记基因")

Idents(sce) <- "celltype"

# 为每种细胞类型找标记基因（排除unknown）
celltype_markers_list <- list()
for(ct in unique(sce$celltype)) {
  if(ct != "unknown") {
    message(paste("分析细胞类型:", ct))
    
    markers <- FindMarkers(
      sce,
      ident.1 = ct,
      ident.2 = NULL,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    
    markers$gene <- rownames(markers)
    markers$celltype <- ct
    celltype_markers_list[[ct]] <- markers
    
    # 保存每种细胞类型的标记基因
    write.csv(
      markers,
      file = file.path(RESULT_DIR, paste0("markers_", ct, ".csv")),
      row.names = FALSE
    )
  }
}

# 合并所有标记基因
all_markers <- do.call(rbind, celltype_markers_list)
write.csv(all_markers, 
          file = file.path(RESULT_DIR, "all_celltype_markers.csv"),
          row.names = FALSE)

# 7. 提取软骨细胞亚群 --------------------------------------------------------
print_progress("步骤6: 提取软骨细胞亚群")

# 筛选Chondrocyte
chondrocyte_cells <- WhichCells(sce, idents = "Chondrocyte")
chondrocyte_subset <- subset(sce, cells = chondrocyte_cells)

message(paste("软骨细胞数量:", ncol(chondrocyte_subset)))
message("分组分布:")
print(table(chondrocyte_subset$group))

# 保存软骨细胞子集
save(chondrocyte_subset, file = file.path(DATA_DIR, "chondrocyte_cells.Rdata"))

# 可视化软骨细胞
p_chondro <- DimPlot(chondrocyte_subset, reduction = "umap", group.by = "group") +
  scale_color_npg() +
  ggtitle("Chondrocyte Subset")

save_plot(p_chondro, "chondrocyte_subset_umap.pdf")

# 8. 保存注释后的完整数据 ----------------------------------------------------
save(sce, file = file.path(DATA_DIR, "scRNA_annotated.Rdata"))

print_progress("细胞类型注释完成！")

# 生成报告
generate_analysis_report(
  sce,
  analysis_name = "细胞类型注释",
  output_file = file.path(RESULT_DIR, "annotation_report.txt")
)

message("\n细胞类型注释完成！")
message("结果保存位置:")
message(paste("- 完整数据:", file.path(DATA_DIR, "scRNA_annotated.Rdata")))
message(paste("- 软骨细胞:", file.path(DATA_DIR, "chondrocyte_cells.Rdata")))
message(paste("- 统计结果:", RESULT_DIR))
