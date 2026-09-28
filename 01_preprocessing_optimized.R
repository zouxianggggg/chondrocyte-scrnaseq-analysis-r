# =============================================================================
# 单细胞测序分析主流程 - 优化版
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")

# 加载包
load_packages()
library(harmony)
library(glmGamPoi)
library(clustree)

print_progress("步骤1: 数据加载与质控")

# 1. 创建Seurat对象 ---------------------------------------------------------
seurat_list <- list()
for(i in 1:nrow(SAMPLE_INFO)) {
  sample_id <- SAMPLE_INFO$sample_id[i]
  data_path <- SAMPLE_INFO$data_path[i]
  
  message(paste("加载样本:", sample_id))
  seurat_list[[sample_id]] <- CreateSeuratObject(
    counts = Read10X(data.dir = data_path),
    project = sample_id
  )
  
  # 计算QC指标
  seurat_list[[sample_id]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_list[[sample_id]], pattern = "^MT-"
  )
  seurat_list[[sample_id]][["percent.rb"]] <- PercentageFeatureSet(
    seurat_list[[sample_id]], 
    pattern = "^(HBB|HBA1|HBA2|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ)$"
  )
}

# 2. 合并数据 ----------------------------------------------------------------
combined <- merge(
  seurat_list[[1]], 
  y = seurat_list[-1],
  add.cell.ids = SAMPLE_INFO$sample_id,
  project = "Combined"
)

# 添加分组信息
combined$group <- SAMPLE_INFO$group[match(combined$orig.ident, SAMPLE_INFO$sample_id)]

message(paste("合并后总细胞数:", ncol(combined)))

# 3. 质控可视化 --------------------------------------------------------------
print_progress("生成质控图")

# 散点图
plot_qc_scatter(combined, save_dir = FIGURE_DIR)

# 小提琴图
p_qc_before <- plot_qc_metrics(combined)
save_plot(p_qc_before, "qc_before_filtering.pdf", width = 30, height = 10)

# 4. 质控过滤 ----------------------------------------------------------------
print_progress("步骤2: 质控过滤")

message("过滤前细胞数:", ncol(combined))
combined <- subset(
  combined, 
  subset = nFeature_RNA > QC_PARAMS$min_features & 
    percent.mt < QC_PARAMS$max_mt_percent & 
    percent.rb < QC_PARAMS$max_rb_percent
)
message("过滤后细胞数:", ncol(combined))

# 质控后可视化
p_qc_after <- plot_qc_metrics(combined)
save_plot(p_qc_after, "qc_after_filtering.pdf", width = 30, height = 10)

# 5. 标准化与降维 ------------------------------------------------------------
print_progress("步骤3: 数据标准化与降维")

combined <- standard_preprocess(
  combined, 
  n_features = 2000,
  vars_to_regress = "percent.mt",
  n_pcs = 50
)

# 查看高变基因
top10_var_genes <- head(VariableFeatures(combined), 10)
message("Top 10高变基因:")
print(top10_var_genes)

# 6. Harmony批次校正 ---------------------------------------------------------
print_progress("步骤4: Harmony批次校正")

combined_harmony <- run_harmony_integration(
  combined, 
  group_by = "orig.ident",
  n_pcs = CLUSTER_PARAMS$n_pcs
)

# 7. 聚类分析 ----------------------------------------------------------------
print_progress("步骤5: 聚类分析")

combined_harmony <- run_clustering(
  combined_harmony,
  reduction = "harmony",
  dims = 1:CLUSTER_PARAMS$n_pcs,
  resolutions = CLUSTER_PARAMS$resolution_range,
  run_umap = TRUE,
  run_tsne = TRUE
)

# 保存中间结果
save(combined_harmony, file = file.path(DATA_DIR, "scRNA_harmony.Rdata"))

# 8. 聚类可视化 --------------------------------------------------------------
print_progress("步骤6: 聚类可视化")

# 聚类树
p_clustree <- clustree(combined_harmony)
save_plot(p_clustree, "clustree.pdf", width = 12, height = 12)

# 选择最佳分辨率
best_resolution <- paste0("RNA_snn_res.", CLUSTER_PARAMS$resolution)
Idents(combined_harmony) <- best_resolution

# 绘制聚类总览
p_overview <- plot_clustering_overview(
  combined_harmony,
  group_by_vars = c("orig.ident", "group"),
  reduction = "umap",
  cluster_col = best_resolution
)
save_plot(p_overview, "clustering_overview.pdf", width = 15, height = 15)

# 分组比较
p_split <- DimPlot(combined_harmony, reduction = "umap", split.by = "group")
save_plot(p_split, "umap_split_by_group.pdf", width = 15, height = 8)

# 保存最终结果
save(combined_harmony, file = file.path(DATA_DIR, "scRNA_harmony_final.Rdata"))

print_progress("预处理完成！")

# 生成报告
generate_analysis_report(
  combined_harmony, 
  analysis_name = "软骨细胞预处理",
  output_file = file.path(RESULT_DIR, "preprocessing_report.txt")
)

message("\n所有结果已保存到:")
message(paste("- 数据:", DATA_DIR))
message(paste("- 图片:", FIGURE_DIR))
