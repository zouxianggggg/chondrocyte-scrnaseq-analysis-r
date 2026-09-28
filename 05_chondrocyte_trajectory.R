# =============================================================================
# 软骨细胞拟时序分析 - 优化版
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")

# 加载必需的包
library(Seurat)
library(tidyverse)
library(monocle)
library(ggsci)
library(ggpubr)

print_progress("软骨细胞拟时序分析")

# 设置参数 -------------------------------------------------------------------
N_CORES <- 16              # 并行核心数
TOP_N_GENES <- 1000        # 用于拟时序的基因数
MIN_EXPR <- 0.1            # 最小表达阈值
MIN_CELLS <- 10            # 最小细胞数
TEXT_SIZE <- 14            # 文本大小
START_CLUSTER <- "0"       # 起始簇（根据实际调整）

# 设置Java内存
options(java.parameters = "-Xmx20g")

# 1. 加载数据 ----------------------------------------------------------------
message("加载软骨细胞数据...")
load(file.path(DATA_DIR, "chondrocyte_final.Rdata"))

sce <- chondrocyte_subset

# 确定使用的分辨率
if("RNA_snn_res.0.5" %in% colnames(sce@meta.data)) {
  resolution <- "RNA_snn_res.0.5"
} else {
  resolution <- "RNA_snn_res.0.7"
}

Idents(sce) <- resolution

message("数据概况:")
message(paste("细胞数:", ncol(sce)))
message(paste("基因数:", nrow(sce)))
message("聚类分布:")
print(table(sce@meta.data[[resolution]]))
message("分组分布:")
print(table(sce$group))

# 2. 构建Monocle CellDataSet对象 ---------------------------------------------
print_progress("步骤1: 构建CellDataSet对象")

# 提取counts矩阵
expr_matrix <- GetAssayData(sce, slot = "counts", assay = "RNA")

# 准备样本信息
sample_sheet <- sce@meta.data

# 准备基因信息
gene_annotation <- data.frame(
  gene_short_name = rownames(sce),
  row.names = rownames(sce)
)

# 创建AnnotatedDataFrame
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

# 创建CellDataSet
cds <- newCellDataSet(
  cellData = expr_matrix,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

message("CellDataSet对象创建完成")
print(cds)

# 3. 数据归一化 --------------------------------------------------------------
print_progress("步骤2: 数据归一化")

message("估计size factors...")
cds <- estimateSizeFactors(cds)

message("估计离散度...")
cds <- estimateDispersions(cds)

# 4. 特征基因选择 ------------------------------------------------------------
print_progress("步骤3: 选择用于拟时序的基因")

# 检测表达基因
message("检测表达基因...")
cds <- detectGenes(cds, min_expr = MIN_EXPR)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= MIN_CELLS))
message(paste("检测到", length(expressed_genes), "个表达基因"))

# 差异基因检测（基于聚类）
message("进行差异基因检测（这可能需要一些时间）...")
diff_test_res <- differentialGeneTest(
  cds[expressed_genes, ],
  fullModelFormulaStr = as.formula(paste0("~", resolution)),
  cores = N_CORES
)

# 排序并选择top基因
diff_test_res <- diff_test_res[order(diff_test_res$pval), ]
ordering_genes <- row.names(diff_test_res)[1:min(TOP_N_GENES, nrow(diff_test_res))]

message(paste("选择了", length(ordering_genes), "个基因用于拟时序分析"))

# 设置ordering基因
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)

# 保存ordering基因
write.csv(
  diff_test_res[ordering_genes, ],
  file = file.path(RESULT_DIR, "trajectory_ordering_genes.csv")
)

# 可视化ordering基因
pdf(file.path(FIGURE_DIR, "trajectory_ordering_genes.pdf"), width = 10, height = 8)
plot_ordering_genes(cds)
dev.off()

# 5. 降维和轨迹构建 ----------------------------------------------------------
print_progress("步骤4: 降维和轨迹构建")

message("使用DDRTree进行降维...")
cds <- reduceDimension(
  cds,
  max_components = 2,
  method = 'DDRTree',
  verbose = TRUE
)

message("排序细胞...")
cds <- orderCells(cds)

# 保存初步结果
save(cds, sce, file = file.path(DATA_DIR, "trajectory_initial.RData"))

# 6. 可视化轨迹（初步） ------------------------------------------------------
print_progress("步骤5: 初步轨迹可视化")

theme_trajectory <- theme(text = element_text(size = TEXT_SIZE))

# 按State着色
p_state <- plot_cell_trajectory(cds, color_by = "State") +
  theme_trajectory +
  scale_color_custom() +
  ggtitle("Cell Trajectory by State")
save_plot(p_state, "trajectory_by_state.pdf")

# 按聚类着色
p_cluster <- plot_cell_trajectory(cds, color_by = resolution) +
  theme_trajectory +
  scale_color_custom() +
  ggtitle("Cell Trajectory by Cluster")
save_plot(p_cluster, "trajectory_by_cluster.pdf")

# 按分组着色
p_group <- plot_cell_trajectory(cds, color_by = "group") +
  theme_trajectory +
  scale_color_custom() +
  ggtitle("Cell Trajectory by Group (SDD vs MDD)")
save_plot(p_group, "trajectory_by_group.pdf")

# 按Pseudotime着色
p_pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  theme_trajectory +
  scale_color_gradientn(colors = c("blue", "cyan", "yellow", "red")) +
  ggtitle("Cell Trajectory by Pseudotime")
save_plot(p_pseudotime, "trajectory_by_pseudotime.pdf")

# 7. 设置轨迹起点 ------------------------------------------------------------
print_progress("步骤6: 设置轨迹起点")

# 定义起点选择函数
GM_state <- function(cds, start_cluster = "0") {
  if(length(unique(pData(cds)$State)) > 1) {
    # 检查起始簇是否存在
    state_cluster_table <- table(pData(cds)$State, pData(cds)[[resolution]])
    
    if(!start_cluster %in% colnames(state_cluster_table)) {
      warning(paste0("警告: 簇 '", start_cluster, "' 不存在"))
      message("可用的簇: ", paste(colnames(state_cluster_table), collapse = ", "))
      return(1)
    }
    
    T0_counts <- state_cluster_table[, start_cluster]
    root_state <- as.numeric(names(T0_counts)[which.max(T0_counts)])
    message(paste("选择State", root_state, "作为起点"))
    return(root_state)
  } else {
    return(1)
  }
}

# 重新排序细胞（指定起点）
message(paste("使用聚类", START_CLUSTER, "作为起始点"))
root_state <- GM_state(cds, start_cluster = START_CLUSTER)
cds <- orderCells(cds, root_state = root_state)

# 8. 最终轨迹可视化 ----------------------------------------------------------
print_progress("步骤7: 最终轨迹可视化")

# 主轨迹图（按Pseudotime）
p_final <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  theme_trajectory +
  scale_color_gradientn(colors = c("blue", "cyan", "green", "yellow", "orange", "red")) +
  ggtitle("Chondrocyte Differentiation Trajectory") +
  theme(legend.position = "right")
save_plot(p_final, "trajectory_final_pseudotime.pdf", width = 10, height = 8)

# 按State分面
p_state_facet <- plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap(~State, nrow = 1) +
  theme_trajectory +
  scale_color_custom()
save_plot(p_state_facet, "trajectory_facet_by_state.pdf", width = 15, height = 6)

# 按聚类分面
p_cluster_facet <- plot_cell_trajectory(cds, color_by = resolution) +
  facet_wrap(as.formula(paste0("~", resolution)), nrow = 2) +
  theme_trajectory +
  scale_color_custom()
save_plot(p_cluster_facet, "trajectory_facet_by_cluster.pdf", width = 20, height = 12)

# 按分组分面
p_group_facet <- plot_cell_trajectory(cds, color_by = "group") +
  facet_wrap(~group, nrow = 1) +
  theme_trajectory +
  scale_color_custom() +
  ggtitle("Trajectory Comparison: SDD vs MDD")
save_plot(p_group_facet, "trajectory_facet_by_group.pdf", width = 15, height = 6)

# 组合图
p_combined <- (p_state + p_cluster) / (p_group + p_pseudotime)
save_plot(p_combined, "trajectory_overview.pdf", width = 15, height = 15)

# 9. Pseudotime分布分析 ------------------------------------------------------
print_progress("步骤8: Pseudotime分布分析")

# 提取Pseudotime数据
pseudotime_df <- data.frame(
  cell = colnames(cds),
  pseudotime = pData(cds)$Pseudotime,
  state = pData(cds)$State,
  cluster = pData(cds)[[resolution]],
  group = pData(cds)$group
)

# 保存Pseudotime数据
write.csv(
  pseudotime_df,
  file = file.path(RESULT_DIR, "pseudotime_values.csv"),
  row.names = FALSE
)

# 按聚类比较Pseudotime分布
p_pseudo_cluster <- ggplot(pseudotime_df, aes(x = cluster, y = pseudotime, fill = cluster)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_custom() +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Pseudotime",
    title = "Pseudotime Distribution by Cluster"
  ) +
  theme(legend.position = "none")
save_plot(p_pseudo_cluster, "pseudotime_by_cluster.pdf", width = 10, height = 6)

# 按分组比较Pseudotime分布
p_pseudo_group <- ggplot(pseudotime_df, aes(x = group, y = pseudotime, fill = group)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_custom() +
  theme_classic() +
  labs(
    x = "Group",
    y = "Pseudotime",
    title = "Pseudotime Distribution: SDD vs MDD"
  ) +
  stat_compare_means(method = "wilcox.test")
save_plot(p_pseudo_group, "pseudotime_by_group.pdf", width = 8, height = 6)

# 密度图
p_density <- ggplot(pseudotime_df, aes(x = pseudotime, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_custom() +
  theme_classic() +
  labs(
    x = "Pseudotime",
    y = "Density",
    title = "Pseudotime Distribution Density"
  )
save_plot(p_density, "pseudotime_density.pdf", width = 10, height = 6)

# 10. 拟时序相关基因分析 -----------------------------------------------------
print_progress("步骤9: 拟时序相关基因分析")

# 软骨相关基因
chondro_genes <- c(
  "SOX9", "SOX5", "SOX6",
  "COL2A1", "COL1A1", "ACAN",
  "COMP", "MIA", "PRG4"
)

# 检查存在的基因
available_genes <- chondro_genes[chondro_genes %in% rownames(cds)]

if(length(available_genes) > 0) {
  message(paste("分析", length(available_genes), "个软骨相关基因"))
  
  # 绘制基因表达趋势
  cds_subset <- cds[available_genes, ]
  
  p_genes <- plot_genes_in_pseudotime(
    cds_subset,
    color_by = "group",
    ncol = 3
  ) +
    theme(text = element_text(size = 12)) +
    scale_color_custom()
  
  save_plot(p_genes, "trajectory_chondrocyte_genes.pdf", width = 15, height = 10)
}

# 11. 差异表达分析（基于Pseudotime） ----------------------------------------
print_progress("步骤10: 基于Pseudotime的差异分析")

message("分析基因表达随拟时序的变化...")
diff_test_pseudotime <- differentialGeneTest(
  cds[expressed_genes[1:min(1000, length(expressed_genes))], ],
  fullModelFormulaStr = "~sm.ns(Pseudotime)",
  cores = N_CORES
)

# 排序
diff_test_pseudotime <- diff_test_pseudotime[order(diff_test_pseudotime$qval), ]

# 保存结果
write.csv(
  diff_test_pseudotime,
  file = file.path(RESULT_DIR, "trajectory_pseudotime_DEGs.csv")
)

message(paste("发现", sum(diff_test_pseudotime$qval < 0.05, na.rm = TRUE), 
              "个显著的拟时序相关基因"))

# 提取top拟时序相关基因
top_pseudotime_genes <- head(row.names(diff_test_pseudotime), 50)

# 绘制热图
if(length(top_pseudotime_genes) > 0) {
  p_heatmap <- plot_pseudotime_heatmap(
    cds[top_pseudotime_genes, ],
    num_clusters = 4,
    cores = N_CORES,
    show_rownames = TRUE
  )
  
  pdf(file.path(FIGURE_DIR, "trajectory_pseudotime_heatmap.pdf"), 
      width = 10, height = 12)
  print(p_heatmap)
  dev.off()
}

# 12. 保存最终结果 -----------------------------------------------------------
save(cds, sce, pseudotime_df, diff_test_pseudotime,
     file = file.path(DATA_DIR, "trajectory_final.RData"))

print_progress("拟时序分析完成！")

# 生成摘要
summary_text <- c(
  "=== 拟时序分析摘要 ===",
  paste("细胞数:", ncol(cds)),
  paste("使用基因数:", length(ordering_genes)),
  paste("State数:", length(unique(pData(cds)$State))),
  paste("起始State:", root_state),
  paste("拟时序相关基因数:", sum(diff_test_pseudotime$qval < 0.05, na.rm = TRUE)),
  "",
  "Pseudotime统计 (by Group):",
  capture.output(
    tapply(pseudotime_df$pseudotime, pseudotime_df$group, summary)
  )
)

writeLines(summary_text, file.path(RESULT_DIR, "trajectory_summary.txt"))

message("\n所有结果已保存于:")
message(paste("- 数据:", file.path(DATA_DIR, "trajectory_final.RData")))
message(paste("- 图片:", FIGURE_DIR))
message(paste("- 表格:", RESULT_DIR))

