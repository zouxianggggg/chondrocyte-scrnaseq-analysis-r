# =============================================================================
# 第15亚群拟时序CYTL1表达分析
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")
load_packages()
library(monocle)
library(ggpubr)
library(RColorBrewer)

print_progress("第15亚群拟时序CYTL1表达分析")

# 参数设置 -------------------------------------------------------------------
TARGET_CLUSTER <- "15"
TARGET_GENES <- c("CYTL1", "COL1A1", "COL1A2", "SOX9")

# 1. 加载拟时序数据 ----------------------------------------------------------
message("加载拟时序数据...")
load(file.path(DATA_DIR, "trajectory_initial.RData"))

# 确定分辨率
if("RNA_snn_res.0.5" %in% colnames(pData(cds))) {
  resolution <- "RNA_snn_res.0.5"
} else {
  resolution <- "RNA_snn_res.0.7"
}

message("拟时序数据概况:")
message(paste("细胞数:", ncol(cds)))
message(paste("基因数:", nrow(cds)))
message("State分布:")
print(table(pData(cds)$State))

# 2. 标记第15亚群 ------------------------------------------------------------
message("\n标记第15亚群...")
pData(cds)$is_cluster15 <- ifelse(
  pData(cds)[[resolution]] == TARGET_CLUSTER,
  "Cluster 15",
  "Other Clusters"
)

cluster15_cells <- sum(pData(cds)$is_cluster15 == "Cluster 15")
message(paste("第15亚群细胞数:", cluster15_cells))

# 提取Pseudotime和群组信息
trajectory_df <- data.frame(
  cell = colnames(cds),
  pseudotime = pData(cds)$Pseudotime,
  state = pData(cds)$State,
  cluster = pData(cds)[[resolution]],
  is_cluster15 = pData(cds)$is_cluster15,
  group = pData(cds)$group
)

# 保存详细信息
write.csv(
  trajectory_df,
  file = file.path(RESULT_DIR, "9_第15亚群_拟时序信息.csv"),
  row.names = FALSE
)

# =============================================================================
# 一、第15亚群在拟时序中的位置
# =============================================================================
print_progress("分析一: 第15亚群在拟时序中的位置")

# 1. 轨迹图高亮第15亚群
p_trajectory_cluster15 <- plot_cell_trajectory(cds, color_by = "is_cluster15") +
  scale_color_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey80")) +
  theme(text = element_text(size = 14))

save_plot(p_trajectory_cluster15, "9_第15亚群_拟时序位置.pdf", width = 10, height = 8)

# 2. 按分组分面展示第15亚群
p_trajectory_group <- plot_cell_trajectory(cds, color_by = "is_cluster15") +
  facet_wrap(~group) +
  scale_color_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey80")) +
  theme(text = element_text(size = 14))

save_plot(p_trajectory_group, "9_第15亚群_拟时序分组位置.pdf", width = 15, height = 7)

# 统计Pseudotime（先统计再绘图）
# 检查第15亚群是否存在
if(cluster15_cells == 0) {
  message("\n警告: 第15亚群没有细胞，跳过后续分析")
  quit(save = "no")
}

pseudotime_stats <- trajectory_df %>%
  group_by(is_cluster15) %>%
  summarise(
    n = n(),
    mean_pseudotime = mean(pseudotime, na.rm = TRUE),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    sd_pseudotime = sd(pseudotime, na.rm = TRUE),
    .groups = "drop"
  )

message("\nPseudotime统计（第15亚群 vs 其他）:")
print(pseudotime_stats)

# 检查是否有两组数据
if(nrow(pseudotime_stats) == 2) {
  wilcox_pseudo <- wilcox.test(pseudotime ~ is_cluster15, data = trajectory_df)
  message(paste("Wilcoxon检验 p-value:", format(wilcox_pseudo$p.value, scientific = TRUE)))
} else {
  message("警告: 只有一个组，无法进行统计检验")
  wilcox_pseudo <- list(p.value = NA)
}

if(nrow(pseudotime_stats) > 0) {
  write.csv(
    pseudotime_stats,
    file = file.path(RESULT_DIR, "9_第15亚群_Pseudotime统计.csv"),
    row.names = FALSE
  )
}

# 3. Pseudotime分布比较图（第15亚群 vs 其他）
if(nrow(pseudotime_stats) == 2) {
  p_pseudo_cluster15 <- ggplot(trajectory_df, aes(x = is_cluster15, y = pseudotime, fill = is_cluster15)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey70")) +
    theme_classic() +
    labs(
      x = "",
      y = "Pseudotime"
    ) +
    stat_compare_means(method = "wilcox.test", size = 5) +
    theme(legend.position = "none")
  
  save_plot(p_pseudo_cluster15, "9_第15亚群_Pseudotime分布.pdf", width = 8, height = 6)
} else {
  message("跳过Pseudotime分布图（只有一个组）")
}

# =============================================================================
# 二、CYTL1在第15亚群中的拟时序表达
# =============================================================================
print_progress("分析二: CYTL1在第15亚群的拟时序表达")

# 检查CYTL1是否存在
if("CYTL1" %in% rownames(cds)) {
  message("✓ CYTL1在数据中")
  
  # 1. CYTL1在拟时序轨迹上的表达
  p_cytl1_trajectory <- plot_cell_trajectory(cds, markers = "CYTL1") +
    theme(text = element_text(size = 14))
  
  save_plot(p_cytl1_trajectory, "9_CYTL1_拟时序轨迹表达.pdf", width = 10, height = 8)
  
  # 2. CYTL1按第15亚群分组的拟时序轨迹（使用patchwork组合）
  # 注意：plot_cell_trajectory在子集数据时可能会有问题，改用整体图+高亮
  # 改为在整体轨迹上用不同颜色标记CYTL1表达
  
  # 提取CYTL1表达数据用于后续分析
  cytl1_expr <- exprs(cds)["CYTL1", ]
  
  # 创建一个自定义的CYTL1表达图，只显示第15亚群
  # 获取UMAP坐标
  umap_coords <- data.frame(
    Component_1 = pData(cds)$State,  # 这里需要从reducedDimS提取实际坐标
    Component_2 = pData(cds)$State,
    CYTL1 = cytl1_expr,
    is_cluster15 = pData(cds)$is_cluster15
  )
  
  # 由于plot_cell_trajectory的子集问题，我们跳过分组轨迹图
  message("注意: 跳过子集轨迹图以避免维度不匹配错误")
  
  # 3. CYTL1表达沿Pseudotime的变化（第15亚群高亮）
  # CYTL1表达数据已在前面提取
  trajectory_df$CYTL1 <- cytl1_expr
  
  # 按第15亚群着色的散点图
  p_cytl1_pseudo_scatter <- ggplot(trajectory_df, 
                                   aes(x = pseudotime, y = CYTL1, color = is_cluster15)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.5) +
    scale_color_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey60")) +
    theme_classic() +
    labs(
      x = "Pseudotime",
      y = "CYTL1 Expression",
      color = ""
    ) +
    theme(text = element_text(size = 14))
  
  save_plot(p_cytl1_pseudo_scatter, "9_CYTL1_Pseudotime表达趋势.pdf", width = 10, height = 6)
  
  # 4. 分面展示CYTL1在第15亚群和其他细胞的表达趋势
  p_cytl1_pseudo_facet <- ggplot(trajectory_df, 
                                 aes(x = pseudotime, y = CYTL1)) +
    geom_point(alpha = 0.3, size = 1, color = "grey70") +
    geom_smooth(method = "loess", se = TRUE, color = "#E64B35", fill = "#E64B35", alpha = 0.3) +
    facet_wrap(~is_cluster15, ncol = 2) +
    theme_bw() +
    labs(
      x = "Pseudotime",
      y = "CYTL1 Expression"
    ) +
    theme(text = element_text(size = 14))
  
  save_plot(p_cytl1_pseudo_facet, "9_CYTL1_Pseudotime分面表达.pdf", width = 14, height = 6)
  
  # 5. CYTL1沿Pseudotime的表达统计（binned分析）
  # 将Pseudotime分成10个区间
  trajectory_df$pseudo_bin <- cut(trajectory_df$pseudotime, 
                                  breaks = 10, 
                                  labels = paste0("Bin", 1:10))
  
  # 计算每个bin的CYTL1平均表达
  cytl1_bin_stats <- trajectory_df %>%
    filter(!is.na(pseudo_bin)) %>%
    group_by(pseudo_bin, is_cluster15) %>%
    summarise(
      mean_cytl1 = mean(CYTL1, na.rm = TRUE),
      median_cytl1 = median(CYTL1, na.rm = TRUE),
      sd_cytl1 = sd(CYTL1, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # 折线图展示CYTL1沿Pseudotime的变化
  p_cytl1_line <- ggplot(cytl1_bin_stats, 
                         aes(x = pseudo_bin, y = mean_cytl1, 
                             group = is_cluster15, color = is_cluster15)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_cytl1 - sd_cytl1, ymax = mean_cytl1 + sd_cytl1), 
                  width = 0.2) +
    scale_color_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey60")) +
    theme_classic() +
    labs(
      x = "Pseudotime Bin",
      y = "Mean CYTL1 Expression",
      color = ""
    ) +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot(p_cytl1_line, "9_CYTL1_Pseudotime分段表达.pdf", width = 10, height = 6)
  
  # 保存统计结果
  write.csv(
    cytl1_bin_stats,
    file = file.path(RESULT_DIR, "9_CYTL1_Pseudotime分段统计.csv"),
    row.names = FALSE
  )
  
  # 6. 检验CYTL1在第15亚群是否显著增加
  # 按早期/中期/晚期分组
  trajectory_df$stage <- cut(
    trajectory_df$pseudotime,
    breaks = quantile(trajectory_df$pseudotime, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE),
    labels = c("Early", "Middle", "Late"),
    include.lowest = TRUE
  )
  
  # 计算第15亚群在各阶段的CYTL1表达
  cytl1_stage_stats <- trajectory_df %>%
    filter(is_cluster15 == "Cluster 15" & !is.na(stage)) %>%
    group_by(stage) %>%
    summarise(
      n = n(),
      mean_cytl1 = mean(CYTL1, na.rm = TRUE),
      median_cytl1 = median(CYTL1, na.rm = TRUE),
      sd_cytl1 = sd(CYTL1, na.rm = TRUE),
      .groups = "drop"
    )
  
  message("\nCYTL1在第15亚群各阶段的表达:")
  print(cytl1_stage_stats)
  
  # Kruskal-Wallis检验
  if(nrow(cytl1_stage_stats) > 1) {
    kruskal_test <- kruskal.test(CYTL1 ~ stage, 
                                  data = trajectory_df[trajectory_df$is_cluster15 == "Cluster 15" & !is.na(trajectory_df$stage), ])
    message(paste("\nKruskal-Wallis检验 p-value:", format(kruskal_test$p.value, scientific = TRUE)))
  }
  
  # 可视化各阶段表达
  p_cytl1_stage <- ggplot(
    trajectory_df[trajectory_df$is_cluster15 == "Cluster 15" & !is.na(trajectory_df$stage), ],
    aes(x = stage, y = CYTL1, fill = stage)
  ) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    labs(
      x = "Trajectory Stage",
      y = "CYTL1 Expression"
    ) +
    stat_compare_means(comparisons = list(c("Early", "Middle"), c("Middle", "Late"), c("Early", "Late"))) +
    theme(legend.position = "none", text = element_text(size = 14))
  
  save_plot(p_cytl1_stage, "9_CYTL1_第15亚群分阶段表达.pdf", width = 8, height = 6)
  
  write.csv(
    cytl1_stage_stats,
    file = file.path(RESULT_DIR, "9_CYTL1_第15亚群分阶段统计.csv"),
    row.names = FALSE
  )
  
} else {
  message("✗ CYTL1不在数据中")
}

# =============================================================================
# 三、多基因拟时序动态表达比较（CYTL1, COL1A1, COL1A2, SOX9）
# =============================================================================
print_progress("分析三: 多基因拟时序动态表达")

# 检查目标基因
available_genes <- TARGET_GENES[TARGET_GENES %in% rownames(cds)]
missing_genes <- TARGET_GENES[!TARGET_GENES %in% rownames(cds)]

if(length(missing_genes) > 0) {
  message("警告: 以下基因不在数据中:")
  print(missing_genes)
}

message(paste("分析", length(available_genes), "个基因:", paste(available_genes, collapse = ", ")))

if(length(available_genes) > 0) {
  
  # 1. 使用monocle的plot_genes_in_pseudotime
  cds_genes <- cds[available_genes, ]
  
  # 所有细胞的基因表达趋势
  p_genes_all <- plot_genes_in_pseudotime(
    cds_genes,
    color_by = "group",
    ncol = 2
  ) +
    scale_color_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000"))
  
  save_plot(p_genes_all, "9_多基因_拟时序表达趋势_全部细胞.pdf", width = 12, height = 10)
  
  # 第15亚群的基因表达趋势
  if(cluster15_cells > 0) {
    cds_genes_cluster15 <- cds_genes[, pData(cds_genes)$is_cluster15 == "Cluster 15"]
    
    p_genes_cluster15 <- plot_genes_in_pseudotime(
      cds_genes_cluster15,
      color_by = "group",
      ncol = 2
    ) +
      scale_color_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000"))
    
    save_plot(p_genes_cluster15, "9_多基因_拟时序表达趋势_第15亚群.pdf", width = 12, height = 10)
  }
  
  # 2. 自定义可视化：四个基因的表达趋势对比
  # 提取基因表达
  genes_expr_df <- trajectory_df
  for(gene in available_genes) {
    genes_expr_df[[gene]] <- exprs(cds)[gene, ]
  }
  
  # 转换为长格式
  genes_expr_long <- genes_expr_df %>%
    dplyr::select(cell, pseudotime, is_cluster15, group, all_of(available_genes)) %>%
    pivot_longer(
      cols = all_of(available_genes),
      names_to = "Gene",
      values_to = "Expression"
    )
  
  # 3. 所有基因在一张图中（按第15亚群着色）
  p_genes_combined <- ggplot(genes_expr_long, 
                             aes(x = pseudotime, y = Expression, color = Gene)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, alpha = 0.8, span = 0.75) +
    facet_wrap(~is_cluster15, ncol = 2) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      x = "Pseudotime",
      y = "Expression"
    ) +
    theme(text = element_text(size = 14))
  
  save_plot(p_genes_combined, "9_多基因_拟时序对比_第15亚群分面.pdf", width = 14, height = 6)
  
  # 4. 每个基因单独一行，按第15亚群分列
  p_genes_grid <- ggplot(genes_expr_long, 
                         aes(x = pseudotime, y = Expression)) +
    geom_point(aes(color = is_cluster15), alpha = 0.2, size = 0.8) +
    geom_smooth(aes(color = is_cluster15, fill = is_cluster15), 
                method = "loess", se = FALSE, linewidth = 1.5, span = 0.75) +
    facet_wrap(~Gene, ncol = 2, scales = "free_y") +
    scale_color_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey60")) +
    scale_fill_manual(values = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey60")) +
    theme_bw() +
    labs(
      x = "Pseudotime",
      y = "Expression",
      color = "",
      fill = ""
    ) +
    theme(text = element_text(size = 14))
  
  save_plot(p_genes_grid, "9_多基因_各基因拟时序趋势.pdf", width = 12, height = 10)
  
  # 5. 计算基因表达相关性（沿Pseudotime）
  genes_cor_matrix <- cor(genes_expr_df[, available_genes], 
                           use = "complete.obs", 
                           method = "spearman")
  
  message("\n基因表达相关性（Spearman）:")
  print(round(genes_cor_matrix, 3))
  
  # 保存相关性矩阵
  write.csv(
    genes_cor_matrix,
    file = file.path(RESULT_DIR, "9_多基因_表达相关性矩阵.csv")
  )
  
  # 相关性热图
  library(pheatmap)
  pdf(file.path(FIGURE_DIR, "9_多基因_表达相关性热图.pdf"), width = 8, height = 7)
  pheatmap(
    genes_cor_matrix,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "",
    fontsize = 12
  )
  dev.off()
  
  # 6. 分段统计各基因的表达模式
  genes_bin_stats <- genes_expr_long %>%
    mutate(pseudo_bin = cut(pseudotime, breaks = 10, labels = paste0("Bin", 1:10))) %>%
    filter(!is.na(pseudo_bin)) %>%
    group_by(Gene, pseudo_bin, is_cluster15) %>%
    summarise(
      mean_expr = mean(Expression, na.rm = TRUE),
      median_expr = median(Expression, na.rm = TRUE),
      sd_expr = sd(Expression, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # 折线图对比
  p_genes_line <- ggplot(genes_bin_stats, 
                         aes(x = pseudo_bin, y = mean_expr, 
                             group = interaction(Gene, is_cluster15),
                             color = Gene, linetype = is_cluster15)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_brewer(palette = "Set1") +
    scale_linetype_manual(values = c("Cluster 15" = "solid", "Other Clusters" = "dashed")) +
    theme_classic() +
    labs(
      x = "Pseudotime Bin",
      y = "Mean Expression",
      color = "Gene",
      linetype = ""
    ) +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot(p_genes_line, "9_多基因_Pseudotime分段折线图.pdf", width = 12, height = 7)
  
  # 保存分段统计
  write.csv(
    genes_bin_stats,
    file = file.path(RESULT_DIR, "9_多基因_Pseudotime分段统计.csv"),
    row.names = FALSE
  )
  
  # 7. 热图展示基因沿Pseudotime的表达模式
  # 选择部分细胞避免过大
  set.seed(123)
  sampled_cells <- sample(colnames(cds), min(1000, ncol(cds)))
  
  genes_matrix <- exprs(cds)[available_genes, sampled_cells]
  cell_order <- order(pData(cds)[sampled_cells, "Pseudotime"])
  genes_matrix_ordered <- genes_matrix[, cell_order]
  
  # 注释
  col_anno <- data.frame(
    Pseudotime = pData(cds)[sampled_cells[cell_order], "Pseudotime"],
    Cluster = pData(cds)[sampled_cells[cell_order], "is_cluster15"],
    Group = pData(cds)[sampled_cells[cell_order], "group"]
  )
  rownames(col_anno) <- colnames(genes_matrix_ordered)
  
  anno_colors <- list(
    Cluster = c("Cluster 15" = "#E64B35", "Other Clusters" = "grey70"),
    Group = c("MDD" = "#00A087", "SDD" = "#DC0000"),
    Pseudotime = colorRampPalette(c("blue", "cyan", "yellow", "red"))(100)
  )
  
  pdf(file.path(FIGURE_DIR, "9_多基因_Pseudotime热图.pdf"), width = 12, height = 6)
  pheatmap(
    genes_matrix_ordered,
    annotation_col = col_anno,
    annotation_colors = anno_colors,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "",
    fontsize_row = 12
  )
  dev.off()
  
}

# =============================================================================
# 总结报告
# =============================================================================
print_progress("生成综合分析报告")

summary_text <- c(
  "=============================================================================",
  "第15亚群拟时序CYTL1表达分析报告",
  "=============================================================================",
  "",
  paste("分析日期:", Sys.time()),
  paste("总细胞数:", ncol(cds)),
  paste("第15亚群细胞数:", cluster15_cells),
  "",
  "--- 拟时序分布 ---",
  paste("第15亚群平均Pseudotime:", 
        round(mean(trajectory_df$pseudotime[trajectory_df$is_cluster15 == "Cluster 15"], na.rm = TRUE), 3)),
  paste("其他亚群平均Pseudotime:", 
        round(mean(trajectory_df$pseudotime[trajectory_df$is_cluster15 == "Other Clusters"], na.rm = TRUE), 3)),
  ifelse(!is.na(wilcox_pseudo$p.value),
         paste("Wilcoxon检验 p-value:", format(wilcox_pseudo$p.value, scientific = TRUE)),
         "Wilcoxon检验: 无法进行（数据不足）"),
  ""
)

if("CYTL1" %in% available_genes) {
  cytl1_cluster15_mean <- mean(trajectory_df$CYTL1[trajectory_df$is_cluster15 == "Cluster 15"], na.rm = TRUE)
  cytl1_others_mean <- mean(trajectory_df$CYTL1[trajectory_df$is_cluster15 == "Other Clusters"], na.rm = TRUE)
  
  summary_text <- c(
    summary_text,
    "--- CYTL1表达 ---",
    paste("第15亚群平均表达:", round(cytl1_cluster15_mean, 3)),
    paste("其他亚群平均表达:", round(cytl1_others_mean, 3)),
    paste("表达倍数变化:", round(cytl1_cluster15_mean / cytl1_others_mean, 2), "倍"),
    ""
  )
}

summary_text <- c(
  summary_text,
  "--- 分析基因 ---",
  paste("可用基因:", paste(available_genes, collapse = ", ")),
  "",
  "=============================================================================",
  "结论:",
  "1. 第15亚群在拟时序中的分布特征已可视化",
  "2. CYTL1在第15亚群沿拟时序的动态表达已分析",
  "3. CYTL1与COL1A1、COL1A2、SOX9的表达趋势已对比",
  "4. 详细结果请查看各分析步骤生成的图表和表格",
  "============================================================================="
)

writeLines(summary_text, file.path(RESULT_DIR, "9_第15亚群拟时序分析总结.txt"))

print_progress("第15亚群拟时序分析完成！")

message("\n生成的结果文件:")
message(paste("- 数据表格:", RESULT_DIR))
message(paste("- 可视化图表:", FIGURE_DIR))
message(paste("- 综合报告:", file.path(RESULT_DIR, "9_第15亚群拟时序分析总结.txt")))
