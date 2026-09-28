# =============================================================================
# 软骨细胞各亚群在MDD和SDD组中的分布分析
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")
load_packages()
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(scales)

print_progress("软骨细胞亚群分组分布分析")

# 1. 加载数据 ----------------------------------------------------------------
message("加载软骨细胞数据...")
load(file.path(DATA_DIR, "chondrocyte_clustered.Rdata"))

# 确定分辨率
if("RNA_snn_res.0.5" %in% colnames(chondrocyte_subset@meta.data)) {
  resolution <- "RNA_snn_res.0.5"
} else {
  resolution <- "RNA_snn_res.0.7"
}

message("数据概况:")
message(paste("总细胞数:", ncol(chondrocyte_subset)))
message(paste("使用分辨率:", resolution))
message(paste("聚类数量:", length(unique(chondrocyte_subset@meta.data[[resolution]]))))

# 2. 统计每个cluster在各组中的分布 ------------------------------------------
print_progress("统计分析")

# 创建统计表
cluster_group_table <- table(chondrocyte_subset@meta.data[[resolution]], 
                              chondrocyte_subset$group)
cluster_group_df <- as.data.frame(cluster_group_table)
colnames(cluster_group_df) <- c("Cluster", "Group", "Count")

# 计算每个cluster的总数
cluster_totals <- cluster_group_df %>%
  group_by(Cluster) %>%
  summarise(Total = sum(Count), .groups = "drop")

# 计算百分比
cluster_group_df <- cluster_group_df %>%
  left_join(cluster_totals, by = "Cluster") %>%
  mutate(Percentage = (Count / Total) * 100)

# 计算每组的总细胞数
group_totals <- cluster_group_df %>%
  group_by(Group) %>%
  summarise(Group_Total = sum(Count), .groups = "drop")

cluster_group_df <- cluster_group_df %>%
  left_join(group_totals, by = "Group") %>%
  mutate(Group_Percentage = (Count / Group_Total) * 100)

# 保存详细统计
write.csv(cluster_group_df, 
          file = file.path(RESULT_DIR, "cluster_group_distribution_detailed.csv"),
          row.names = FALSE)

message("\n各亚群在组内的分布比例:")
print(cluster_group_df)

# 3. 可视化分析 -------------------------------------------------------------
print_progress("生成可视化图表")

# 准备颜色
n_clusters <- length(unique(cluster_group_df$Cluster))
if(n_clusters <= 12) {
  cluster_colors <- brewer.pal(max(3, n_clusters), "Set3")
} else {
  base_colors <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
  cluster_colors <- colorRampPalette(base_colors)(n_clusters)
}

group_colors <- c("MDD" = "#00A087", "SDD" = "#DC0000")

# (1) 堆叠柱状图 - 每个cluster的组成
p1 <- ggplot(cluster_group_df, aes(x = Cluster, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Cell Count",
    title = "Cell Distribution: Each Cluster by Group",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top"
  ) +
  geom_text(aes(label = Count), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "white")

save_plot(p1, "cluster_distribution_stacked_count.pdf", width = 12, height = 8)

# (2) 百分比堆叠柱状图 - 每个cluster的组成比例
p2 <- ggplot(cluster_group_df, aes(x = Cluster, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Percentage (%)",
    title = "Cell Proportion: Each Cluster by Group",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top"
  ) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "white")

save_plot(p2, "cluster_distribution_stacked_percentage.pdf", width = 12, height = 8)

# (3) 分组柱状图 - 对比各cluster在两组中的数量
p3 <- ggplot(cluster_group_df, aes(x = Cluster, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Cell Count",
    title = "Cell Count Comparison: MDD vs SDD by Cluster",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top"
  ) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3)

save_plot(p3, "cluster_distribution_dodged.pdf", width = 14, height = 8)

# (4) 按组分面的柱状图
p4 <- ggplot(cluster_group_df, aes(x = Cluster, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Group, scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Cell Count",
    title = "Cluster Distribution in Each Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3)

save_plot(p4, "cluster_distribution_faceted.pdf", width = 14, height = 8)

# (5) 热图 - cluster x group
cluster_group_matrix <- cluster_group_df %>%
  select(Cluster, Group, Count) %>%
  pivot_wider(names_from = Group, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

library(pheatmap)

pdf(file.path(FIGURE_DIR, "cluster_group_heatmap_count.pdf"), width = 6, height = 10)
pheatmap(
  cluster_group_matrix,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
  main = "Cluster-Group Distribution (Count)",
  angle_col = 0
)
dev.off()

# 百分比热图（按cluster总数归一化）
cluster_group_pct_matrix <- cluster_group_df %>%
  select(Cluster, Group, Percentage) %>%
  pivot_wider(names_from = Group, values_from = Percentage, values_fill = 0) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

pdf(file.path(FIGURE_DIR, "cluster_group_heatmap_percentage.pdf"), width = 6, height = 10)
pheatmap(
  cluster_group_pct_matrix,
  display_numbers = TRUE,
  number_format = "%.1f",
  number_color = "black",
  fontsize_number = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
  main = "Cluster-Group Distribution (%)",
  angle_col = 0
)
dev.off()

# (6) 气泡图
p5 <- ggplot(cluster_group_df, aes(x = Group, y = Cluster)) +
  geom_point(aes(size = Count, color = Percentage)) +
  scale_size_continuous(range = c(3, 20), name = "Cell Count") +
  scale_color_gradientn(
    colors = c("lightblue", "yellow", "orange", "red"),
    name = "Percentage (%)"
  ) +
  theme_minimal() +
  labs(
    x = "Group",
    y = "Cluster",
    title = "Cluster Distribution Bubble Plot"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  geom_text(aes(label = Count), size = 3, color = "black")

save_plot(p5, "cluster_group_bubble.pdf", width = 8, height = 12)

# 4. 统计检验 ---------------------------------------------------------------
print_progress("统计显著性检验")

# 卡方检验 - 检验cluster分布与group是否独立
chisq_result <- chisq.test(cluster_group_table)

message("\n卡方检验结果:")
message(paste("Chi-squared =", round(chisq_result$statistic, 4)))
message(paste("df =", chisq_result$parameter))
message(paste("p-value =", format(chisq_result$p.value, scientific = TRUE)))

if(chisq_result$p.value < 0.05) {
  message("结论: cluster分布在MDD和SDD组之间存在显著差异")
} else {
  message("结论: cluster分布在MDD和SDD组之间无显著差异")
}

# 计算每个cluster的富集情况（相对于期望值）
expected_counts <- chisq_result$expected
observed_counts <- cluster_group_table

enrichment_df <- data.frame(
  Cluster = rep(rownames(observed_counts), ncol(observed_counts)),
  Group = rep(colnames(observed_counts), each = nrow(observed_counts)),
  Observed = as.vector(observed_counts),
  Expected = as.vector(expected_counts),
  Enrichment = as.vector(observed_counts) / as.vector(expected_counts)
)

enrichment_df$Log2_Enrichment <- log2(enrichment_df$Enrichment)
enrichment_df$Enrichment_Type <- ifelse(
  enrichment_df$Enrichment > 1.2, "Enriched",
  ifelse(enrichment_df$Enrichment < 0.8, "Depleted", "Neutral")
)

write.csv(enrichment_df,
          file = file.path(RESULT_DIR, "cluster_group_enrichment.csv"),
          row.names = FALSE)

message("\n富集分析结果:")
print(enrichment_df)

# (7) 富集热图
enrichment_matrix <- enrichment_df %>%
  select(Cluster, Group, Log2_Enrichment) %>%
  pivot_wider(names_from = Group, values_from = Log2_Enrichment) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

pdf(file.path(FIGURE_DIR, "cluster_group_enrichment_heatmap.pdf"), width = 6, height = 10)
pheatmap(
  enrichment_matrix,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-max(abs(enrichment_matrix)), max(abs(enrichment_matrix)), length.out = 101),
  main = "Cluster Enrichment in Groups\n(Log2 Observed/Expected)",
  angle_col = 0
)
dev.off()

# 5. 识别差异显著的clusters -----------------------------------------------
print_progress("识别组间差异显著的clusters")

# 对每个cluster进行Fisher精确检验或卡方检验
cluster_tests <- data.frame()

for(clust in unique(cluster_group_df$Cluster)) {
  # 构建2x2列联表：当前cluster vs 其他clusters, MDD vs SDD
  current_cluster <- cluster_group_df %>% filter(Cluster == clust)
  other_clusters <- cluster_group_df %>% filter(Cluster != clust) %>%
    group_by(Group) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  # 创建列联表
  contingency <- matrix(
    c(current_cluster$Count[current_cluster$Group == "MDD"],
      current_cluster$Count[current_cluster$Group == "SDD"],
      other_clusters$Count[other_clusters$Group == "MDD"],
      other_clusters$Count[other_clusters$Group == "SDD"]),
    nrow = 2,
    dimnames = list(c("Current", "Others"), c("MDD", "SDD"))
  )
  
  # Fisher精确检验
  fisher_result <- fisher.test(contingency)
  
  # 计算比例
  mdd_prop <- current_cluster$Count[current_cluster$Group == "MDD"] / 
    sum(current_cluster$Count)
  sdd_prop <- current_cluster$Count[current_cluster$Group == "SDD"] / 
    sum(current_cluster$Count)
  
  cluster_tests <- rbind(cluster_tests, data.frame(
    Cluster = clust,
    MDD_Count = current_cluster$Count[current_cluster$Group == "MDD"],
    SDD_Count = current_cluster$Count[current_cluster$Group == "SDD"],
    Total = sum(current_cluster$Count),
    MDD_Percentage = mdd_prop * 100,
    SDD_Percentage = sdd_prop * 100,
    Odds_Ratio = fisher_result$estimate,
    P_value = fisher_result$p.value
  ))
}

# 多重检验校正
cluster_tests$P_adjusted <- p.adjust(cluster_tests$P_value, method = "BH")
cluster_tests$Significant <- cluster_tests$P_adjusted < 0.05

# 标记富集方向
cluster_tests$Enriched_in <- ifelse(
  cluster_tests$SDD_Percentage > cluster_tests$MDD_Percentage, 
  "SDD", 
  "MDD"
)

# 按p值排序
cluster_tests <- cluster_tests %>% arrange(P_adjusted)

write.csv(cluster_tests,
          file = file.path(RESULT_DIR, "cluster_group_significance_tests.csv"),
          row.names = FALSE)

message("\n显著富集的clusters (FDR < 0.05):")
sig_clusters <- cluster_tests %>% filter(Significant)
print(sig_clusters)

# (8) 显著性可视化
cluster_tests$log10_padj <- -log10(cluster_tests$P_adjusted)
cluster_tests$Difference <- cluster_tests$SDD_Percentage - cluster_tests$MDD_Percentage

p6 <- ggplot(cluster_tests, aes(x = Difference, y = log10_padj)) +
  geom_point(aes(size = Total, color = Enriched_in), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = group_colors) +
  scale_size_continuous(range = c(3, 10), name = "Total Cells") +
  theme_classic() +
  labs(
    x = "Difference (SDD% - MDD%)",
    y = "-log10(Adjusted P-value)",
    title = "Cluster Enrichment Significance",
    color = "Enriched in"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  geom_text_repel(
    data = cluster_tests %>% filter(Significant),
    aes(label = Cluster),
    size = 4,
    max.overlaps = 20
  )

save_plot(p6, "cluster_enrichment_volcano.pdf", width = 10, height = 8)

# 6. 生成总结报告 -----------------------------------------------------------
print_progress("生成总结报告")

summary_text <- c(
  "=============================================================================",
  "软骨细胞亚群分组分布分析报告",
  "=============================================================================",
  "",
  paste("分析日期:", Sys.time()),
  paste("数据文件: chondrocyte_clustered.Rdata"),
  paste("使用分辨率:", resolution),
  "",
  "--- 基本统计 ---",
  paste("总细胞数:", ncol(chondrocyte_subset)),
  paste("聚类数量:", n_clusters),
  paste("MDD组总细胞数:", sum(chondrocyte_subset$group == "MDD")),
  paste("SDD组总细胞数:", sum(chondrocyte_subset$group == "SDD")),
  "",
  "--- 卡方检验结果 ---",
  paste("Chi-squared =", round(chisq_result$statistic, 4)),
  paste("df =", chisq_result$parameter),
  paste("P-value =", format(chisq_result$p.value, scientific = TRUE)),
  ifelse(chisq_result$p.value < 0.05,
         "结论: cluster分布在MDD和SDD组之间存在显著差异 ***",
         "结论: cluster分布在MDD和SDD组之间无显著差异"),
  "",
  "--- 显著富集的clusters (FDR < 0.05) ---",
  paste("显著cluster数量:", nrow(sig_clusters)),
  ""
)

if(nrow(sig_clusters) > 0) {
  for(i in 1:nrow(sig_clusters)) {
    summary_text <- c(
      summary_text,
      paste0("Cluster ", sig_clusters$Cluster[i], ":"),
      paste0("  富集方向: ", sig_clusters$Enriched_in[i]),
      paste0("  MDD: ", sig_clusters$MDD_Count[i], " (", 
             round(sig_clusters$MDD_Percentage[i], 1), "%)"),
      paste0("  SDD: ", sig_clusters$SDD_Count[i], " (", 
             round(sig_clusters$SDD_Percentage[i], 1), "%)"),
      paste0("  Odds Ratio: ", round(sig_clusters$Odds_Ratio[i], 3)),
      paste0("  Adjusted P-value: ", format(sig_clusters$P_adjusted[i], scientific = TRUE)),
      ""
    )
  }
} else {
  summary_text <- c(summary_text, "无显著富集的clusters")
}

summary_text <- c(
  summary_text,
  "",
  "=============================================================================",
  "输出文件:",
  paste("- 详细统计:", file.path(RESULT_DIR, "cluster_group_distribution_detailed.csv")),
  paste("- 富集分析:", file.path(RESULT_DIR, "cluster_group_enrichment.csv")),
  paste("- 显著性检验:", file.path(RESULT_DIR, "cluster_group_significance_tests.csv")),
  paste("- 可视化图表:", FIGURE_DIR),
  "============================================================================="
)

writeLines(summary_text, file.path(RESULT_DIR, "cluster_group_distribution_summary.txt"))

print_progress("分析完成！")

message("\n生成的文件:")
message(paste("- 统计结果:", RESULT_DIR))
message(paste("- 可视化图表:", FIGURE_DIR))
message(paste("- 总结报告:", file.path(RESULT_DIR, "cluster_group_distribution_summary.txt")))
