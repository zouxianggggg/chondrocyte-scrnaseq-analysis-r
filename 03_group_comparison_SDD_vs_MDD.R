# =============================================================================
# SDD vs MDD 组间差异分析脚本
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")
load_packages()
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)

print_progress("SDD vs MDD 组间差异分析")

# 1. 加载数据 ----------------------------------------------------------------
load(file.path(DATA_DIR, "scRNA_annotated.Rdata"))

# 确保分组信息正确
message("分组信息:")
print(table(sce$group))
print(table(sce$orig.ident, sce$group))

# 2. 细胞类型比例分析 --------------------------------------------------------
print_progress("步骤1: 细胞类型比例分析")

# 统计各细胞类型在两组中的数量
celltype_group_table <- table(sce$celltype, sce$group)
print(celltype_group_table)

# 转换为比例
celltype_prop <- prop.table(celltype_group_table, margin = 2) * 100

# 保存统计结果
write.csv(celltype_group_table, 
          file = file.path(RESULT_DIR, "celltype_count_by_group.csv"))
write.csv(celltype_prop,
          file = file.path(RESULT_DIR, "celltype_proportion_by_group.csv"))

# 可视化比例差异
library(reshape2)
prop_df <- melt(celltype_prop)
colnames(prop_df) <- c("CellType", "Group", "Percentage")

p_prop <- ggplot(prop_df, aes(x = CellType, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_custom() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type", y = "Percentage (%)", 
       title = "Cell Type Composition: SDD vs MDD") +
  geom_text(aes(label = round(Percentage, 1)), 
            position = position_dodge(0.9), vjust = -0.5, size = 3)

save_plot(p_prop, "celltype_proportion_comparison.pdf", width = 12, height = 8)

# 卡方检验
chisq_results <- list()
for(ct in rownames(celltype_group_table)) {
  if(ct != "unknown") {
    ct_table <- matrix(c(
      celltype_group_table[ct, "MDD"],
      sum(celltype_group_table[, "MDD"]) - celltype_group_table[ct, "MDD"],
      celltype_group_table[ct, "SDD"],
      sum(celltype_group_table[, "SDD"]) - celltype_group_table[ct, "SDD"]
    ), nrow = 2)
    
    test_result <- chisq.test(ct_table)
    chisq_results[[ct]] <- data.frame(
      CellType = ct,
      ChiSquare = test_result$statistic,
      PValue = test_result$p.value
    )
  }
}

chisq_df <- do.call(rbind, chisq_results)
chisq_df$Padj <- p.adjust(chisq_df$PValue, method = "BH")
chisq_df <- chisq_df[order(chisq_df$PValue), ]

write.csv(chisq_df, 
          file = file.path(RESULT_DIR, "celltype_proportion_chisq_test.csv"),
          row.names = FALSE)

message("显著差异的细胞类型 (P < 0.05):")
print(chisq_df[chisq_df$Padj < 0.05, ])

# 3. 各细胞类型内的差异基因分析 ----------------------------------------------
print_progress("步骤2: 各细胞类型内的差异基因分析")

# 获取主要细胞类型（排除unknown和细胞数太少的）
min_cells <- 50
celltypes_to_analyze <- names(which(table(sce$celltype) >= min_cells))
celltypes_to_analyze <- celltypes_to_analyze[celltypes_to_analyze != "unknown"]

message(paste("将分析以下细胞类型:", paste(celltypes_to_analyze, collapse = ", ")))

# 为每种细胞类型进行差异分析
deg_results_list <- list()

for(ct in celltypes_to_analyze) {
  message(paste("\n=== 分析细胞类型:", ct, "==="))
  
  # 提取该细胞类型
  ct_cells <- WhichCells(sce, idents = ct)
  ct_subset <- subset(sce, cells = ct_cells)
  
  # 检查两组中是否都有足够的细胞
  group_table <- table(ct_subset$group)
  message(paste("MDD:", group_table["MDD"], "细胞"))
  message(paste("SDD:", group_table["SDD"], "细胞"))
  
  if(all(group_table >= 10)) {
    # 设置分组为ident
    Idents(ct_subset) <- "group"
    
    # 差异分析
    tryCatch({
      degs <- FindMarkers(
        ct_subset,
        ident.1 = "SDD",
        ident.2 = "MDD",
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        verbose = FALSE
      )
      
      if(nrow(degs) > 0) {
        degs$gene <- rownames(degs)
        degs$celltype <- ct
        degs$change <- ifelse(degs$avg_log2FC > 0, "Up_in_SDD", "Down_in_SDD")
        
        deg_results_list[[ct]] <- degs
        
        # 保存每个细胞类型的结果
        write.csv(degs, 
                  file = file.path(RESULT_DIR, paste0("DEG_", ct, "_SDD_vs_MDD.csv")),
                  row.names = FALSE)
        
        message(paste("发现", nrow(degs), "个差异基因"))
        message(paste("  上调:", sum(degs$avg_log2FC > 0)))
        message(paste("  下调:", sum(degs$avg_log2FC < 0)))
        message(paste("  显著(p_adj<0.05):", sum(degs$p_val_adj < 0.05)))
      }
    }, error = function(e) {
      message(paste("分析失败:", e$message))
    })
  } else {
    message("细胞数不足，跳过该细胞类型")
  }
}

# 合并所有细胞类型的结果
if(length(deg_results_list) > 0) {
  all_degs <- do.call(rbind, deg_results_list)
  write.csv(all_degs,
            file = file.path(RESULT_DIR, "DEG_all_celltypes_SDD_vs_MDD.csv"),
            row.names = FALSE)
  
  # 统计每种细胞类型的DEG数量
  deg_summary <- all_degs %>%
    filter(p_val_adj < 0.05) %>%
    group_by(celltype, change) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = change, values_from = count, values_fill = 0)
  
  write.csv(deg_summary,
            file = file.path(RESULT_DIR, "DEG_summary_by_celltype.csv"),
            row.names = FALSE)
  
  # 可视化DEG数量
  deg_summary_long <- all_degs %>%
    filter(p_val_adj < 0.05) %>%
    group_by(celltype, change) %>%
    summarise(count = n(), .groups = "drop")
  
  p_deg_count <- ggplot(deg_summary_long, aes(x = celltype, y = count, fill = change)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Up_in_SDD" = "#E64B35", "Down_in_SDD" = "#4DBBD5")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cell Type", y = "Number of DEGs", 
         title = "Differential Genes: SDD vs MDD",
         fill = "Direction") +
    geom_text(aes(label = count), position = position_dodge(0.9), vjust = -0.5)
  
  save_plot(p_deg_count, "DEG_count_by_celltype.pdf", width = 12, height = 8)
}

# 4. Pseudobulk分析（更严格的方法）-------------------------------------------
print_progress("步骤3: Pseudobulk差异分析")

# 以软骨细胞为例
if("Chondrocyte" %in% celltypes_to_analyze) {
  message("对软骨细胞进行Pseudobulk分析...")
  
  # 提取软骨细胞
  chondro_cells <- WhichCells(sce, expression = celltype == "Chondrocyte")
  chondro_subset <- subset(sce, cells = chondro_cells)
  
  # 聚合每个样本的表达
  chondro_subset$sample <- chondro_subset$orig.ident
  
  # 创建pseudobulk计数矩阵
  pseudobulk_matrix <- AggregateExpression(
    chondro_subset,
    group.by = "sample",
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )$RNA
  
  # 创建样本信息表
  sample_info <- unique(chondro_subset@meta.data[, c("sample", "group")])
  rownames(sample_info) <- sample_info$sample
  sample_info <- sample_info[colnames(pseudobulk_matrix), ]
  
  # 运行DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk_matrix,
    colData = sample_info,
    design = ~ group
  )
  
  # 过滤低表达基因
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # 运行差异分析
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", "SDD", "MDD"))
  
  # 转换为数据框
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$padj), ]
  
  # 保存结果
  write.csv(res_df,
            file = file.path(RESULT_DIR, "DEG_Chondrocyte_pseudobulk_DESeq2.csv"),
            row.names = FALSE)
  
  # 可视化（火山图）
  p_volcano_pb <- plot_volcano(
    res_df,
    logfc_col = "log2FoldChange",
    pval_col = "padj",
    title = "Chondrocyte: SDD vs MDD (Pseudobulk)"
  )
  save_plot(p_volcano_pb, "volcano_Chondrocyte_pseudobulk.pdf")
  
  # MA图
  pdf(file.path(FIGURE_DIR, "MA_plot_Chondrocyte_pseudobulk.pdf"), 
      width = 8, height = 8)
  DESeq2::plotMA(res, alpha = 0.05)
  dev.off()
  
  message(paste("Pseudobulk分析发现", sum(res$padj < 0.05, na.rm = TRUE), "个显著DEG"))
}

# 5. 功能富集分析 ------------------------------------------------------------
print_progress("步骤4: 功能富集分析")

# 对软骨细胞的DEG进行GO富集
if("Chondrocyte" %in% names(deg_results_list)) {
  chondro_degs <- deg_results_list[["Chondrocyte"]]
  sig_degs <- chondro_degs %>% filter(p_val_adj < 0.05)
  
  if(nrow(sig_degs) > 10) {
    # 分别分析上调和下调基因
    up_genes <- sig_degs %>% filter(avg_log2FC > 0) %>% pull(gene)
    down_genes <- sig_degs %>% filter(avg_log2FC < 0) %>% pull(gene)
    
    # GO富集 - 上调基因
    if(length(up_genes) >= 5) {
      tryCatch({
        ego_up <- enrichGO(
          gene = up_genes,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        if(!is.null(ego_up) && nrow(ego_up) > 0) {
          write.csv(as.data.frame(ego_up),
                    file = file.path(RESULT_DIR, "GO_enrichment_Chondrocyte_up_in_SDD.csv"),
                    row.names = FALSE)
          
          p_go_up <- dotplot(ego_up, showCategory = 20) + 
            ggtitle("GO Enrichment: Upregulated in SDD (Chondrocyte)")
          save_plot(p_go_up, "GO_Chondrocyte_up.pdf", width = 12, height = 10)
        }
      }, error = function(e) {
        message("GO富集分析失败（上调基因）:", e$message)
      })
    }
    
    # GO富集 - 下调基因
    if(length(down_genes) >= 5) {
      tryCatch({
        ego_down <- enrichGO(
          gene = down_genes,
          OrgDb = org.Hs.eg.db,
          keyType = "SYMBOL",
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        
        if(!is.null(ego_down) && nrow(ego_down) > 0) {
          write.csv(as.data.frame(ego_down),
                    file = file.path(RESULT_DIR, "GO_enrichment_Chondrocyte_down_in_SDD.csv"),
                    row.names = FALSE)
          
          p_go_down <- dotplot(ego_down, showCategory = 20) + 
            ggtitle("GO Enrichment: Downregulated in SDD (Chondrocyte)")
          save_plot(p_go_down, "GO_Chondrocyte_down.pdf", width = 12, height = 10)
        }
      }, error = function(e) {
        message("GO富集分析失败（下调基因）:", e$message)
      })
    }
    
    # KEGG富集
    tryCatch({
      # 转换为Entrez ID
      gene_entrez <- bitr(sig_degs$gene, 
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)
      
      kegg_result <- enrichKEGG(
        gene = gene_entrez$ENTREZID,
        organism = "hsa",
        pvalueCutoff = 0.05
      )
      
      if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
        write.csv(as.data.frame(kegg_result),
                  file = file.path(RESULT_DIR, "KEGG_enrichment_Chondrocyte.csv"),
                  row.names = FALSE)
        
        p_kegg <- dotplot(kegg_result, showCategory = 20) + 
          ggtitle("KEGG Pathway Enrichment (Chondrocyte)")
        save_plot(p_kegg, "KEGG_Chondrocyte.pdf", width = 12, height = 10)
      }
    }, error = function(e) {
      message("KEGG富集分析失败:", e$message)
    })
  }
}

# 6. 热图展示关键基因 --------------------------------------------------------
print_progress("步骤5: 关键基因热图")

if("Chondrocyte" %in% names(deg_results_list)) {
  chondro_degs <- deg_results_list[["Chondrocyte"]]
  top_degs <- chondro_degs %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(50) %>%
    pull(gene)
  
  if(length(top_degs) > 0) {
    # 提取软骨细胞
    chondro_cells <- WhichCells(sce, expression = celltype == "Chondrocyte")
    chondro_subset <- subset(sce, cells = chondro_cells)
    
    # 提取表达矩阵
    expr_matrix <- extract_expression_matrix(
      chondro_subset, 
      features = top_degs,
      scale = TRUE
    )
    
    # 准备注释
    col_anno <- data.frame(
      Group = chondro_subset$group,
      row.names = colnames(chondro_subset)
    )
    
    # 绘制热图
    pdf(file.path(FIGURE_DIR, "heatmap_Chondrocyte_top50_DEGs.pdf"),
        width = 12, height = 10)
    
    pheatmap(
      expr_matrix,
      color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      annotation_col = col_anno,
      annotation_colors = list(Group = c("MDD" = "#00A087", "SDD" = "#DC0000")),
      fontsize_row = 6,
      main = "Top 50 DEGs in Chondrocyte: SDD vs MDD"
    )
    
    dev.off()
  }
}

print_progress("SDD vs MDD 组间差异分析完成")

message("\n结果文件保存于:")
message(paste("- ", RESULT_DIR))

