# =============================================================================
# 单细胞分析工具函数库
# =============================================================================

# 通用颜色方案函数 -----------------------------------------------------------
# 生成适用于任意数量类别的颜色
get_color_palette <- function(n) {
  if(n <= 12) {
    return(brewer.pal(max(3, n), "Set3"))
  } else if(n <= 24) {
    base_colors <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
    return(base_colors[1:n])
  } else {
    base_colors <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
    return(colorRampPalette(base_colors)(n))
  }
}

# Seurat兼容的颜色scale函数
scale_color_custom <- function() {
  scale_color_manual(values = get_color_palette(100))
}

scale_fill_custom <- function() {
  scale_fill_manual(values = get_color_palette(100))
}

# 质控可视化函数 -------------------------------------------------------------
plot_qc_metrics <- function(seurat_obj, group_by = "orig.ident", 
                             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                             save_path = NULL) {
  
  p <- VlnPlot(seurat_obj, features = features, 
               group.by = group_by, layer = "count", ncol = length(features))
  p <- p + scale_fill_custom()
  
  if(!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 10 * length(features), height = 10, units = "in")
  }
  
  return(p)
}

# 批量绘制散点图 -------------------------------------------------------------
plot_qc_scatter <- function(seurat_obj, feature_pairs = NULL, save_dir = NULL) {
  
  if(is.null(feature_pairs)) {
    feature_pairs <- list(
      c("nCount_RNA", "percent.mt"),
      c("nCount_RNA", "nFeature_RNA")
    )
  }
  
  plots <- list()
  for(i in seq_along(feature_pairs)) {
    pair <- feature_pairs[[i]]
    p <- FeatureScatter(seurat_obj, feature1 = pair[1], feature2 = pair[2]) +
      scale_color_custom()
    plots[[i]] <- p
    
    if(!is.null(save_dir)) {
      filename <- paste0(save_dir, "/scatter_", pair[1], "_vs_", pair[2], ".pdf")
      ggsave(filename, plot = p, width = 10, height = 10, units = "in")
    }
  }
  
  return(plots)
}

# 标准化预处理流程 -----------------------------------------------------------
standard_preprocess <- function(seurat_obj, n_features = 2000, 
                                 vars_to_regress = "percent.mt", 
                                 n_pcs = 50, verbose = TRUE) {
  
  if(verbose) message("开始标准化...")
  seurat_obj <- NormalizeData(seurat_obj)
  
  if(verbose) message("查找高变基因...")
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                     nfeatures = n_features)
  
  if(verbose) message("数据缩放...")
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = vars_to_regress)
  
  if(verbose) message("运行PCA...")
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
  
  return(seurat_obj)
}

# Harmony整合函数 ------------------------------------------------------------
run_harmony_integration <- function(seurat_obj, group_by = "orig.ident", 
                                    n_pcs = 20, verbose = TRUE) {
  
  if(verbose) message("运行Harmony整合...")
  
  seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = verbose
  )
  
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  
  return(seurat_obj)
}

# 聚类分析流程 ---------------------------------------------------------------
run_clustering <- function(seurat_obj, reduction = "pca", dims = 1:20,
                           resolutions = seq(0.1, 1.0, 0.1), 
                           run_umap = TRUE, run_tsne = TRUE) {
  
  message("查找邻居...")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction, dims = dims)
  
  message("聚类分析...")
  seurat_obj <- FindClusters(seurat_obj, resolution = resolutions)
  
  if(run_umap) {
    message("运行UMAP...")
    seurat_obj <- RunUMAP(seurat_obj, dims = dims, reduction = reduction)
  }
  
  if(run_tsne) {
    message("运行tSNE...")
    seurat_obj <- RunTSNE(seurat_obj, dims = dims, reduction = reduction)
  }
  
  return(seurat_obj)
}

# 绘制聚类结果总览 -----------------------------------------------------------
plot_clustering_overview <- function(seurat_obj, group_by_vars = c("orig.ident", "group"),
                                    reduction = "umap", cluster_col = "seurat_clusters",
                                    save_path = NULL) {
  
  plots <- list()
  
  # 按不同分组变量绘图
  for(var in group_by_vars) {
    if(var %in% colnames(seurat_obj@meta.data)) {
      plots[[var]] <- DimPlot(seurat_obj, reduction = reduction, group.by = var) +
        ggtitle(var) +
        scale_color_custom()
    }
  }
  
  # 按聚类标签绘图
  plots[["clusters"]] <- DimPlot(seurat_obj, reduction = reduction, 
                                 group.by = cluster_col, label = TRUE)
  
  # 组合图
  combined <- wrap_plots(plots, ncol = 2)
  
  if(!is.null(save_path)) {
    ggsave(save_path, plot = combined, width = 15, height = 15, units = "in")
  }
  
  return(combined)
}

# 细胞类型注释函数 -----------------------------------------------------------
annotate_cell_types <- function(seurat_obj, cluster_col, annotation_df) {
  # annotation_df: 包含 ClusterID 和 celltype 列的数据框
  
  seurat_obj@meta.data$celltype <- "unknown"
  
  for(i in 1:nrow(annotation_df)) {
    cluster_id <- as.character(annotation_df$ClusterID[i])
    cell_type <- annotation_df$celltype[i]
    
    cells_idx <- which(as.character(seurat_obj@meta.data[[cluster_col]]) == cluster_id)
    seurat_obj@meta.data$celltype[cells_idx] <- cell_type
  }
  
  message("细胞类型注释完成:")
  print(table(seurat_obj@meta.data$celltype))
  
  return(seurat_obj)
}

# 差异基因分析函数 -----------------------------------------------------------
find_degs <- function(seurat_obj, ident_1, ident_2 = NULL, 
                     test_use = "wilcox", logfc_threshold = 0.25,
                     min_pct = 0.1, only_pos = FALSE) {
  
  message(paste("运行差异分析:", test_use))
  
  degs <- FindMarkers(
    object = seurat_obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = test_use,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    only.pos = only_pos
  )
  
  degs$gene <- rownames(degs)
  degs <- degs %>% arrange(p_val_adj)
  
  message(paste("发现", nrow(degs), "个差异基因"))
  
  return(degs)
}

# 绘制火山图 -----------------------------------------------------------------
plot_volcano <- function(deg_df, logfc_col = "avg_log2FC", 
                        pval_col = "p_val_adj",
                        logfc_threshold = 0.25, 
                        pval_threshold = 0.05,
                        title = "Volcano Plot",
                        save_path = NULL) {
  
  # 添加分类标签
  deg_df$sig <- case_when(
    deg_df[[pval_col]] < pval_threshold & deg_df[[logfc_col]] > logfc_threshold ~ "Up",
    deg_df[[pval_col]] < pval_threshold & deg_df[[logfc_col]] < -logfc_threshold ~ "Down",
    TRUE ~ "Not Sig"
  )
  
  n_up <- sum(deg_df$sig == "Up")
  n_down <- sum(deg_df$sig == "Down")
  
  p <- ggplot(deg_df, aes_string(x = logfc_col, y = paste0("-log10(", pval_col, ")"), 
                                 color = "sig")) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "Not Sig" = "grey80")) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", color = "grey30") +
    labs(
      title = title,
      subtitle = paste0("Up: ", n_up, " | Down: ", n_down),
      x = "log2(Fold Change)",
      y = "-log10(Adjusted P-value)",
      color = "Expression"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  
  if(!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 10, height = 8, units = "in")
  }
  
  return(p)
}

# 提取表达矩阵的辅助函数 -----------------------------------------------------
extract_expression_matrix <- function(seurat_obj, features = NULL, 
                                     slot = "data", scale = FALSE) {
  
  if(is.null(features)) {
    features <- rownames(seurat_obj)
  }
  
  expr_matrix <- GetAssayData(seurat_obj, slot = slot)[features, ]
  
  if(scale) {
    expr_matrix <- t(scale(t(as.matrix(expr_matrix))))
  }
  
  return(expr_matrix)
}

# 批量保存数据 ---------------------------------------------------------------
save_analysis_data <- function(obj_list, file_prefix, save_dir = "./临时数据") {
  
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- file.path(save_dir, paste0(file_prefix, "_", timestamp, ".Rdata"))
  
  save(list = names(obj_list), file = filename, envir = list2env(obj_list))
  message(paste("数据已保存:", filename))
  
  return(filename)
}

# 生成分析报告 ---------------------------------------------------------------
generate_analysis_report <- function(seurat_obj, analysis_name = "Analysis",
                                     output_file = "analysis_report.txt") {
  
  report <- c(
    paste(rep("=", 70), collapse = ""),
    paste("分析报告:", analysis_name),
    paste("生成时间:", Sys.time()),
    paste(rep("=", 70), collapse = ""),
    "",
    "数据概况:",
    paste("- 总细胞数:", ncol(seurat_obj)),
    paste("- 总基因数:", nrow(seurat_obj)),
    "",
    "聚类信息:"
  )
  
  # 添加聚类统计
  if("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cluster_table <- table(seurat_obj@meta.data$seurat_clusters)
    report <- c(report, capture.output(print(cluster_table)))
  }
  
  # 添加细胞类型统计
  if("celltype" %in% colnames(seurat_obj@meta.data)) {
    report <- c(report, "", "细胞类型分布:")
    celltype_table <- table(seurat_obj@meta.data$celltype)
    report <- c(report, capture.output(print(celltype_table)))
  }
  
  report <- c(report, "", paste(rep("=", 70), collapse = ""))
  
  writeLines(report, output_file)
  message(paste("报告已生成:", output_file))
  
  return(invisible(report))
}

message("工具函数库加载完成")
