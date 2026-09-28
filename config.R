# =============================================================================
# 单细胞测序分析全局配置文件
# =============================================================================

# 项目路径设置 ---------------------------------------------------------------
PROJECT_ROOT <- "d:/Work/高天单细胞/shengxinfenxi/1116"
DATA_DIR <- file.path(PROJECT_ROOT, "临时数据")
FIGURE_DIR <- file.path(PROJECT_ROOT, "pics")
RESULT_DIR <- file.path(PROJECT_ROOT, "差异基因结果")

# 创建目录
for(dir in c(DATA_DIR, FIGURE_DIR, RESULT_DIR)) {
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# 分析参数设置 ---------------------------------------------------------------
# 质控参数
QC_PARAMS <- list(
  min_features = 200,        # 最小基因数
  max_mt_percent = 20,       # 最大线粒体比例
  max_rb_percent = 1         # 最大红细胞基因比例
)

# 聚类参数
CLUSTER_PARAMS <- list(
  n_pcs = 20,                # PCA维度数
  resolution = 0.5,          # 聚类分辨率
  resolution_range = seq(0.1, 2.0, 0.1)  # 多分辨率测试
)

# 差异分析参数
DEG_PARAMS <- list(
  logfc_threshold = 0.25,    # log2FC阈值
  pval_threshold = 0.05,     # p值阈值
  min_pct = 0.1,             # 最小表达比例
  top_n_genes = 50           # top基因数
)

# 拟时序参数
TRAJECTORY_PARAMS <- list(
  n_top_genes = 1000,        # 用于拟时序的基因数
  min_expr = 0.1,            # 最小表达量
  min_cells = 10,            # 最小细胞数
  start_cluster = "0"        # 起始簇
)

# 计算资源参数
COMPUTE_PARAMS <- list(
  n_cores = 16,              # 并行核心数
  java_memory = "20g"        # Java内存
)

# 图形参数设置 ---------------------------------------------------------------
FIGURE_PARAMS <- list(
  width = 10,
  height = 10,
  text_size = 14,
  point_size = 0.5,
  dpi = 300
)

# 样本信息 -------------------------------------------------------------------
SAMPLE_INFO <- data.frame(
  sample_id = c("ZBIVDD-1", "ZBIVDD-2", "ZBIVDD-3", 
                "ZBCON-2", "ZBCON-3", "ZBCON-4"),
  group = c("SDD", "SDD", "SDD", "MDD", "MDD", "MDD"),
  data_path = c(
    "../Chondrocytes/rawData/ZBIVDD-1/filtered_feature_bc_matrix",
    "../Chondrocytes/rawData/ZBIVDD-2/filtered_feature_bc_matrix",
    "../Chondrocytes/rawData/ZBIVDD-3/filtered_feature_bc_matrix",
    "../Chondrocytes/rawData/ZBCon-2/filtered_feature_bc_matrix",
    "../Chondrocytes/rawData/ZBCon-3/filtered_feature_bc_matrix",
    "../Chondrocytes/rawData/ZBCon-4/filtered_feature_bc_matrix"
  ),
  stringsAsFactors = FALSE
)

# 细胞类型标记基因 -----------------------------------------------------------
CELL_TYPE_MARKERS <- list(
  Chondrocyte = c("SOX5", "SOX6", "SOX9", "COL2A1", "ACAN", "COMP", "MIA", "COL9A1"),
  NP_cell = c("PAX1", "FOXF1", "KRT8", "KRT18", "CA12", "LGALS3", "CD24"),
  Notochord = c("KRT19", "SHH", "TBXT", "FOXA2", "CHRD", "NOG", "CDH2"),
  AF_cel = c("MKX", "SCX", "COL12A1", "SFRP2", "COL5A1", "COL1A1", "TNMD"),
  T_Cell = c("CD3E", "TRAC"),
  B_Plasma_Cell = c("CD79A", "MS4A1", "MZB1"),
  Endothelia = c("PECAM1", "VWF", "CDH5"),
  SMC = c("ACTA2", "MYH11", "TAGLN"),
  Erythrocyte = c("HBB", "HBA1", "CA1")
)

# 辅助函数 -------------------------------------------------------------------
# 加载常用包
load_packages <- function() {
  packages <- c("Seurat", "tidyverse", "dplyr", "ggplot2", "patchwork", 
                "harmony", "RColorBrewer", "glmGamPoi", "clustree")
  
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("安装", pkg, "..."))
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# 保存图片的标准函数
save_plot <- function(plot, filename, width = NULL, height = NULL, subdir = "pics") {
  if(is.null(width)) width <- FIGURE_PARAMS$width
  if(is.null(height)) height <- FIGURE_PARAMS$height
  
  filepath <- file.path(PROJECT_ROOT, subdir, filename)
  ggsave(filepath, plot = plot, width = width, height = height, 
         units = "in", dpi = FIGURE_PARAMS$dpi)
  message(paste("图片已保存:", filepath))
}

# 打印分析进度
print_progress <- function(message) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(message, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

message("配置文件加载完成")
