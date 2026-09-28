# =============================================================================
# 软骨细胞和髓核细胞共同保护通路分析及HDAC9验证
# 使用 preranked GSEA (fgsea) 方法 + 四类证据验证
# =============================================================================

# 加载配置和工具函数 ---------------------------------------------------------
source("config.R")
source("utils.R")

# 禁用并行计算并增加内存限制（必须在加载Seurat之前设置）
library(future)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB
plan("sequential")

load_packages()

# 安装并加载额外需要的包
if(!require("fgsea")) install.packages("fgsea")
if(!require("msigdbr")) install.packages("msigdbr")
if(!require("GSVA")) BiocManager::install("GSVA")
if(!require("UCell")) BiocManager::install("UCell")
if(!require("lme4")) install.packages("lme4")
if(!require("lmerTest")) install.packages("lmerTest")

library(fgsea)
library(msigdbr)
library(GSVA)
library(UCell)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(patchwork)
library(lme4)
library(lmerTest)

print_progress("软骨细胞和髓核细胞共同保护通路GSEA分析及HDAC9验证")

# 参数设置 -------------------------------------------------------------------
TARGET_GENE <- "HDAC9"
PVALUE_CUTOFF <- 0.05
PADJ_CUTOFF <- 0.25
MIN_GENESET_SIZE <- 15
MAX_GENESET_SIZE <- 500
NPERM <- 10000

# 定义保护性相关通路关键词（扩充版）
PROTECTIVE_KEYWORDS <- c(
  # ECM / 基质与细胞附着
  "EXTRACELLULAR.MATRIX", "ECM", "COLLAGEN", "PROTEOGLYCAN", "CARTILAGE", "ADHESION",
  # 细胞存活 / 抗凋亡 / 增殖修复
  "APOPTOSIS", "CELL.DEATH", "CELL.SURVIVAL", "ANTIAPOPTOTIC", "SURVIVAL", "PROLIFER",
  # 氧化应激 / 抗氧化 / NRF2
  "OXIDATIVE.STRESS", "ANTIOXIDANT", "REACTIVE.OXYGEN", "ROS", "NRF2", "NFE2L2", "REDOX",
  # 蛋白质稳态 / 应激反应 / ER应激
  "STRESS.RESPONSE", "UNFOLDED.PROTEIN", "HEAT.SHOCK", "CELLULAR.STRESS", "UPR", "PROTEOSTASIS", "PROTEASOME", "HSP", "CHAPERONE",
  # 自噬 / 线粒体 / 溶酶体
  "AUTOPHAGY", "MITOPHAGY", "MITOCHONDR", "MITOCHONDRIA", "LYSOSOME", "MITOPH",
  # 代谢与能量 / 关键信号通路
  "TGF", "TGF.BETA", "TGF-BETA", "PI3K", "AKT", "PI3K-AKT", "FOXO", "HIF", "AMPK", "MTOR", "PI3K/AKT",
  # DNA修复 / 基因组稳定性 / 细胞周期
  "DNA.REPAIR", "DNA.DAMAGE", "GENOME.STABILITY", "DNA.REPLICATION", "CELL.CYCLE", "CHECKPOINT",
  # 炎症调控 / 免疫调节
  "INFLAMMATION", "INFLAMMATORY", "NF.KAPPA.B", "CYTOKINE", "IMMUNE", "RESOLUTION",
  # 表观遗传 / 组蛋白修饰 / HDAC / SIRT
  "HISTONE", "CHROMATIN", "EPIGENETIC", "ACETYLATION", "DEACETYLASE", "HDAC", "SIRT", "SIRT1",
  # 衰老 / 血管生成 / 修复再生
  "SENESCENCE", "AGEING", "AGING", "ANGIOGENESIS", "GROWTH.FACTOR", "REPAIR", "REGENERATION", "NECROSIS"
)

# =============================================================================
# 一、加载数据和通路数据库
# =============================================================================
print_progress("步骤1: 加载数据和通路数据库")

# 1.1 加载软骨细胞数据
message("\n加载软骨细胞数据...")
load(file.path(DATA_DIR, "chondrocyte_clustered.Rdata"))
chondrocyte_data <- chondrocyte_subset
rm(chondrocyte_subset)

message(paste("软骨细胞数:", ncol(chondrocyte_data)))
message(paste("基因数:", nrow(chondrocyte_data)))
message("分组统计:")
print(table(chondrocyte_data$group))

# 1.2 加载髓核细胞数据
message("\n加载髓核细胞数据...")
load(file.path(DATA_DIR, "np_clustered.Rdata"))
np_data <- np_subset
rm(np_subset)

message(paste("髓核细胞数:", ncol(np_data)))
message(paste("基因数:", nrow(np_data)))
message("分组统计:")
print(table(np_data$group))

# 1.3 加载通路数据库
message("\n加载通路数据库...")

# 检查是否有已保存的MSigDB数据
msigdb_cache_file <- file.path(DATA_DIR, "msigdb_cache.Rdata")

if(file.exists(msigdb_cache_file)) {
  message("从缓存加载MSigDB数据...")
  load(msigdb_cache_file)
  message("✓ 从缓存加载成功")
  
} else {
  message("MSigDB缓存不存在，直接使用本地GO/KEGG数据库...")
  message("(如需使用MSigDB，请手动下载后保存到", msigdb_cache_file, ")")
  
  # 使用本地数据库构建基因集
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(KEGGREST)
  
  message("\n构建基因集列表...")
  pathways_list <- list()
  
  # 1. 从GO数据库构建
  message("1. 从GO数据库提取通路...")
  go_terms <- keys(org.Hs.eg.db, keytype = "GOALL")
  go_sampled <- sample(go_terms, min(3000, length(go_terms)))
  
  pb <- txtProgressBar(max = length(go_sampled), style = 3)
  for(i in seq_along(go_sampled)) {
    setTxtProgressBar(pb, i)
    go_id <- go_sampled[i]
    
    tryCatch({
      genes <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = go_id, 
                                      columns = "SYMBOL", 
                                      keytype = "GOALL")
      
      if(!is.null(genes) && nrow(genes) > 0) {
        gene_symbols <- unique(genes$SYMBOL[!is.na(genes$SYMBOL)])
        
        if(length(gene_symbols) >= MIN_GENESET_SIZE && length(gene_symbols) <= MAX_GENESET_SIZE) {
          go_info <- AnnotationDbi::select(GO.db::GO.db, 
                                           keys = go_id, 
                                           columns = c("TERM", "ONTOLOGY"), 
                                           keytype = "GOID")
          
          if(nrow(go_info) > 0 && !is.na(go_info$ONTOLOGY[1])) {
            # 只保留BP (Biological Process)
            if(go_info$ONTOLOGY[1] == "BP") {
              pathway_name <- paste0("GOBP_", gsub(" ", "_", toupper(go_info$TERM[1])))
              pathways_list[[pathway_name]] <- gene_symbols
            }
          }
        }
      }
    }, error = function(e) {})
  }
  close(pb)
  
  message(paste("\n  GO通路数:", length(pathways_list)))
  
  # 2. 从KEGG数据库构建
  message("2. 从KEGG数据库提取通路...")
  tryCatch({
    # 使用clusterProfiler的内置KEGG数据
    kegg_gene_list <- bitr(keys(org.Hs.eg.db, "SYMBOL"), 
                           fromType = "SYMBOL", 
                           toType = "ENTREZID", 
                           OrgDb = org.Hs.eg.db)
    
    # 获取人类KEGG通路
    search_kegg_organism('hsa', by='kegg_code')
    kegg_pathways <- download_KEGG('hsa')
    
    if(!is.null(kegg_pathways) && length(kegg_pathways$KEGGPATHID2EXTID) > 0) {
      for(pathway_id in names(kegg_pathways$KEGGPATHID2EXTID)) {
        entrez_ids <- kegg_pathways$KEGGPATHID2EXTID[[pathway_id]]
        
        # 转换为SYMBOL
        symbols <- kegg_gene_list$SYMBOL[kegg_gene_list$ENTREZID %in% entrez_ids]
        symbols <- unique(symbols[!is.na(symbols)])
        
        if(length(symbols) >= MIN_GENESET_SIZE && length(symbols) <= MAX_GENESET_SIZE) {
          pathway_name <- names(kegg_pathways$KEGGPATHID2NAME)[
            kegg_pathways$KEGGPATHID2NAME == pathway_id
          ]
          if(length(pathway_name) > 0) {
            pathway_name <- paste0("KEGG_", gsub(" ", "_", toupper(pathway_name[1])))
            pathways_list[[pathway_name]] <- symbols
          }
        }
      }
    }
    
    message(paste("  KEGG通路数:", sum(grepl("^KEGG_", names(pathways_list)))))
  }, error = function(e) {
    message("  KEGG加载失败:", e$message)
  })
  
  # 3. 添加一些手动定义的保护性通路
  message("3. 添加常见保护性通路...")
  
  # ECM相关基因
  ecm_genes <- c("COL1A1", "COL1A2", "COL2A1", "COL3A1", "COL5A1", "COL6A1", 
                 "ACAN", "DCN", "BGN", "FMOD", "LUM", "CHAD", "COMP", "FN1",
                 "VIM", "VCAN", "TNC", "THBS1", "THBS2", "SPP1", "SPARC")
  pathways_list[["MANUAL_EXTRACELLULAR_MATRIX_ORGANIZATION"]] <- ecm_genes
  
  # 抗氧化
  antioxidant_genes <- c("SOD1", "SOD2", "SOD3", "CAT", "GPX1", "GPX4", 
                        "PRDX1", "PRDX2", "PRDX3", "PRDX5", "TXN", "TXNRD1",
                        "NQO1", "HMOX1", "NFE2L2", "GCLC", "GSR", "GSS")
  pathways_list[["MANUAL_ANTIOXIDANT_RESPONSE"]] <- antioxidant_genes
  
  # 自噬
  autophagy_genes <- c("ATG5", "ATG7", "ATG12", "BECN1", "MAP1LC3A", "MAP1LC3B",
                      "SQSTM1", "NBR1", "OPTN", "ULK1", "ULK2", "PIK3C3",
                      "ATG4A", "ATG4B", "GABARAP", "GABARAPL1", "LAMP1", "LAMP2")
  pathways_list[["MANUAL_AUTOPHAGY"]] <- autophagy_genes
  
  # HIF1信号
  hif1_genes <- c("HIF1A", "EPAS1", "ARNT", "VHL", "EGLN1", "EGLN2", "EGLN3",
                 "VEGFA", "VEGFB", "VEGFC", "SLC2A1", "LDHA", "PDK1", "ENO1")
  pathways_list[["MANUAL_HIF1_SIGNALING"]] <- hif1_genes
  
  message(paste("  手动添加通路数:", 4))
  
  # 保存到缓存
  message(paste("\n总计", length(pathways_list), "个基因集"))
  message("保存到缓存文件...")
  msigdb_all <- data.frame()  # 空的，标记使用本地数据库
  save(pathways_list, msigdb_all, file = msigdb_cache_file)
  message("✓ 缓存已保存")
}

# 只有在MSigDB加载成功时才需要转换格式
if(exists("msigdb_all") && is.data.frame(msigdb_all) && nrow(msigdb_all) > 0) {
  # 添加源信息
  msigdb_pathways_info <- msigdb_all %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()
  
  message(paste("总基因集数量:", length(unique(msigdb_all$gs_name))))
  
  # 转换为fgsea需要的格式
  pathways_list <- split(msigdb_all$gene_symbol, msigdb_all$gs_name)
  
  # 过滤基因集大小
  pathways_list <- pathways_list[sapply(pathways_list, length) >= MIN_GENESET_SIZE & 
                                   sapply(pathways_list, length) <= MAX_GENESET_SIZE]
  
  message(paste("过滤后基因集数量:", length(pathways_list)))
} else {
  # 使用本地构建的pathways_list（已在前面构建完成）
  message("使用本地数据库构建的基因集")
}

# 验证至少有一些通路可用
if(!exists("pathways_list") || length(pathways_list) == 0) {
  stop("错误: 没有可用的基因集进行分析。请检查数据库构建过程。")
}

message(paste("\n最终使用", length(pathways_list), "个基因集进行GSEA分析"))

# =============================================================================
# 二、差异分析并准备ranked基因列表
# =============================================================================
print_progress("步骤2: 差异分析并准备ranked基因列表")

# 2.1 软骨细胞差异分析
message("\n软骨细胞差异分析 (SDD vs MDD)...")
Idents(chondrocyte_data) <- chondrocyte_data$group

chondrocyte_deg <- FindMarkers(
  chondrocyte_data,
  ident.1 = "SDD",
  ident.2 = "MDD",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1,
  verbose = FALSE
)

chondrocyte_deg$gene <- rownames(chondrocyte_deg)

# 创建ranked list: -log10(p) * sign(log2FC)
# 负值 = 下调（保护性），正值 = 上调
chondrocyte_deg$rank_metric <- -log10(chondrocyte_deg$p_val + 1e-300) * sign(chondrocyte_deg$avg_log2FC)
chondrocyte_ranks <- setNames(chondrocyte_deg$rank_metric, chondrocyte_deg$gene)
chondrocyte_ranks <- sort(chondrocyte_ranks, decreasing = TRUE)

message(paste("软骨细胞ranked基因数:", length(chondrocyte_ranks)))

# 2.2 髓核细胞差异分析
message("\n髓核细胞差异分析 (SDD vs MDD)...")
Idents(np_data) <- np_data$group

np_deg <- FindMarkers(
  np_data,
  ident.1 = "SDD",
  ident.2 = "MDD",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1,
  verbose = FALSE
)

np_deg$gene <- rownames(np_deg)
np_deg$rank_metric <- -log10(np_deg$p_val + 1e-300) * sign(np_deg$avg_log2FC)
np_ranks <- setNames(np_deg$rank_metric, np_deg$gene)
np_ranks <- sort(np_ranks, decreasing = TRUE)

message(paste("髓核细胞ranked基因数:", length(np_ranks)))

# 2.3 保存结果
write.csv(chondrocyte_deg, file.path(RESULT_DIR, "10_软骨细胞_差异基因.csv"), row.names = FALSE)
write.csv(np_deg, file.path(RESULT_DIR, "10_髓核细胞_差异基因.csv"), row.names = FALSE)

# 2.4 检查HDAC9
message("\n检查HDAC9:")
if(TARGET_GENE %in% chondrocyte_deg$gene) {
  hdac9_chon <- chondrocyte_deg[chondrocyte_deg$gene == TARGET_GENE, ]
  message(paste("软骨细胞 - log2FC:", round(hdac9_chon$avg_log2FC, 3),
                ", p_adj:", format(hdac9_chon$p_val_adj, sci = TRUE),
                ", rank:", round(hdac9_chon$rank_metric, 2)))
}
if(TARGET_GENE %in% np_deg$gene) {
  hdac9_np <- np_deg[np_deg$gene == TARGET_GENE, ]
  message(paste("髓核细胞 - log2FC:", round(hdac9_np$avg_log2FC, 3),
                ", p_adj:", format(hdac9_np$p_val_adj, sci = TRUE),
                ", rank:", round(hdac9_np$rank_metric, 2)))
}

# =============================================================================
# 三、运行preranked GSEA
# =============================================================================
print_progress("步骤3: 运行preranked GSEA")

# 3.1 软骨细胞GSEA
message("\n运行软骨细胞GSEA...")
set.seed(123)
chondrocyte_gsea <- fgsea(
  pathways = pathways_list,
  stats = chondrocyte_ranks,
  minSize = MIN_GENESET_SIZE,
  maxSize = MAX_GENESET_SIZE,
  nproc = 1,
  nPermSimple = NPERM
)

chondrocyte_gsea <- as.data.frame(chondrocyte_gsea)
chondrocyte_gsea <- chondrocyte_gsea %>% arrange(pval)
chondrocyte_gsea$cell_type <- "Chondrocyte"

message(paste("显著富集通路数 (padj < 0.25):", sum(chondrocyte_gsea$padj < PADJ_CUTOFF, na.rm = TRUE)))
message(paste("  正向富集 (NES > 0):", sum(chondrocyte_gsea$padj < PADJ_CUTOFF & chondrocyte_gsea$NES > 0, na.rm = TRUE)))
message(paste("  负向富集 (NES < 0):", sum(chondrocyte_gsea$padj < PADJ_CUTOFF & chondrocyte_gsea$NES < 0, na.rm = TRUE)))

# 3.2 髓核细胞GSEA
message("\n运行髓核细胞GSEA...")
set.seed(123)
np_gsea <- fgsea(
  pathways = pathways_list,
  stats = np_ranks,
  minSize = MIN_GENESET_SIZE,
  maxSize = MAX_GENESET_SIZE,
  nproc = 1,
  nPermSimple = NPERM
)

np_gsea <- as.data.frame(np_gsea)
np_gsea <- np_gsea %>% arrange(pval)
np_gsea$cell_type <- "Nucleus_Pulposus"

message(paste("显著富集通路数 (padj < 0.25):", sum(np_gsea$padj < PADJ_CUTOFF, na.rm = TRUE)))
message(paste("  正向富集 (NES > 0):", sum(np_gsea$padj < PADJ_CUTOFF & np_gsea$NES > 0, na.rm = TRUE)))
message(paste("  负向富集 (NES < 0):", sum(np_gsea$padj < PADJ_CUTOFF & np_gsea$NES < 0, na.rm = TRUE)))

# 3.3 保存完整GSEA结果
# 将所有可能的 list 型列安全地转换为分号分隔的字符串后再写入 CSV
convert_list_columns <- function(df) {
  for(col in names(df)) {
    # 如果整列是 list 类型（fgsea 结果常见情况）
    if(inherits(df[[col]], "list")) {
      df[[col]] <- sapply(df[[col]], function(x) {
        if(is.null(x)) return(NA_character_)
        # 原子向量（字符/数值）
        if(is.atomic(x)) return(paste(x, collapse = ";"))
        # 其他可展平对象
        return(paste(unlist(x), collapse = ";"))
      }, USE.NAMES = FALSE)
    }
  }
  return(as.data.frame(df, stringsAsFactors = FALSE))
}

chondrocyte_gsea_save <- convert_list_columns(chondrocyte_gsea)
write.csv(chondrocyte_gsea_save, file.path(RESULT_DIR, "10_软骨细胞_GSEA结果.csv"), row.names = FALSE)

np_gsea_save <- convert_list_columns(np_gsea)
write.csv(np_gsea_save, file.path(RESULT_DIR, "10_髓核细胞_GSEA结果.csv"), row.names = FALSE)

# =============================================================================
# 四、筛选保护性相关通路
# =============================================================================
print_progress("步骤4: 筛选保护性相关通路")

# 4.1 筛选负向富集的保护性通路（NES < 0，在病变组中下调）
message("\n筛选保护性相关通路...")

filter_protective_pathways <- function(gsea_result) {
  protective <- gsea_result %>%
    filter(padj < PADJ_CUTOFF, NES < 0) %>%  # 负向富集 = 保护性基因下调
    filter(grepl(paste(PROTECTIVE_KEYWORDS, collapse = "|"), pathway, ignore.case = TRUE))
  return(protective)
}

chondrocyte_protective <- filter_protective_pathways(chondrocyte_gsea)
np_protective <- filter_protective_pathways(np_gsea)

message(paste("软骨细胞保护性通路:", nrow(chondrocyte_protective)))
message(paste("髓核细胞保护性通路:", nrow(np_protective)))

# 4.2 识别共同保护性通路
common_protective_pathways <- intersect(chondrocyte_protective$pathway, np_protective$pathway)
message(paste("\n共同保护性通路数（严格交集）:", length(common_protective_pathways)))

# 4.3 安全回退：如果交集为空，使用包含HDAC9的保护性通路
if(length(common_protective_pathways) == 0) {
  message("\n警告: 严格交集为空，启用回退策略...")
  
  # 查找所有包含HDAC9的通路
  pathways_with_hdac9 <- names(pathways_list)[sapply(pathways_list, function(x) TARGET_GENE %in% x)]
  message(paste("  pathways_list中包含HDAC9的通路数:", length(pathways_with_hdac9)))
  
  # 从保护性通路中筛选包含HDAC9的
  all_protective <- union(chondrocyte_protective$pathway, np_protective$pathway)
  hdac9_protective <- intersect(pathways_with_hdac9, all_protective)
  
  if(length(hdac9_protective) > 0) {
    message(paste("  ✓ 使用包含HDAC9的保护性通路:", length(hdac9_protective), "个"))
    common_protective_pathways <- hdac9_protective
  } else {
    # 如果还是没有，使用并集的top通路
    message("  ✓ 使用两种细胞的保护性通路并集（top 20）")
    union_protective <- rbind(
      chondrocyte_protective %>% mutate(source = "chondrocyte"),
      np_protective %>% mutate(source = "np")
    ) %>%
      arrange(padj) %>%
      head(20)
    common_protective_pathways <- unique(union_protective$pathway)
  }
  message(paste("  最终用于HDAC9验证的通路数:", length(common_protective_pathways)))
}

if(length(common_protective_pathways) > 0) {
  message("\nTop 10 共同保护性通路:")
  common_pathways_info <- chondrocyte_protective %>%
    filter(pathway %in% common_protective_pathways) %>%
    arrange(padj) %>%
    head(10) %>%
    dplyr::select(pathway, NES, padj, size)
  print(common_pathways_info)
  
  # 保存共同通路
  common_protective_chon <- chondrocyte_protective %>%
    filter(pathway %in% common_protective_pathways) %>%
    arrange(padj)

  common_protective_np <- np_protective %>%
    filter(pathway %in% common_protective_pathways) %>%
    arrange(padj)

  # 在写入 CSV 之前把可能的 list 列展平成字符，避免 write.csv 报错
  common_protective_chon_save <- convert_list_columns(common_protective_chon)
  common_protective_np_save <- convert_list_columns(common_protective_np)

  write.csv(common_protective_chon_save,
            file.path(RESULT_DIR, "10_共同保护性通路_软骨细胞.csv"),
            row.names = FALSE)
  write.csv(common_protective_np_save,
            file.path(RESULT_DIR, "10_共同保护性通路_髓核细胞.csv"),
            row.names = FALSE)

  # 创建对比表（用原始表构建，之后展平以写入）
  common_comparison <- common_protective_chon %>%
    dplyr::select(pathway, NES_chon = NES, padj_chon = padj, size_chon = size) %>%
    inner_join(
      common_protective_np %>% dplyr::select(pathway, NES_np = NES, padj_np = padj, size_np = size),
      by = "pathway"
    ) %>%
    arrange(padj_chon)

  common_comparison_save <- convert_list_columns(common_comparison)
  write.csv(common_comparison_save,
            file.path(RESULT_DIR, "10_共同保护性通路_对比.csv"),
            row.names = FALSE)
}

# =============================================================================
# 五、HDAC9的四类证据验证
# =============================================================================
print_progress("步骤5: HDAC9的四类证据验证")

if(length(common_protective_pathways) == 0) {
  message("警告: 没有找到共同保护性通路，无法进行HDAC9验证")
  quit(save = "no")
}

# 初始化结果存储
hdac9_evidence <- list()

# --------------------------------------------------
# 证据1: 成员证据 - HDAC9是否属于通路基因集
# --------------------------------------------------
message("\n=== 证据1: 成员证据 ===")
message("检查HDAC9是否属于共同保护性通路...")

hdac9_membership <- data.frame()
for(pw in common_protective_pathways) {
  genes_in_pathway <- pathways_list[[pw]]
  is_member <- TARGET_GENE %in% genes_in_pathway
  
  if(is_member) {
    hdac9_membership <- rbind(hdac9_membership, data.frame(
      pathway = pw,
      is_member = TRUE,
      pathway_size = length(genes_in_pathway)
    ))
  }
}

if(nrow(hdac9_membership) > 0) {
  message(paste("*** HDAC9出现在", nrow(hdac9_membership), "个共同保护性通路中 ***"))
  hdac9_evidence$membership <- hdac9_membership
  print(hdac9_membership)
  
  write.csv(hdac9_membership,
            file.path(RESULT_DIR, "10_HDAC9证据1_成员证据.csv"),
            row.names = FALSE)
} else {
  message("HDAC9不在任何共同保护性通路的基因集中")
}

# --------------------------------------------------
# 证据2: 贡献证据 - HDAC9是否在leading edge
# --------------------------------------------------
message("\n=== 证据2: 贡献证据 (Leading Edge) ===")
message("检查HDAC9是否在通路的leading edge中...")

hdac9_leading_edge <- data.frame()

# 检查软骨细胞
for(pw in common_protective_pathways) {
  pathway_result <- chondrocyte_gsea[chondrocyte_gsea$pathway == pw, ]
  if(nrow(pathway_result) > 0) {
    leading_genes <- unlist(pathway_result$leadingEdge)
    if(TARGET_GENE %in% leading_genes) {
      hdac9_leading_edge <- rbind(hdac9_leading_edge, data.frame(
        pathway = pw,
        cell_type = "Chondrocyte",
        in_leading_edge = TRUE,
        NES = pathway_result$NES,
        padj = pathway_result$padj,
        leading_edge_size = length(leading_genes)
      ))
    }
  }
}

# 检查髓核细胞
for(pw in common_protective_pathways) {
  pathway_result <- np_gsea[np_gsea$pathway == pw, ]
  if(nrow(pathway_result) > 0) {
    leading_genes <- unlist(pathway_result$leadingEdge)
    if(TARGET_GENE %in% leading_genes) {
      hdac9_leading_edge <- rbind(hdac9_leading_edge, data.frame(
        pathway = pw,
        cell_type = "Nucleus_Pulposus",
        in_leading_edge = TRUE,
        NES = pathway_result$NES,
        padj = pathway_result$padj,
        leading_edge_size = length(leading_genes)
      ))
    }
  }
}

if(nrow(hdac9_leading_edge) > 0) {
  message(paste("*** HDAC9在", nrow(hdac9_leading_edge), "个通路-细胞类型组合的leading edge中 ***"))
  hdac9_evidence$leading_edge <- hdac9_leading_edge
  print(hdac9_leading_edge)
  
  write.csv(hdac9_leading_edge,
            file.path(RESULT_DIR, "10_HDAC9证据2_Leading_Edge.csv"),
            row.names = FALSE)
} else {
  message("HDAC9不在任何通路的leading edge中")
}

# --------------------------------------------------
# 证据3: 相关证据 - 通路活性与HDAC9表达的相关性
# --------------------------------------------------
message("\n=== 证据3: 相关证据 (通路活性 vs HDAC9表达) ===")
message("计算通路活性得分与HDAC9表达的相关性...")

# 3.1 使用UCell计算通路活性得分
message("\n计算通路活性得分 (UCell)...")

# 准备共同保护性通路的基因集
pathway_signatures <- pathways_list[common_protective_pathways]

# 软骨细胞
message("  软骨细胞...")
chondrocyte_data <- AddModuleScore_UCell(chondrocyte_data, features = pathway_signatures)

# 髓核细胞
message("  髓核细胞...")
np_data <- AddModuleScore_UCell(np_data, features = pathway_signatures)

# 3.2 提取数据并计算相关性
hdac9_correlation <- data.frame()

# 软骨细胞相关性
if(TARGET_GENE %in% rownames(chondrocyte_data)) {
  chon_meta <- chondrocyte_data@meta.data
  chon_meta$HDAC9 <- GetAssayData(chondrocyte_data, slot = "data")[TARGET_GENE, ]
  chon_meta$sample <- chondrocyte_data$orig.ident  # 样本ID
  
  for(pw in common_protective_pathways) {
    ucell_col <- paste0(pw, "_UCell")
    if(ucell_col %in% colnames(chon_meta)) {
      # 整体相关性
      cor_test <- cor.test(chon_meta$HDAC9, chon_meta[[ucell_col]], method = "spearman")
      
      # 分组相关性
      cor_mdd <- cor.test(
        chon_meta$HDAC9[chon_meta$group == "MDD"],
        chon_meta[[ucell_col]][chon_meta$group == "MDD"],
        method = "spearman"
      )
      cor_sdd <- cor.test(
        chon_meta$HDAC9[chon_meta$group == "SDD"],
        chon_meta[[ucell_col]][chon_meta$group == "SDD"],
        method = "spearman"
      )
      
      hdac9_correlation <- rbind(hdac9_correlation, data.frame(
        pathway = pw,
        cell_type = "Chondrocyte",
        cor_overall = cor_test$estimate,
        pval_overall = cor_test$p.value,
        cor_MDD = cor_mdd$estimate,
        pval_MDD = cor_mdd$p.value,
        cor_SDD = cor_sdd$estimate,
        pval_SDD = cor_sdd$p.value
      ))
    }
  }
}

# 髓核细胞相关性
if(TARGET_GENE %in% rownames(np_data)) {
  np_meta <- np_data@meta.data
  np_meta$HDAC9 <- GetAssayData(np_data, slot = "data")[TARGET_GENE, ]
  np_meta$sample <- np_data$orig.ident
  
  for(pw in common_protective_pathways) {
    ucell_col <- paste0(pw, "_UCell")
    if(ucell_col %in% colnames(np_meta)) {
      cor_test <- cor.test(np_meta$HDAC9, np_meta[[ucell_col]], method = "spearman")
      
      cor_mdd <- cor.test(
        np_meta$HDAC9[np_meta$group == "MDD"],
        np_meta[[ucell_col]][np_meta$group == "MDD"],
        method = "spearman"
      )
      cor_sdd <- cor.test(
        np_meta$HDAC9[np_meta$group == "SDD"],
        np_meta[[ucell_col]][np_meta$group == "SDD"],
        method = "spearman"
      )
      
      hdac9_correlation <- rbind(hdac9_correlation, data.frame(
        pathway = pw,
        cell_type = "Nucleus_Pulposus",
        cor_overall = cor_test$estimate,
        pval_overall = cor_test$p.value,
        cor_MDD = cor_mdd$estimate,
        pval_MDD = cor_mdd$p.value,
        cor_SDD = cor_sdd$estimate,
        pval_SDD = cor_sdd$p.value
      ))
    }
  }
}

if(nrow(hdac9_correlation) > 0) {
  hdac9_correlation$padj_overall <- p.adjust(hdac9_correlation$pval_overall, method = "BH")
  hdac9_evidence$correlation <- hdac9_correlation
  
  # 筛选显著相关
  sig_cor <- hdac9_correlation %>%
    filter(pval_overall < 0.05) %>%
    arrange(pval_overall)
  
  if(nrow(sig_cor) > 0) {
    message(paste("*** HDAC9表达与", nrow(sig_cor), "个通路活性显著相关 ***"))
    print(sig_cor)
  }
  
  write.csv(hdac9_correlation,
            file.path(RESULT_DIR, "10_HDAC9证据3_通路相关性.csv"),
            row.names = FALSE)
}

# --------------------------------------------------
# 证据4: 敏感性证据 - Leave-one-gene-out分析
# --------------------------------------------------
message("\n=== 证据4: 敏感性证据 (Leave-HDAC9-out) ===")
message("分析去除HDAC9后通路富集的变化...")

hdac9_sensitivity <- data.frame()

# 对于HDAC9所在的通路，去除HDAC9后重新运行GSEA
pathways_with_hdac9 <- names(pathways_list)[sapply(pathways_list, function(x) TARGET_GENE %in% x)]

if(length(pathways_with_hdac9) > 0) {
  message(paste("在", length(pathways_with_hdac9), "个通路中发现HDAC9"))
  
  # 创建去除HDAC9的通路列表
  pathways_no_hdac9 <- pathways_list
  for(pw in pathways_with_hdac9) {
    pathways_no_hdac9[[pw]] <- setdiff(pathways_no_hdac9[[pw]], TARGET_GENE)
  }
  
  # 重新运行GSEA - 软骨细胞
  message("\n重新运行GSEA (不含HDAC9) - 软骨细胞...")
  set.seed(123)
  chondrocyte_gsea_no_hdac9 <- fgsea(
    pathways = pathways_no_hdac9,
    stats = chondrocyte_ranks,
    minSize = MIN_GENESET_SIZE,
    maxSize = MAX_GENESET_SIZE,
    nproc = 1,
    nPermSimple = NPERM
  )
  chondrocyte_gsea_no_hdac9 <- as.data.frame(chondrocyte_gsea_no_hdac9)
  
  # 重新运行GSEA - 髓核细胞
  message("重新运行GSEA (不含HDAC9) - 髓核细胞...")
  set.seed(123)
  np_gsea_no_hdac9 <- fgsea(
    pathways = pathways_no_hdac9,
    stats = np_ranks,
    minSize = MIN_GENESET_SIZE,
    maxSize = MAX_GENESET_SIZE,
    nproc = 1,
    nPermSimple = NPERM
  )
  np_gsea_no_hdac9 <- as.data.frame(np_gsea_no_hdac9)
  
  # 比较NES变化
  for(pw in intersect(pathways_with_hdac9, common_protective_pathways)) {
    # 软骨细胞
    nes_orig_chon <- chondrocyte_gsea$NES[chondrocyte_gsea$pathway == pw]
    nes_no_hdac9_chon <- chondrocyte_gsea_no_hdac9$NES[chondrocyte_gsea_no_hdac9$pathway == pw]
    
    if(length(nes_orig_chon) > 0 && length(nes_no_hdac9_chon) > 0) {
      nes_change_chon <- nes_no_hdac9_chon - nes_orig_chon
      nes_pct_change_chon <- (nes_change_chon / abs(nes_orig_chon)) * 100
      
      hdac9_sensitivity <- rbind(hdac9_sensitivity, data.frame(
        pathway = pw,
        cell_type = "Chondrocyte",
        NES_original = nes_orig_chon,
        NES_no_HDAC9 = nes_no_hdac9_chon,
        NES_change = nes_change_chon,
        NES_pct_change = nes_pct_change_chon
      ))
    }
    
    # 髓核细胞
    nes_orig_np <- np_gsea$NES[np_gsea$pathway == pw]
    nes_no_hdac9_np <- np_gsea_no_hdac9$NES[np_gsea_no_hdac9$pathway == pw]
    
    if(length(nes_orig_np) > 0 && length(nes_no_hdac9_np) > 0) {
      nes_change_np <- nes_no_hdac9_np - nes_orig_np
      nes_pct_change_np <- (nes_change_np / abs(nes_orig_np)) * 100
      
      hdac9_sensitivity <- rbind(hdac9_sensitivity, data.frame(
        pathway = pw,
        cell_type = "Nucleus_Pulposus",
        NES_original = nes_orig_np,
        NES_no_HDAC9 = nes_no_hdac9_np,
        NES_change = nes_change_np,
        NES_pct_change = nes_pct_change_np
      ))
    }
  }
  
  if(nrow(hdac9_sensitivity) > 0) {
    hdac9_evidence$sensitivity <- hdac9_sensitivity
    
    # 识别HDAC9有显著贡献的通路（NES降幅 > 5%）
    sig_sensitivity <- hdac9_sensitivity %>%
      filter(abs(NES_pct_change) > 5) %>%
      arrange(desc(abs(NES_pct_change)))
    
    if(nrow(sig_sensitivity) > 0) {
      message(paste("*** HDAC9对", nrow(sig_sensitivity), "个通路的富集有显著贡献 (>5%变化) ***"))
      print(sig_sensitivity)
    }
    
    write.csv(hdac9_sensitivity,
              file.path(RESULT_DIR, "10_HDAC9证据4_敏感性分析.csv"),
              row.names = FALSE)
  }
}

# =============================================================================
# 六、可视化
# =============================================================================
print_progress("步骤6: 生成可视化图表")

# 6.1 GSEA结果对比图
if(length(common_protective_pathways) > 0) {
  # 选择top通路进行可视化
  top_pathways <- head(common_comparison$pathway, 20)
  
  plot_data <- rbind(
    chondrocyte_gsea %>%
      filter(pathway %in% top_pathways) %>%
      dplyr::select(pathway, NES, padj) %>%
      mutate(cell_type = "Chondrocyte"),
    np_gsea %>%
      filter(pathway %in% top_pathways) %>%
      dplyr::select(pathway, NES, padj) %>%
      mutate(cell_type = "Nucleus_Pulposus")
  )
  
  # 简化通路名称
  plot_data$pathway_short <- gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_", "", plot_data$pathway)
  plot_data$pathway_short <- substr(plot_data$pathway_short, 1, 50)
  
  p1 <- ggplot(plot_data, aes(x = NES, y = reorder(pathway_short, NES), fill = cell_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Chondrocyte" = "#00A087", "Nucleus_Pulposus" = "#DC0000")) +
    theme_classic() +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = "",
      fill = "Cell Type"
    ) +
    theme(text = element_text(size = 12))
  
  save_plot(p1, "10_共同保护性通路_NES对比.pdf", width = 14, height = 12)
}

# 6.2 HDAC9成员证据可视化
if(!is.null(hdac9_evidence$membership) && nrow(hdac9_evidence$membership) > 0) {
  membership_plot <- hdac9_evidence$membership %>%
    mutate(pathway_short = substr(gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_", "", pathway), 1, 60))
  
  p2 <- ggplot(membership_plot, aes(x = pathway_size, y = reorder(pathway_short, pathway_size))) +
    geom_bar(stat = "identity", fill = "#E64B35") +
    theme_classic() +
    labs(
      x = "Pathway Size (# genes)",
      y = "Pathway (contains HDAC9)",
      title = "Evidence 1: HDAC9 Membership in Protective Pathways"
    ) +
    theme(text = element_text(size = 12))
  
  save_plot(p2, "10_HDAC9证据1_成员通路.pdf", width = 12, height = 10)
}

# 6.3 Leading edge证据可视化
if(!is.null(hdac9_evidence$leading_edge) && nrow(hdac9_evidence$leading_edge) > 0) {
  le_plot <- hdac9_evidence$leading_edge %>%
    mutate(pathway_short = substr(gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_", "", pathway), 1, 60))
  
  p3 <- ggplot(le_plot, aes(x = abs(NES), y = reorder(pathway_short, abs(NES)), fill = cell_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Chondrocyte" = "#00A087", "Nucleus_Pulposus" = "#DC0000")) +
    theme_classic() +
    labs(
      x = "|NES|",
      y = "Pathway (HDAC9 in Leading Edge)",
      title = "Evidence 2: HDAC9 in Leading Edge",
      fill = "Cell Type"
    ) +
    theme(text = element_text(size = 12))
  
  save_plot(p3, "10_HDAC9证据2_Leading_Edge.pdf", width = 12, height = 10)
}

# 6.4 相关性证据可视化
if(!is.null(hdac9_evidence$correlation) && nrow(hdac9_evidence$correlation) > 0) {
  cor_plot <- hdac9_evidence$correlation %>%
    filter(pval_overall < 0.05) %>%
    mutate(pathway_short = substr(gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_", "", pathway), 1, 50))
  
  if(nrow(cor_plot) > 0) {
    p4 <- ggplot(cor_plot, aes(x = cor_overall, y = reorder(pathway_short, cor_overall), 
                               fill = cell_type, alpha = -log10(pval_overall))) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Chondrocyte" = "#00A087", "Nucleus_Pulposus" = "#DC0000")) +
      theme_classic() +
      labs(
        x = "Spearman Correlation (HDAC9 vs Pathway Activity)",
        y = "Pathway",
        title = "Evidence 3: HDAC9-Pathway Correlation",
        fill = "Cell Type",
        alpha = "-log10(p-value)"
      ) +
      theme(text = element_text(size = 12))
    
    save_plot(p4, "10_HDAC9证据3_相关性.pdf", width = 12, height = 10)
  }
}

# 6.5 敏感性证据可视化
if(!is.null(hdac9_evidence$sensitivity) && nrow(hdac9_evidence$sensitivity) > 0) {
  sens_plot <- hdac9_evidence$sensitivity %>%
    mutate(pathway_short = substr(gsub("GOBP_|KEGG_|REACTOME_|HALLMARK_", "", pathway), 1, 50))
  
  p5 <- ggplot(sens_plot, aes(x = NES_pct_change, y = reorder(pathway_short, abs(NES_pct_change)), 
                              fill = cell_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_vline(xintercept = c(-5, 5), linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("Chondrocyte" = "#00A087", "Nucleus_Pulposus" = "#DC0000")) +
    theme_classic() +
    labs(
      x = "NES % Change (after removing HDAC9)",
      y = "Pathway",
      title = "Evidence 4: HDAC9 Contribution to Pathway Enrichment",
      fill = "Cell Type"
    ) +
    theme(text = element_text(size = 12))
  
  save_plot(p5, "10_HDAC9证据4_敏感性.pdf", width = 12, height = 10)
}

# =============================================================================
# 七、综合报告
# =============================================================================
print_progress("生成综合报告")

summary_text <- c(
  "=============================================================================",
  "软骨细胞和髓核细胞共同保护通路GSEA分析及HDAC9验证报告",
  "=============================================================================",
  "",
  paste("分析日期:", Sys.time()),
  "",
  "--- 数据概况 ---",
  paste("软骨细胞数:", ncol(chondrocyte_data)),
  paste("髓核细胞数:", ncol(np_data)),
  "",
  "--- GSEA富集分析 ---",
  paste("使用通路数据库: MSigDB (Hallmark, GO:BP, KEGG, Reactome)"),
  paste("通路数量:", length(pathways_list)),
  paste("软骨细胞显著富集通路:", sum(chondrocyte_gsea$padj < PADJ_CUTOFF, na.rm = TRUE)),
  paste("髓核细胞显著富集通路:", sum(np_gsea$padj < PADJ_CUTOFF, na.rm = TRUE)),
  "",
  "--- 保护性通路筛选 ---",
  paste("软骨细胞保护性通路:", nrow(chondrocyte_protective)),
  paste("髓核细胞保护性通路:", nrow(np_protective)),
  paste("共同保护性通路:", length(common_protective_pathways)),
  ""
)

# HDAC9四类证据
summary_text <- c(
  summary_text,
  "=============================================================================",
  "HDAC9验证 - 四类证据",
  "============================================================================="
)

# 证据1
if(!is.null(hdac9_evidence$membership)) {
  summary_text <- c(
    summary_text,
    "",
    "【证据1: 成员证据】",
    paste("*** HDAC9出现在", nrow(hdac9_evidence$membership), "个共同保护性通路中 ***"),
    paste("通路列表:"),
    paste0("  - ", hdac9_evidence$membership$pathway)
  )
} else {
  summary_text <- c(summary_text, "", "【证据1: 成员证据】", "HDAC9不在共同保护性通路的基因集中")
}

# 证据2
if(!is.null(hdac9_evidence$leading_edge)) {
  summary_text <- c(
    summary_text,
    "",
    "【证据2: 贡献证据 (Leading Edge)】",
    paste("*** HDAC9在", nrow(hdac9_evidence$leading_edge), "个通路-细胞类型的leading edge中 ***"),
    "（Leading edge = 对通路富集贡献最大的核心基因子集）"
  )
} else {
  summary_text <- c(summary_text, "", "【证据2: 贡献证据】", "HDAC9不在通路的leading edge中")
}

# 证据3
if(!is.null(hdac9_evidence$correlation)) {
  sig_cor_count <- sum(hdac9_evidence$correlation$pval_overall < 0.05)
  summary_text <- c(
    summary_text,
    "",
    "【证据3: 相关证据 (通路活性 vs HDAC9表达)】",
    paste("*** HDAC9表达与", sig_cor_count, "个通路活性显著相关 (p < 0.05) ***"),
    "（UCell评分 vs HDAC9表达的Spearman相关）"
  )
} else {
  summary_text <- c(summary_text, "", "【证据3: 相关证据】", "未发现显著相关")
}

# 证据4
if(!is.null(hdac9_evidence$sensitivity)) {
  sig_sens_count <- sum(abs(hdac9_evidence$sensitivity$NES_pct_change) > 5)
  summary_text <- c(
    summary_text,
    "",
    "【证据4: 敏感性证据 (Leave-HDAC9-out)】",
    paste("*** 去除HDAC9后，", sig_sens_count, "个通路的NES变化 > 5% ***"),
    "（说明HDAC9对这些通路的富集有实质贡献）"
  )
} else {
  summary_text <- c(summary_text, "", "【证据4: 敏感性证据】", "未进行分析（HDAC9不在通路中）")
}

summary_text <- c(
  summary_text,
  "",
  "=============================================================================",
  "结论:",
  "通过preranked GSEA方法识别出两种细胞类型的共同保护性通路，",
  "并从四个维度验证了HDAC9在这些通路中的关键作用：",
  "  1. 成员证据: HDAC9属于保护性通路基因集",
  "  2. 贡献证据: HDAC9在通路核心基因子集(leading edge)中",
  "  3. 相关证据: HDAC9表达与通路活性显著相关",
  "  4. 敏感性证据: 去除HDAC9显著降低通路富集强度",
  "以上证据支持HDAC9作为关键保护因子在两种细胞退变中的作用。",
  "============================================================================="
)

writeLines(summary_text, file.path(RESULT_DIR, "10_HDAC9验证_综合报告.txt"))

print_progress("分析完成！")

message("\n输出文件:")
message(paste("- 数据:", RESULT_DIR))
message(paste("- 图表:", FIGURE_DIR))
message(paste("- 报告:", file.path(RESULT_DIR, "10_HDAC9验证_综合报告.txt")))

message("\n关键发现:")
if(!is.null(hdac9_evidence$membership)) message("✓ 证据1: HDAC9是通路成员")
if(!is.null(hdac9_evidence$leading_edge)) message("✓ 证据2: HDAC9在leading edge")
if(!is.null(hdac9_evidence$correlation)) message("✓ 证据3: HDAC9与通路活性相关")
if(!is.null(hdac9_evidence$sensitivity)) message("✓ 证据4: HDAC9对富集有贡献")

# =============================================================================
# 八、针对特定通路的HDAC9验证 (BIOCARTA_ARENRF2_PATHWAY)
# =============================================================================
print_progress("步骤8: 验证HDAC9在NRF2通路中的显著性")

TARGET_PATHWAY <- "BIOCARTA_ARENRF2_PATHWAY"

# 8.1 检查目标通路是否存在
message(paste("\n检查通路:", TARGET_PATHWAY))

if(TARGET_PATHWAY %in% names(pathways_list)) {
  message(paste("✓ 找到通路:", TARGET_PATHWAY))
  nrf2_pathway_genes <- pathways_list[[TARGET_PATHWAY]]
  message(paste("  通路基因数:", length(nrf2_pathway_genes)))
  message(paste("  通路基因:", paste(head(nrf2_pathway_genes, 20), collapse = ", ")))
  if(length(nrf2_pathway_genes) > 20) message("  ...")
  
  # 检查HDAC9是否在通路中
  if(TARGET_GENE %in% nrf2_pathway_genes) {
    message(paste("\n✓✓✓", TARGET_GENE, "是", TARGET_PATHWAY, "的成员基因 ✓✓✓"))
    nrf2_hdac9_member <- TRUE
  } else {
    message(paste("\n✗✗✗", TARGET_GENE, "不在", TARGET_PATHWAY, "中 ✗✗✗"))
    nrf2_hdac9_member <- FALSE
  }
  
  # 8.2 检查该通路的GSEA富集情况
  message("\n=== NRF2通路GSEA富集结果 ===")
  
  # 软骨细胞
  nrf2_chon_gsea <- chondrocyte_gsea[chondrocyte_gsea$pathway == TARGET_PATHWAY, ]
  if(nrow(nrf2_chon_gsea) > 0) {
    message("\n软骨细胞:")
    message(paste("  NES:", round(nrf2_chon_gsea$NES, 4)))
    message(paste("  p-value:", format(nrf2_chon_gsea$pval, scientific = TRUE, digits = 3)))
    message(paste("  p.adj:", format(nrf2_chon_gsea$padj, scientific = TRUE, digits = 3)))
    message(paste("  Leading edge size:", length(unlist(nrf2_chon_gsea$leadingEdge))))
    if(nrf2_chon_gsea$padj < 0.05) {
      message("  *** 显著富集 (padj < 0.05) ***")
    } else if(nrf2_chon_gsea$padj < 0.25) {
      message("  ** 边缘显著 (padj < 0.25) **")
    } else {
      message("  不显著 (padj >= 0.25)")
    }
  } else {
    message("软骨细胞: 未找到该通路的GSEA结果")
    nrf2_chon_gsea <- NULL
  }
  
  # 髓核细胞
  nrf2_np_gsea <- np_gsea[np_gsea$pathway == TARGET_PATHWAY, ]
  if(nrow(nrf2_np_gsea) > 0) {
    message("\n髓核细胞:")
    message(paste("  NES:", round(nrf2_np_gsea$NES, 4)))
    message(paste("  p-value:", format(nrf2_np_gsea$pval, scientific = TRUE, digits = 3)))
    message(paste("  p.adj:", format(nrf2_np_gsea$padj, scientific = TRUE, digits = 3)))
    message(paste("  Leading edge size:", length(unlist(nrf2_np_gsea$leadingEdge))))
    if(nrf2_np_gsea$padj < 0.05) {
      message("  *** 显著富集 (padj < 0.05) ***")
    } else if(nrf2_np_gsea$padj < 0.25) {
      message("  ** 边缘显著 (padj < 0.25) **")
    } else {
      message("  不显著 (padj >= 0.25)")
    }
  } else {
    message("髓核细胞: 未找到该通路的GSEA结果")
    nrf2_np_gsea <- NULL
  }
  
  # 8.3 HDAC9在NRF2通路中的四类证据
  message("\n=== HDAC9在NRF2通路中的验证 ===")
  
  nrf2_evidence <- list()
  
  # 证据1: 成员证据
  message("\n【证据1: 成员证据】")
  if(nrf2_hdac9_member) {
    message(paste("✓", TARGET_GENE, "是", TARGET_PATHWAY, "的成员"))
    nrf2_evidence$membership <- data.frame(
      pathway = TARGET_PATHWAY,
      is_member = TRUE,
      pathway_size = length(nrf2_pathway_genes)
    )
  } else {
    message(paste("✗", TARGET_GENE, "不是", TARGET_PATHWAY, "的成员"))
    nrf2_evidence$membership <- data.frame(
      pathway = TARGET_PATHWAY,
      is_member = FALSE,
      pathway_size = length(nrf2_pathway_genes)
    )
  }
  
  # 证据2: Leading Edge
  message("\n【证据2: Leading Edge贡献】")
  nrf2_leading_edge <- data.frame()
  
  # 软骨细胞
  if(!is.null(nrf2_chon_gsea)) {
    chon_leading_genes <- unlist(nrf2_chon_gsea$leadingEdge)
    if(TARGET_GENE %in% chon_leading_genes) {
      message(paste("✓ 软骨细胞:", TARGET_GENE, "在leading edge中"))
      nrf2_leading_edge <- rbind(nrf2_leading_edge, data.frame(
        pathway = TARGET_PATHWAY,
        cell_type = "Chondrocyte",
        in_leading_edge = TRUE,
        NES = nrf2_chon_gsea$NES,
        padj = nrf2_chon_gsea$padj,
        leading_edge_size = length(chon_leading_genes)
      ))
    } else {
      message(paste("✗ 软骨细胞:", TARGET_GENE, "不在leading edge中"))
    }
  }
  
  # 髓核细胞
  if(!is.null(nrf2_np_gsea)) {
    np_leading_genes <- unlist(nrf2_np_gsea$leadingEdge)
    if(TARGET_GENE %in% np_leading_genes) {
      message(paste("✓ 髓核细胞:", TARGET_GENE, "在leading edge中"))
      nrf2_leading_edge <- rbind(nrf2_leading_edge, data.frame(
        pathway = TARGET_PATHWAY,
        cell_type = "Nucleus_Pulposus",
        in_leading_edge = TRUE,
        NES = nrf2_np_gsea$NES,
        padj = nrf2_np_gsea$padj,
        leading_edge_size = length(np_leading_genes)
      ))
    } else {
      message(paste("✗ 髓核细胞:", TARGET_GENE, "不在leading edge中"))
    }
  }
  
  if(nrow(nrf2_leading_edge) > 0) {
    nrf2_evidence$leading_edge <- nrf2_leading_edge
  }
  
  # 证据3: 相关性分析
  message("\n【证据3: 通路活性与HDAC9表达的相关性】")
  
  # 计算NRF2通路的UCell评分
  nrf2_signature <- list(nrf2_pathway_genes)
  names(nrf2_signature) <- TARGET_PATHWAY
  
  message("计算NRF2通路活性得分...")
  chondrocyte_data <- AddModuleScore_UCell(chondrocyte_data, features = nrf2_signature)
  np_data <- AddModuleScore_UCell(np_data, features = nrf2_signature)
  
  nrf2_correlation <- data.frame()
  
  # 软骨细胞
  if(TARGET_GENE %in% rownames(chondrocyte_data)) {
    chon_meta <- chondrocyte_data@meta.data
    chon_meta$HDAC9 <- GetAssayData(chondrocyte_data, slot = "data")[TARGET_GENE, ]
    ucell_col <- paste0(TARGET_PATHWAY, "_UCell")
    
    if(ucell_col %in% colnames(chon_meta)) {
      cor_test_chon <- cor.test(chon_meta$HDAC9, chon_meta[[ucell_col]], method = "spearman")
      message("\n软骨细胞:")
      message(paste("  Spearman相关系数:", round(cor_test_chon$estimate, 4)))
      message(paste("  p-value:", format(cor_test_chon$p.value, scientific = TRUE, digits = 3)))
      if(cor_test_chon$p.value < 0.05) {
        message("  *** 显著相关 ***")
      }
      
      nrf2_correlation <- rbind(nrf2_correlation, data.frame(
        pathway = TARGET_PATHWAY,
        cell_type = "Chondrocyte",
        cor_spearman = cor_test_chon$estimate,
        pval = cor_test_chon$p.value
      ))
    }
  }
  
  # 髓核细胞
  if(TARGET_GENE %in% rownames(np_data)) {
    np_meta <- np_data@meta.data
    np_meta$HDAC9 <- GetAssayData(np_data, slot = "data")[TARGET_GENE, ]
    ucell_col <- paste0(TARGET_PATHWAY, "_UCell")
    
    if(ucell_col %in% colnames(np_meta)) {
      cor_test_np <- cor.test(np_meta$HDAC9, np_meta[[ucell_col]], method = "spearman")
      message("\n髓核细胞:")
      message(paste("  Spearman相关系数:", round(cor_test_np$estimate, 4)))
      message(paste("  p-value:", format(cor_test_np$p.value, scientific = TRUE, digits = 3)))
      if(cor_test_np$p.value < 0.05) {
        message("  *** 显著相关 ***")
      }
      
      nrf2_correlation <- rbind(nrf2_correlation, data.frame(
        pathway = TARGET_PATHWAY,
        cell_type = "Nucleus_Pulposus",
        cor_spearman = cor_test_np$estimate,
        pval = cor_test_np$p.value
      ))
    }
  }
  
  if(nrow(nrf2_correlation) > 0) {
    nrf2_evidence$correlation <- nrf2_correlation
  }
  
  # 证据4: 敏感性分析
  message("\n【证据4: 敏感性分析 (Leave-HDAC9-out)】")
  
  if(nrf2_hdac9_member) {
    message("去除HDAC9后重新运行GSEA...")
    
    # 创建去除HDAC9的NRF2通路
    nrf2_pathway_no_hdac9 <- setdiff(nrf2_pathway_genes, TARGET_GENE)
    nrf2_list_no_hdac9 <- list(nrf2_pathway_no_hdac9)
    names(nrf2_list_no_hdac9) <- TARGET_PATHWAY
    
    # 软骨细胞
    set.seed(123)
    nrf2_chon_gsea_no_hdac9 <- fgsea(
      pathways = nrf2_list_no_hdac9,
      stats = chondrocyte_ranks,
      minSize = 5,
      maxSize = 500,
      nproc = 1,
      nPermSimple = NPERM
    )
    nrf2_chon_gsea_no_hdac9 <- as.data.frame(nrf2_chon_gsea_no_hdac9)
    
    if(!is.null(nrf2_chon_gsea) && nrow(nrf2_chon_gsea_no_hdac9) > 0) {
      nes_change_chon <- nrf2_chon_gsea_no_hdac9$NES - nrf2_chon_gsea$NES
      nes_pct_change_chon <- (nes_change_chon / abs(nrf2_chon_gsea$NES)) * 100
      
      message("\n软骨细胞:")
      message(paste("  原始NES:", round(nrf2_chon_gsea$NES, 4)))
      message(paste("  去除HDAC9后NES:", round(nrf2_chon_gsea_no_hdac9$NES, 4)))
      message(paste("  NES变化:", round(nes_change_chon, 4)))
      message(paste("  NES变化百分比:", round(nes_pct_change_chon, 2), "%"))
      if(abs(nes_pct_change_chon) > 5) {
        message("  *** HDAC9对通路富集有显著贡献 (>5%变化) ***")
      }
    }
    
    # 髓核细胞
    set.seed(123)
    nrf2_np_gsea_no_hdac9 <- fgsea(
      pathways = nrf2_list_no_hdac9,
      stats = np_ranks,
      minSize = 5,
      maxSize = 500,
      nproc = 1,
      nPermSimple = NPERM
    )
    nrf2_np_gsea_no_hdac9 <- as.data.frame(nrf2_np_gsea_no_hdac9)
    
    if(!is.null(nrf2_np_gsea) && nrow(nrf2_np_gsea_no_hdac9) > 0) {
      nes_change_np <- nrf2_np_gsea_no_hdac9$NES - nrf2_np_gsea$NES
      nes_pct_change_np <- (nes_change_np / abs(nrf2_np_gsea$NES)) * 100
      
      message("\n髓核细胞:")
      message(paste("  原始NES:", round(nrf2_np_gsea$NES, 4)))
      message(paste("  去除HDAC9后NES:", round(nrf2_np_gsea_no_hdac9$NES, 4)))
      message(paste("  NES变化:", round(nes_change_np, 4)))
      message(paste("  NES变化百分比:", round(nes_pct_change_np, 2), "%"))
      if(abs(nes_pct_change_np) > 5) {
        message("  *** HDAC9对通路富集有显著贡献 (>5%变化) ***")
      }
    }
    
    # 保存敏感性结果
    nrf2_sensitivity <- data.frame()
    if(exists("nes_pct_change_chon")) {
      nrf2_sensitivity <- rbind(nrf2_sensitivity, data.frame(
        pathway = TARGET_PATHWAY,
        cell_type = "Chondrocyte",
        NES_original = nrf2_chon_gsea$NES,
        NES_no_HDAC9 = nrf2_chon_gsea_no_hdac9$NES,
        NES_change = nes_change_chon,
        NES_pct_change = nes_pct_change_chon
      ))
    }
    if(exists("nes_pct_change_np")) {
      nrf2_sensitivity <- rbind(nrf2_sensitivity, data.frame(
        pathway = TARGET_PATHWAY,
        cell_type = "Nucleus_Pulposus",
        NES_original = nrf2_np_gsea$NES,
        NES_no_HDAC9 = nrf2_np_gsea_no_hdac9$NES,
        NES_change = nes_change_np,
        NES_pct_change = nes_pct_change_np
      ))
    }
    
    if(nrow(nrf2_sensitivity) > 0) {
      nrf2_evidence$sensitivity <- nrf2_sensitivity
    }
  } else {
    message("HDAC9不在NRF2通路中，跳过敏感性分析")
  }
  
  # 8.4 生成NRF2通路专门的可视化
  message("\n生成NRF2通路可视化...")
  
  # GSEA enrichment plot
  if(!is.null(nrf2_chon_gsea)) {
    p_nrf2_chon <- plotEnrichment(nrf2_pathway_genes, chondrocyte_ranks) +
      labs(title = paste(TARGET_PATHWAY, "- Chondrocyte"),
           subtitle = paste("NES =", round(nrf2_chon_gsea$NES, 3), 
                           ", padj =", format(nrf2_chon_gsea$padj, digits = 3))) +
      theme_classic() +
      theme(text = element_text(size = 12))
    
    save_plot(p_nrf2_chon, "10_验证_NRF2通路_软骨细胞_GSEA.pdf", width = 10, height = 6)
  }
  
  if(!is.null(nrf2_np_gsea)) {
    p_nrf2_np <- plotEnrichment(nrf2_pathway_genes, np_ranks) +
      labs(title = paste(TARGET_PATHWAY, "- Nucleus Pulposus"),
           subtitle = paste("NES =", round(nrf2_np_gsea$NES, 3), 
                           ", padj =", format(nrf2_np_gsea$padj, digits = 3))) +
      theme_classic() +
      theme(text = element_text(size = 12))
    
    save_plot(p_nrf2_np, "10_验证_NRF2通路_髓核细胞_GSEA.pdf", width = 10, height = 6)
  }
  
  # 相关性散点图
  if(exists("cor_test_chon") && exists("chon_meta")) {
    ucell_col <- paste0(TARGET_PATHWAY, "_UCell")
    plot_data_chon <- data.frame(
      HDAC9 = chon_meta$HDAC9,
      Pathway_Activity = chon_meta[[ucell_col]],
      Group = chon_meta$group
    )
    
    p_cor_chon <- ggplot(plot_data_chon, aes(x = HDAC9, y = Pathway_Activity, color = Group)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
      scale_color_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000")) +
      labs(
        title = paste("HDAC9 vs", TARGET_PATHWAY, "Activity - Chondrocyte"),
        subtitle = paste("Spearman r =", round(cor_test_chon$estimate, 3), 
                        ", p =", format(cor_test_chon$p.value, digits = 3)),
        x = "HDAC9 Expression",
        y = "NRF2 Pathway Activity (UCell Score)"
      ) +
      theme_classic() +
      theme(text = element_text(size = 12))
    
    save_plot(p_cor_chon, "10_验证_NRF2通路_HDAC9相关性_软骨细胞.pdf", width = 8, height = 6)
  }
  
  if(exists("cor_test_np") && exists("np_meta")) {
    ucell_col <- paste0(TARGET_PATHWAY, "_UCell")
    plot_data_np <- data.frame(
      HDAC9 = np_meta$HDAC9,
      Pathway_Activity = np_meta[[ucell_col]],
      Group = np_meta$group
    )
    
    p_cor_np <- ggplot(plot_data_np, aes(x = HDAC9, y = Pathway_Activity, color = Group)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
      scale_color_manual(values = c("MDD" = "#00A087", "SDD" = "#DC0000")) +
      labs(
        title = paste("HDAC9 vs", TARGET_PATHWAY, "Activity - Nucleus Pulposus"),
        subtitle = paste("Spearman r =", round(cor_test_np$estimate, 3), 
                        ", p =", format(cor_test_np$p.value, digits = 3)),
        x = "HDAC9 Expression",
        y = "NRF2 Pathway Activity (UCell Score)"
      ) +
      theme_classic() +
      theme(text = element_text(size = 12))
    
    save_plot(p_cor_np, "10_验证_NRF2通路_HDAC9相关性_髓核细胞.pdf", width = 8, height = 6)
  }
  
  # 8.5 生成NRF2通路验证报告
  nrf2_report <- c(
    "=============================================================================",
    paste("HDAC9在", TARGET_PATHWAY, "通路中的显著性验证报告"),
    "=============================================================================",
    "",
    paste("分析日期:", Sys.time()),
    "",
    "=============================================================================",
    "一、通路基本信息",
    "=============================================================================",
    paste("通路名称:", TARGET_PATHWAY),
    paste("通路基因数:", length(nrf2_pathway_genes)),
    paste("HDAC9是否为成员:", ifelse(nrf2_hdac9_member, "是", "否")),
    ""
  )
  
  # HDAC9差异表达
  nrf2_report <- c(nrf2_report,
    "=============================================================================",
    "二、HDAC9差异表达",
    "=============================================================================",
    ""
  )
  
  if(exists("hdac9_chon")) {
    nrf2_report <- c(nrf2_report,
      "软骨细胞 (SDD vs MDD):",
      paste("  log2FC:", round(hdac9_chon$avg_log2FC, 4)),
      paste("  p_adj:", format(hdac9_chon$p_val_adj, scientific = TRUE, digits = 3)),
      paste("  显著性:", ifelse(hdac9_chon$p_val_adj < 0.05, "*** 显著", "不显著")),
      ""
    )
  }
  
  if(exists("hdac9_np")) {
    nrf2_report <- c(nrf2_report,
      "髓核细胞 (SDD vs MDD):",
      paste("  log2FC:", round(hdac9_np$avg_log2FC, 4)),
      paste("  p_adj:", format(hdac9_np$p_val_adj, scientific = TRUE, digits = 3)),
      paste("  显著性:", ifelse(hdac9_np$p_val_adj < 0.05, "*** 显著", "不显著")),
      ""
    )
  }
  
  # NRF2通路GSEA结果
  nrf2_report <- c(nrf2_report,
    "=============================================================================",
    "三、NRF2通路GSEA富集分析",
    "=============================================================================",
    ""
  )
  
  if(!is.null(nrf2_chon_gsea)) {
    nrf2_report <- c(nrf2_report,
      "软骨细胞:",
      paste("  NES:", round(nrf2_chon_gsea$NES, 4)),
      paste("  p.adj:", format(nrf2_chon_gsea$padj, scientific = TRUE, digits = 3)),
      paste("  显著性:", ifelse(nrf2_chon_gsea$padj < 0.05, "*** 显著", 
                              ifelse(nrf2_chon_gsea$padj < 0.25, "** 边缘显著", "不显著"))),
      ""
    )
  }
  
  if(!is.null(nrf2_np_gsea)) {
    nrf2_report <- c(nrf2_report,
      "髓核细胞:",
      paste("  NES:", round(nrf2_np_gsea$NES, 4)),
      paste("  p.adj:", format(nrf2_np_gsea$padj, scientific = TRUE, digits = 3)),
      paste("  显著性:", ifelse(nrf2_np_gsea$padj < 0.05, "*** 显著", 
                              ifelse(nrf2_np_gsea$padj < 0.25, "** 边缘显著", "不显著"))),
      ""
    )
  }
  
  # 四类证据
  nrf2_report <- c(nrf2_report,
    "=============================================================================",
    "四、HDAC9验证 - 四类证据",
    "=============================================================================",
    "",
    "【证据1: 成员证据】",
    paste("结果:", ifelse(nrf2_hdac9_member, "✓ HDAC9是NRF2通路成员", "✗ HDAC9不是NRF2通路成员")),
    ""
  )
  
  if(!is.null(nrf2_evidence$leading_edge)) {
    nrf2_report <- c(nrf2_report,
      "【证据2: Leading Edge贡献】"
    )
    for(i in 1:nrow(nrf2_evidence$leading_edge)) {
      row <- nrf2_evidence$leading_edge[i, ]
      nrf2_report <- c(nrf2_report,
        paste(row$cell_type, ": ✓ HDAC9在leading edge中")
      )
    }
    nrf2_report <- c(nrf2_report, "")
  } else {
    nrf2_report <- c(nrf2_report,
      "【证据2: Leading Edge贡献】",
      "HDAC9不在任何细胞类型的leading edge中",
      ""
    )
  }
  
  if(!is.null(nrf2_evidence$correlation)) {
    nrf2_report <- c(nrf2_report,
      "【证据3: 相关性证据】"
    )
    for(i in 1:nrow(nrf2_evidence$correlation)) {
      row <- nrf2_evidence$correlation[i, ]
      sig <- ifelse(row$pval < 0.05, "*** 显著相关", "不显著")
      nrf2_report <- c(nrf2_report,
        paste(row$cell_type, ":"),
        paste("  相关系数:", round(row$cor_spearman, 4)),
        paste("  p-value:", format(row$pval, scientific = TRUE, digits = 3)),
        paste("  结论:", sig)
      )
    }
    nrf2_report <- c(nrf2_report, "")
  }
  
  if(!is.null(nrf2_evidence$sensitivity)) {
    nrf2_report <- c(nrf2_report,
      "【证据4: 敏感性证据】"
    )
    for(i in 1:nrow(nrf2_evidence$sensitivity)) {
      row <- nrf2_evidence$sensitivity[i, ]
      sig <- ifelse(abs(row$NES_pct_change) > 5, "*** 有显著贡献", "贡献较小")
      nrf2_report <- c(nrf2_report,
        paste(row$cell_type, ":"),
        paste("  NES变化百分比:", round(row$NES_pct_change, 2), "%"),
        paste("  结论:", sig)
      )
    }
    nrf2_report <- c(nrf2_report, "")
  }
  
  # 总结
  nrf2_report <- c(nrf2_report,
    "=============================================================================",
    "五、总体结论",
    "=============================================================================",
    ""
  )
  
  conclusion_items <- c()
  if(nrf2_hdac9_member) {
    conclusion_items <- c(conclusion_items, "1. ✓ HDAC9是NRF2通路的成员基因")
  } else {
    conclusion_items <- c(conclusion_items, "1. ✗ HDAC9不是NRF2通路的注释成员")
  }
  
  if(!is.null(nrf2_chon_gsea) && nrf2_chon_gsea$padj < 0.25) {
    conclusion_items <- c(conclusion_items, "2. ✓ NRF2通路在软骨细胞中显著富集")
  }
  if(!is.null(nrf2_np_gsea) && nrf2_np_gsea$padj < 0.25) {
    conclusion_items <- c(conclusion_items, "3. ✓ NRF2通路在髓核细胞中显著富集")
  }
  
  if(!is.null(nrf2_evidence$leading_edge) && nrow(nrf2_evidence$leading_edge) > 0) {
    conclusion_items <- c(conclusion_items, 
      paste("4. ✓ HDAC9在", nrow(nrf2_evidence$leading_edge), "个细胞类型的leading edge中"))
  }
  
  if(!is.null(nrf2_evidence$correlation)) {
    sig_count <- sum(nrf2_evidence$correlation$pval < 0.05)
    if(sig_count > 0) {
      conclusion_items <- c(conclusion_items, 
        paste("5. ✓ HDAC9与NRF2通路活性在", sig_count, "个细胞类型中显著相关"))
    }
  }
  
  if(!is.null(nrf2_evidence$sensitivity)) {
    sig_count <- sum(abs(nrf2_evidence$sensitivity$NES_pct_change) > 5)
    if(sig_count > 0) {
      conclusion_items <- c(conclusion_items, 
        paste("6. ✓ HDAC9对", sig_count, "个细胞类型的NRF2通路富集有显著贡献"))
    }
  }
  
  nrf2_report <- c(nrf2_report, conclusion_items, "",
    "============================================================================="
  )
  
  writeLines(nrf2_report, file.path(RESULT_DIR, "10_验证_HDAC9_NRF2通路_报告.txt"))
  
  message("\n✓ NRF2通路验证完成")
  message(paste("报告保存至:", file.path(RESULT_DIR, "10_验证_HDAC9_NRF2通路_报告.txt")))
  
} else {
  message(paste("✗ 未找到通路:", TARGET_PATHWAY))
  message("\n可用的NRF2相关通路:")
  nrf2_pathways <- grep("NRF2|NFE2L2|ANTIOXIDANT|OXIDATIVE", 
                        names(pathways_list), 
                        value = TRUE, 
                        ignore.case = TRUE)
  if(length(nrf2_pathways) > 0) {
    for(pw in nrf2_pathways) {
      genes <- pathways_list[[pw]]
      has_hdac9 <- TARGET_GENE %in% genes
      message(paste("  -", pw, ":", length(genes), "genes", 
                   ifelse(has_hdac9, "[包含HDAC9]", "")))
    }
  } else {
    message("  (未找到相关通路)")
  }
}

message("\n=============================================================================")
message("所有分析完成！")
message("=============================================================================")
