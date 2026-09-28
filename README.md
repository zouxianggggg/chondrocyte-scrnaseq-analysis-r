# Chondrocyte scRNA-seq analysis scripts

R scripts for single-cell RNA-seq processing and downstream analysis of chondrocytes and related intervertebral-disc cell populations. This repository archives the scripts supplied for the SDD versus MDD analysis.

## Contents

- `01_preprocessing_optimized.R` and `02_cell_annotation_optimized.R`: preprocessing, integration, clustering, and cell-type annotation.
- `03_group_comparison_SDD_vs_MDD.R`: cell composition and differential-expression analyses between the two study groups.
- `04_chondrocyte_subanalysis.R` through `09_cluster15_trajectory_analysis.R`: chondrocyte subclustering, trajectory, differential-expression, and cluster 15 analyses.
- `10_共同保护通路_HDAC9验证.R` and `10_共同保护通路_HDAC9验证_GSEA.R`: shared pathway and HDAC9 analyses involving chondrocytes and nucleus pulposus cells.
- `config.R` and `utils.R`: sample definitions, analysis settings, and shared functions.

The numbered files reflect the analysis stages; the two files numbered `10` are separate analyses.

## Inputs and use

The raw 10x matrices, intermediate R objects, and generated results are **not included**. `config.R` contains the original Windows project path and relative input locations for six samples; adapt these paths to your local environment before running the scripts. Run scripts from the directory containing `config.R` and `utils.R` so their relative `source()` calls resolve. Package requirements appear in the individual scripts and `config.R`.

## Archive status

Four supplied files (`03_group_comparison_SDD_vs_MDD.R` through `06_chondrocyte_cluster15_DEG.R`) contained invalid UTF-8 text and truncated Chinese strings. Their affected comments and display strings were repaired, including missing closing quotes. All 13 R files now pass a syntax-only `parse()` check. The workflow has not been executed here, so this archive should not be described as a fully reproduced analysis.
