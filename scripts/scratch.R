# Script for creating Fig. 2 in manuscript

LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(Seurat)
  library(tidyverse)
  library(pheatmap)
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')
pl_so = readRDS("analysis/final_rds/subset/pl_obsctr_final.rds")
DefaultAssay(pl_so) = "RNA"

FeaturePlot(pl_so, feature="Rxra", label=T, order=T)
DotPlot(pl_so, feature="Rxra", scale=F)

Idents(pl_so)
FeaturePlot(seurat_integrated_pl_rename, feature="Rxra", label=T, order=T)
DotPlot(seurat_integrated_pl_rename, feature="Rxra", scale=F)

