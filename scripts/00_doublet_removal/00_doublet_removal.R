#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# This script starts with raw cellranger files and removes doublets using DoubletFinder
# It takes three arguments: Rscript 00_doublet_removal.R sample_name cellranger_directory output_path

# HPC R library
LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

library(Seurat)
library(tidyverse)
library(sctransform)
library(DoubletFinder)

# test if there are three arguments
if (length(args)==0) {
  stop("Three argument must be supplied", call.=FALSE)
} 

sample = args[1]
cellranger_path=args[2]
outpath=args[3]

print(sample)

# run seurat workflow ----
so_data <- Read10X(data.dir = paste0(cellranger_path,sample,"/outs/filtered_feature_bc_matrix/"))
so <- CreateSeuratObject(counts = so_data, min.features =100)
so <- NormalizeData(so)
so <- SCTransform(so)
so <- RunPCA(so, verbose = FALSE)
so <- RunUMAP(so, dims = 1:50, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:50, verbose = FALSE)
so <- FindClusters(so, verbose = FALSE)

# run doublet finder  ----
# see vignette as an example https://github.com/chris-mcginnis-ucsf/DoubletFinder
sweep = paramSweep_v3(so, PCs = 1:50, sct = TRUE)
sweep.stats <- summarizeSweep(sweep, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
pK <- as.numeric(levels(pK))[pK]
annotations <- so@meta.data$SCT_snn_res.0.8
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.04*nrow(so@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
so <- doubletFinder_v3(so, PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

# add to seurat metadata
DF.name = colnames(so@meta.data)[grepl("DF.classification", colnames(so@meta.data))]
so@meta.data$DF.classification = so@meta.data[[DF.name]]

saveRDS(so, 
        paste0(outpath, sample,"_with_doublets.rds"))

# make plots ----
#so = readRDS("analysis/doublet_removal/pre_doublet_removal/AEb2-11_P48_CTR_PL_pre_doublet_removal.rds")

## make a plot of DF classification
Idents(so)<-"DF.classification"

p = DimPlot(so, label=F)
dpi = 300
png(file=paste0("analysis/doublet_removal/doublet_removal_plots/",sample,"_dimplot.png"), 
    width = dpi*5, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

# make a composite identity ----
doublet_ident = paste0(so@meta.data$SCT_snn_res.0.8,"-",so@meta.data$DF.classification)
names(doublet_ident) <- colnames(x = so)
so <- AddMetaData(
  object = so,
  metadata = doublet_ident,
  col.name = 'doublet_ident'
)

# sort this so it levels are in order
Idents(so) = doublet_ident
levels(so) = sort(levels(so))

# make plot of two sex genes 
sex_genes=c("Xist","Ddx3y")
p = DotPlot(so, features=sex_genes)
dpi = 300
png(file=paste0("analysis/doublet_removal/doublet_removal_plots/",sample,"_dotplot.png"), 
    width = dpi*5, height = dpi*10, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

# select singlet only and save
so <- subset(so, subset = DF.classification == "Singlet")

saveRDS(so, 
        paste0(outpath, sample,"_doublet_removed.rds"))

