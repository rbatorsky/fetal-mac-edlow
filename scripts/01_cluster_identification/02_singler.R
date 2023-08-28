# use SingleR package to identify clusters

# HPC R  library
LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(tidyverse)
  library(remotes)
  library(celldex)
  library(SummarizedExperiment)
  library(SingleR)
  library(Seurat)
})

# run singleR on dataset with mouseRNAseq database ----

so_file_name="analysis/final_rds/mandf_br_no128-4each_mf_obsctr_no115ref_res_0.4_subset_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg.rds"
so = readRDS(so_file_name)

# Set resolution
res = "integrated_snn_res.0.2"
Idents(object = so) <- res

# subsample for speed
so.subsampled <- so[, sample(colnames(so), size =10000, replace=F)]

# extract counts from RNA assay
counts <- GetAssayData(object = so.subsampled, assay="RNA", slot = "counts")

# Load mouse RNAseq reference
ref <- celldex::MouseRNAseqData()

# do prediction
pred <- SingleR(test = counts, ref = ref, labels = ref$label.main)

# add prediction labels to seurat object
so.subsampled[["SingleR.main.labels"]] <- pred$labels

# display prediction as table
df = data.frame(table(so.subsampled@meta.data[[res]], so.subsampled$SingleR.main.labels))
write.csv(df, paste0("analysis/", outstring, ".csv"))

