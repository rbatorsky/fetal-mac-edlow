# This script creates an integrated Seurat Object starting with either cellranger outputs or individual seurat objects  
# The workflow uses Seurat SCT integration, with the option of specifying reference samples

LIB='/cluster/tufts/patralab/rbator01/R_libs/4.3.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')


# placenta data ----
seurat_integrated_pl = readRDS("analysis/final_rds/full/mandf_pl_no128-4each_mf_obsctr_no115ref_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg_singler_rename.rds")
Idents(object = seurat_integrated_pl) <- "integrated_snn_res.0.4"
DefaultAssay(seurat_integrated_pl) = "RNA"

# old, now use HBC, PAMM subsets below
# pl_data =tidyr::tribble(
#   ~number, ~name,
#   "0","Macrophage",
#   "1","Monocyte",
#   "2","Macrophage",
#   "3","Macrophage",
#   "4","Monocyte",
#   "5","Monocyte",
#   "6","Fibroblast",
#   "7","Monocyte",
#   "8","Endothelial cell",
#   "9","Fibroblast",
#   "10","Fibroblast",
#   "11","Granulocyte",
#   "12","Trophoblast",
#   "13","Macrophage",
#   "14","Endothelial cell",
#   "15","Macrophage",
#   "16","Trophoblast",
#   "17","Granulocyte",
#   "18","Fibroblast",
#   "19","Monocyte",
#   "20","high-mt-cluster",
#   "21","Monocyte",
#   "22","Trophoblast",
#   "23","Monocyte",
#   "24","Macrophage",
#   "25","Endothelial cell",
#   "26","Granulocyte",
#   "27","NK cell",
#   "28","Erythrocyte",
#   "29","NK cell",
#   "30","Monocyte",
#   "31","Fibroblast",
#   "32","Granulocyte"
# )
# 
#  use HBC, PAMM subsets
pl_data =tidyr::tribble(
  ~number, ~name,
  "0","PAMM",
  "1","PAMM",
  "2","HBC",
  "3","HBC",
  "4","PAMM",
  "5","PAMM",
  "6","Fibroblast",
  "7","PAMM",
  "8","Endothelial cell",
  "9","Fibroblast",
  "10","Fibroblast",
  "11","PAMM",
  "12","Trophoblast",
  "13","PAMM",
  "14","Endothelial cell",
  "15","HBC",
  "16","Trophoblast",
  "17","Granulocyte",
  "18","Fibroblast",
  "19","Fetal Monocyte",
  "20","high-mt-cluster",
  "21","PAMM",
  "22","Trophoblast",
  "23","PAMM",
  "24","PAMM",
  "25","Endothelial cell",
  "26","Granulocyte",
  "27","NK cell",
  "28","Erythrocyte",
  "29","NK cell",
  "30","PAMM",
  "31","Fibroblast",
  "32","Granulocyte"
)

new_pl_idents= deframe(pl_data[,1:2])
seurat_integrated_pl_rename<- RenameIdents(object = seurat_integrated_pl,new_pl_idents)
seurat_integrated_pl_rename$cell_type<- Idents(seurat_integrated_pl_rename)

#saveRDS(seurat_integrated_pl_rename, "mandf_pl_no128-4each_mf_obsctr_no115ref_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg_singler_rename.rds")

# p = DimPlot(seurat_integrated_pl_rename, label=T)
# show(p)
# ggsave(p, filename = "analysis/plots/pl_obsctr_allcluster_subsetnames.pdf",
#        device = cairo_pdf,
#        width = 7, height = 5, units = "in")
# 

# calc pl markers ---- 
DefaultAssay(seurat_integrated_pl_rename) = "RNA"
Idents(seurat_integrated_pl_rename) = "cell_type"
markers <- FindAllMarkers(seurat_integrated_pl_rename, only.pos = TRUE)
write.csv(markers, "analysis/markers/obsctr/mandf_pl_alldata_named_cluster_markers.csv")

# calc pl deg -----
# Idents(seurat_integrated_pl_rename) = "cell_type"
# results_lb2 = NULL
# 
# DefaultAssay(seurat_integrated_pl_rename)
# # calc markers ---- 
# DefaultAssay(seurat_integrated_pl_rename) = "RNA"
# DefaultAssay(seurat_integrated_pl_rename)
# 
# for (cl in unique(Idents(seurat_integrated_pl_rename))){
#   so.subset = subset(seurat_integrated_pl_rename, idents=cl)
#   DefaultAssay(so.subset) = "RNA"
#   Idents(so.subset)<-"group"
#   print(cl)
# 
#   de_lb2 <- FindMarkers(so.subset, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("prep_batch_2","pair", "sex"))
#   de_lb2$cluster = cl
#   de_lb2$tissue = "Placenta"
#   de_lb2$sex = 'both'
#   de_lb2$gene = rownames(de_lb2)
#   de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
#   results_lb2 = rbind(results_lb2, de_lb2)
# }
# 
# write.csv( results_lb2, "analysis/deg_seurat/obsctr/mergesome/mfcombined/mandf_pl_alldata_named_cluster_deg.csv", row.names = FALSE)