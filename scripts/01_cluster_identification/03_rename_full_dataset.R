# This script creates an integrated Seurat Object starting with either cellranger outputs or individual seurat objects  
# The workflow uses Seurat SCT integration, with the option of specifying reference samples

LIB='/cluster/tufts/bio/tools/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
        library(tidyverse)
        library(Seurat)
})

#  brain data ----
save_string_br="mandf_br_no128-4each_mf_obsctr_no115ref_w35_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg_singler.rds"
seurat_integrated_br = readRDS(paste0("final_rds/full/",save_string_br))

# set resolution
res = "integrated_snn_res.0.4"
Idents(object = seurat_integrated_br) <- res

# add cell type labels 
br_data =tidyr::tribble(
        ~number, ~name,
        "0","Microglia",
        "1","Macrophage",
        "2","Macrophage",
        "3","Neuron",
        "4","Microglia",
        "5","Macrophage",
        "6","Neuron",
        "7","Microglia",
        "8","Macrophage",
        "9","Macrophage",
        "10","Neuron",
        "11","Microglia",
        "12","Neuron",
        "13","Neuron",
        "14","Monocyte",
        "15","Macrophage",
        "16","Endothelial cell",
        "17","Fibroblast",
        "18","Macrophage",
        "19","high-mt cluster",
        "20","Fibroblast",
        "21","Fibroblast",
        "22","Granulocyte"
)

new_br_idents= deframe(br_data[,1:2])
seurat_integrated_br_rename<- RenameIdents(object = seurat_integrated_br,new_br_idents)
seurat_integrated_br_rename$cell_type<- Idents(seurat_integrated_br_rename)

p = DimPlot(seurat_integrated_br_rename, label=T)
show(p)
ggsave(p, filename = "analysis/plots/br_obsctr_allcluster.pdf",
       device = cairo_pdf,
       width = 7, height = 5, units = "in")


# placenta data ----

save_string_pl="mandf_pl_no128-4each_mf_obsctr_no115ref_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg_singler.rds"
seurat_integrated_pl = readRDS(paste0("final_rds/full/",save_string_pl))

res = "integrated_snn_res.0.4"
Idents(object = seurat_integrated_pl) <- res

pl_data =tidyr::tribble(
        ~number, ~name,
        "0","Macrophage",
        "1","Monocyte",
        "2","Macrophage",
        "3","Macrophage",
        "4","Monocyte",
        "5","Monocyte",
        "6","Fibroblast",
        "7","Monocyte",
        "8","Endothelial cell",
        "9","Fibroblast",
        "10","Fibroblast",
        "11","Granulocyte",
        "12","Trophoblast",
        "13","Macrophage",
        "14","Endothelial cell",
        "15","Macrophage",
        "16","Trophoblast",
        "17","Granulocyte",
        "18","Fibroblast",
        "19","Monocyte",
        "20","high-mt-cluster",
        "21","Monocyte",
        "22","Trophoblast",
        "23","Monocyte",
        "24","Macrophage",
        "25","Endothelial cell",
        "26","Granulocyte",
        "27","NK cell",
        "28","Erythrocyte",
        "29","NK cell",
        "30","Monocyte",
        "31","Fibroblast",
        "32","Granulocyte"
)

new_pl_idents= deframe(pl_data[,1:2])
seurat_integrated_pl_rename<- RenameIdents(object = seurat_integrated_pl,new_pl_idents)
seurat_integrated_pl_rename$cell_type<- Idents(seurat_integrated_pl_rename)
saveRDS(seurat_integrated_pl_rename, "mandf_pl_no128-4each_mf_obsctr_no115ref_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg_singler_rename.rds")
