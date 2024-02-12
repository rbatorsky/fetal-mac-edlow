LIB='/cluster/tufts/patralab/rbator01/R_libs/4.3.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})


# read in our data
setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')

# read in the data ----
br_so = readRDS("analysis/final_rds/subset/br_obsctr_final.rds")
pl_so = readRDS("analysis/final_rds/subset/pl_obsctr_final.rds")

Idents(br_so) = "cluster_name"
Idents(pl_so) = "cluster_name"

# # 4.2.0 has a bug in avg_log2FC
# # 4.1.1 and 4.3.0 are ok
# 
# # MF combined ----
# print("Brain")
# 
# results_lb2 = NULL
# DefaultAssay(br_so)
# for (cl in unique(Idents(br_so))){
#   so.subset = subset(br_so, idents=cl)
#   Idents(so.subset)<-"group"
#   DefaultAssay(so.subset) = "RNA"
#   print(cl)
# 
#   de_lb2 <- FindMarkers(so.subset, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("prep_batch_2","sex","pair"))
#   de_lb2$cluster = cl
#   de_lb2$tissue = "Brain"
#   de_lb2$sex = 'both'
#   de_lb2$gene = rownames(de_lb2)
#   de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
# 
#   if(is.null(results_lb2)){
#     results_lb2 = de_lb2
#   }else{
#     results_lb2 = rbind(results_lb2, de_lb2)
#   }
# 
# }
# 
# 
# rm(so.subset)
# 
# ## Placenta
# for (cl in unique(Idents(pl_so))){
#   so.subset = subset(pl_so, idents=cl)
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
# 
# 
# }
# 
# write.csv( results_lb2, "analysis/deg_seurat/obsctr/mergesome/de.mfcombined.obs_vs_ctr.latent_pair_lb2_sex_mast_4apr23.csv", row.names = FALSE)
# 

# GROUP CLUSTERS by cell type ----

br_so@meta.data$cluster_type = ifelse(grepl('^Mg_YSI',br_so@meta.data$cluster_name),'MgYSI',
                                                            ifelse(grepl('^Mg_',br_so@meta.data$cluster_name),'Mg',
                                                                   ifelse(grepl('^HBC',br_so@meta.data$cluster_name),'HBC',
                                                                          ifelse(grepl('^PAMM',br_so@meta.data$cluster_name),'PAMM',
                                                                                 ifelse(grepl('^Mono', br_so@meta.data$cluster_name),'Mono','NA')))))

Idents(br_so) = "cluster_type"

pl_so@meta.data$cluster_type = ifelse(grepl('^Mg_YSI',pl_so@meta.data$cluster_name),'MgYSI',
                                                            ifelse(grepl('^Mg_',pl_so@meta.data$cluster_name),'Mg',
                                                                   ifelse(grepl('^HBC',pl_so@meta.data$cluster_name),'HBC',
                                                                          ifelse(grepl('^PAMM',pl_so@meta.data$cluster_name),'PAMM',
                                                                                 ifelse(grepl('^Mono', pl_so@meta.data$cluster_name),'Mono','NA')))))

Idents(pl_so) = "cluster_type"
# MF together,  combined cluster -----
print("Brain")

results_lb2 = NULL
DefaultAssay(br_so)
for (cl in unique(Idents(br_so))){
  so.subset = subset(br_so, idents=cl)
  Idents(so.subset)<-"group"
  DefaultAssay(so.subset) = "RNA"
  print(cl)
  
  de_lb2 <- FindMarkers(so.subset, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("prep_batch_2","sex","pair"))
  de_lb2$cluster = cl
  de_lb2$tissue = "Brain"
  de_lb2$sex = 'both'
  de_lb2$gene = rownames(de_lb2)
  de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
  
  if(is.null(results_lb2)){
    results_lb2 = de_lb2
  }else{
    results_lb2 = rbind(results_lb2, de_lb2)
  }
}


rm(so.subset)

## Placenta
for (cl in unique(Idents(pl_so))){
  so.subset = subset(pl_so, idents=cl)
  DefaultAssay(so.subset) = "RNA"
  Idents(so.subset)<-"group"
  print(cl)
  
  de_lb2 <- FindMarkers(so.subset, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("prep_batch_2","pair", "sex"))
  de_lb2$cluster = cl
  de_lb2$tissue = "Placenta"
  de_lb2$sex = 'both'
  de_lb2$gene = rownames(de_lb2)
  de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
  results_lb2 = rbind(results_lb2, de_lb2)
  
}

write.csv( results_lb2, "analysis/deg_seurat/obsctr/mergesome/de.mfcombined.cc.obs_vs_ctr.latent_pair_lb2_sex_mast_31jan24.csv", row.names = FALSE)

# remove latent vars ----
print("Brain")

results_lb2 = NULL
DefaultAssay(br_so)
for (cl in unique(Idents(br_so))){
  so.subset = subset(br_so, idents=cl)
  Idents(so.subset)<-"group"
  DefaultAssay(so.subset) = "RNA"
  print(cl)
  
  de_lb2 <- FindMarkers(so.subset, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("pair"))
  de_lb2$cluster = cl
  de_lb2$tissue = "Brain"
  de_lb2$sex = 'both'
  de_lb2$gene = rownames(de_lb2)
  de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
  
  if(is.null(results_lb2)){
    results_lb2 = de_lb2
  }else{
    results_lb2 = rbind(results_lb2, de_lb2)
  }
}


rm(so.subset)

## Placenta
for (cl in unique(Idents(pl_so))){
  so.subset = subset(pl_so, idents=cl)
  DefaultAssay(so.subset) = "RNA"
  Idents(so.subset)<-"group"
  print(cl)
  
  de_lb2 <- FindMarkers(so.subset, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("pair"))
  de_lb2$cluster = cl
  de_lb2$tissue = "Placenta"
  de_lb2$sex = 'both'
  de_lb2$gene = rownames(de_lb2)
  de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
  results_lb2 = rbind(results_lb2, de_lb2)
  
  
}

write.csv( results_lb2, "analysis/deg_seurat/obsctr/mergesome/de.mfcombined.cc.obs_vs_ctr.latent_pair_mast_31jan24.csv", row.names = FALSE)



# Separate M and F combined cluster -----
# 
# print("Brain")
# 
# results = NULL
# results_lb2 = NULL
# 
# 
# for (cl in unique(Idents(br_so))){
#   so.subset = subset(br_so, idents=cl)
#   Idents(so.subset)<-"group"
#   print(cl)
#   for (sex_i in c('F','M')){
#     print(sex_i)
#     so.subset.2 = subset(so.subset, sex == sex_i)
#     DefaultAssay(so.subset.2) = "RNA"
# 
#     de <- FindMarkers(so.subset.2, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("pair"))
# 
#     de$cluster = cl
#     de$tissue = "Brain"
#     de$sex = sex_i
#     de$gene = rownames(de)
#     de = de %>% dplyr::select(c('gene',everything()))
# 
# 
#     if (sex_i == "F"){
#       print("inner")
#       print(sex_i)
#       de_lb2 <- FindMarkers(so.subset.2, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("prep_batch_2","pair"))
#       de_lb2$cluster = cl
#       de_lb2$tissue = "Brain"
#       de_lb2$sex = sex_i
#       de_lb2$gene = rownames(de_lb2)
#       de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
# 
#       if(is.null(results_lb2)){
#         results_lb2 = de_lb2
#       }else{
#         results_lb2 = rbind(results_lb2, de_lb2)
#       }
#     }
# 
#     if(is.null(results)){
#       results = de
#     }else{
#       results = rbind(results, de)
#     }
#   }
# }
# 
# 
# rm(so.subset)
# rm(so.subset.2)
# 
# ## Placenta
# for (cl in unique(Idents(pl_so))){
#   so.subset = subset(pl_so, idents=cl)
#   Idents(so.subset)<-"group"
#   print(cl)
#   for (sex_i in c('F','M')){
#     print(sex_i)
#     so.subset.2 = subset(so.subset, sex == sex_i)
#     DefaultAssay(so.subset.2) = "RNA"
# 
#     de <- FindMarkers(so.subset.2, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("pair"))
#     de$cluster = cl
#     de$tissue = "Placenta"
#     de$sex = sex_i
#     de$gene = rownames(de)
#     de = de %>% dplyr::select(c('gene',everything()))
#     results = rbind(results, de)
# 
#     if(sex_i == "F"){
#       print("inner")
#       print(sex_i)
#       de_lb2 <- FindMarkers(so.subset.2, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("prep_batch_2","pair"))
#       de_lb2$cluster = cl
#       de_lb2$tissue = "Placenta"
#       de_lb2$sex = sex_i
#       de_lb2$gene = rownames(de_lb2)
#       de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
#       results_lb2 = rbind(results_lb2, de_lb2)
#     }
#   }
# }
# 
# write.csv( results, "analysis/deg_seurat/obsctr/mergesome/de.mfseperate.cc.obs_vs_ctr.latent_pair_mast_4apr23.csv", row.names = FALSE)
# write.csv( results_lb2, "analysis/deg_seurat/obsctr/mergesome/de.mfseparate.cc.obs_vs_ctr.latent_pair_lb2_mast_4apr23.csv", row.names = FALSE)