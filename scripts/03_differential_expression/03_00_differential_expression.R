LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})


# read in our data
source("~/scripts/load_obsctr_subset_fix35_data.R")

# 4.2.0 has a bug in avg_log2FC
# 4.1.1 and 4.3.0 are ok

## Output file name and directory
outdir="analysis/deg_seurat/obsctr/mergesome/"
out_string="obsctr_fix35"


# MF combined ----
print("Brain")

results_lb2 = NULL
DefaultAssay(seurat_integrated_br_rename)
for (cl in unique(Idents(seurat_integrated_br_rename))){
  so.subset = subset(seurat_integrated_br_rename, idents=cl)
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
for (cl in unique(Idents(seurat_integrated_pl_rename))){
  so.subset = subset(seurat_integrated_pl_rename, idents=cl)
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

write.csv( results_lb2, paste0(outdir,"de.mfcombined.",out_string,".obs_vs_ctr.latent_pair_lb2_sex_mast_4apr23.csv"), row.names = FALSE)


# GROUP CLUSTERS by cell type ----

seurat_integrated_br_rename@meta.data$cluster_type = ifelse(grepl('^Mg_YSI',seurat_integrated_br_rename@meta.data$cluster_name),'MgYSI',
                                                            ifelse(grepl('^Mg_',seurat_integrated_br_rename@meta.data$cluster_name),'Mg',
                                                                   ifelse(grepl('^HBC',seurat_integrated_br_rename@meta.data$cluster_name),'HBC',
                                                                          ifelse(grepl('^PAMM',seurat_integrated_br_rename@meta.data$cluster_name),'PAMM',
                                                                                 ifelse(grepl('^Mono', seurat_integrated_br_rename@meta.data$cluster_name),'Mono','NA')))))

Idents(seurat_integrated_br_rename) = "cluster_type"

seurat_integrated_pl_rename@meta.data$cluster_type = ifelse(grepl('^Mg_YSI',seurat_integrated_pl_rename@meta.data$cluster_name),'MgYSI',
                                                            ifelse(grepl('^Mg_',seurat_integrated_pl_rename@meta.data$cluster_name),'Mg',
                                                                   ifelse(grepl('^HBC',seurat_integrated_pl_rename@meta.data$cluster_name),'HBC',
                                                                          ifelse(grepl('^PAMM',seurat_integrated_pl_rename@meta.data$cluster_name),'PAMM',
                                                                                 ifelse(grepl('^Mono', seurat_integrated_pl_rename@meta.data$cluster_name),'Mono','NA')))))

Idents(seurat_integrated_pl_rename) = "cluster_type"

# Separate M and F combined cluster -----

print("Brain")

results = NULL
results_lb2 = NULL


for (cl in unique(Idents(seurat_integrated_br_rename))){
  so.subset = subset(seurat_integrated_br_rename, idents=cl)
  Idents(so.subset)<-"group"
  print(cl)
  for (sex_i in c('F','M')){
    print(sex_i)
    so.subset.2 = subset(so.subset, sex == sex_i)
    DefaultAssay(so.subset.2) = "RNA"

    de <- FindMarkers(so.subset.2, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("pair"))

    de$cluster = cl
    de$tissue = "Brain"
    de$sex = sex_i
    de$gene = rownames(de)
    de = de %>% dplyr::select(c('gene',everything()))


    if (sex_i == "F"){
      print("inner")
      print(sex_i)
      de_lb2 <- FindMarkers(so.subset.2, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("prep_batch_2","pair"))
      de_lb2$cluster = cl
      de_lb2$tissue = "Brain"
      de_lb2$sex = sex_i
      de_lb2$gene = rownames(de_lb2)
      de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))

      if(is.null(results_lb2)){
        results_lb2 = de_lb2
      }else{
        results_lb2 = rbind(results_lb2, de_lb2)
      }
    }

    if(is.null(results)){
      results = de
    }else{
      results = rbind(results, de)
    }
  }
}


rm(so.subset)
rm(so.subset.2)

## Placenta
for (cl in unique(Idents(seurat_integrated_pl_rename))){
  so.subset = subset(seurat_integrated_pl_rename, idents=cl)
  Idents(so.subset)<-"group"
  print(cl)
  for (sex_i in c('F','M')){
    print(sex_i)
    so.subset.2 = subset(so.subset, sex == sex_i)
    DefaultAssay(so.subset.2) = "RNA"

    de <- FindMarkers(so.subset.2, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("pair"))
    de$cluster = cl
    de$tissue = "Placenta"
    de$sex = sex_i
    de$gene = rownames(de)
    de = de %>% dplyr::select(c('gene',everything()))
    results = rbind(results, de)

    if(sex_i == "F"){
      print("inner")
      print(sex_i)
      de_lb2 <- FindMarkers(so.subset.2, ident.1 = "OBS", ident.2 = "CTR", test.use="MAST", latent.vars = c("prep_batch_2","pair"))
      de_lb2$cluster = cl
      de_lb2$tissue = "Placenta"
      de_lb2$sex = sex_i
      de_lb2$gene = rownames(de_lb2)
      de_lb2 = de_lb2 %>% dplyr::select(c('gene',everything()))
      results_lb2 = rbind(results_lb2, de_lb2)
    }
  }
}

write.csv( results, paste0(outdir,"de.mfseperate.cc.",out_string,".obs_vs_ctr.latent_pair_mast_4apr23.csv"), row.names = FALSE)
write.csv( results_lb2, paste0(outdir,"de.mfseparate.cc.",out_string,".obs_vs_ctr.latent_pair_lb2_mast_4apr23.csv"), row.names = FALSE)