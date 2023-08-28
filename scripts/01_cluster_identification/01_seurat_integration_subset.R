# After examining the clusters, remove clusters to keep only mac/mono cell types ----

# pl
input_sample_string_subset = "mandf_pl_no128-4each_mf_obsctr_no115ref_nop6b2_res_0.4_subset"
seurat_integrated_file = "analysis/mandf_pl_no128-4each_mf_obsctr_no115ref_nop6b2_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg.rds"
select_clusters=c(0,1,2,3,4,5,8,10,11,15,17,19,23,24,29)

# br
# input_sample_string_subset = "mandf_br_no128-4each_mf_obsctr_no115ref_w35_res_0.4_subset"
# seurat_integrated_file = "analysis/mandf_br_no128-4each_mf_obsctr_no115ref_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg.rds"
# select_clusters=c(0,1,2,4,5,7,8,9,11,14,15,18,20,22)

save_string_subset=paste(input_sample_string_subset, integration_method, "nfeatures", nfeatures, "npcadim", npcadim, sep="_")

seurat_integrated = readRDS(seurat_integrated_file)
Idents(object = seurat_integrated) <- res

seurat_subset <- subset(seurat_integrated, idents = c(select_clusters))
DefaultAssay(seurat_subset) <- "RNA"
seurat_subset <- DietSeurat(seurat_subset, assays = "RNA")

# reintegrate subset ----
split_seurat = single_cell_transform(seurat_subset)

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = nfeatures)

# RunPCA
split_seurat <- lapply(X = split_seurat,
                       FUN = RunPCA,
                       verbose = FALSE,
                       features = integ_features)


# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)

print("LOG find anchors")

if(use_ref == FALSE){
  
  print("not using reference dataset")
  integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                          normalization.method = "SCT",
                                          anchor.features = integ_features,
                                          reduction = integration_method)
}else if (use_ref == TRUE){
  
  # Read in input samples
  sample_list=read.csv(input_sample_file)
  reference_sample_list = sample_list %>% dplyr::filter(ref == "Y")
  
  print("using reference dataset")
  print(reference_sample_list$name)
  
  if (length(reference_sample_list) == 0){
    print("ERROR - use_ref is TRUE but no reference samples specified")
    quit(status=1)
  }
  
  reference_grep_string = str_c(reference_sample_list$name, sep = "", collapse = "|")
  reference_dataset <- grep(reference_grep_string,names(split_seurat),value=FALSE)
  
  # Find Neighbors
  integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                          normalization.method = "SCT",
                                          anchor.features = integ_features,
                                          reduction = integration_method,
                                          reference = reference_dataset)
}


seurat_integrated <- IntegrateData(anchorset = integ_anchors,
                                   normalization.method = "SCT")
rm(split_seurat)

# PCA, UMAP, cluster subset ----

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated,
                             dims = 1:npcadim,
                             reduction = "pca")

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated,
                                   dims = 1:npcadim)

# Determine the clusters for various resolutions
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = resolutions)



saveRDS(seurat_integrated, paste("analysis/final_rds/subset/", save_string_subset,"_integrated_seurat_noccreg.rds",sep=""))

# find cluster markers subset ----
for (r in resolutions){
  results=NULL
  for (cl in clusters){
    de_cl <- FindMarkers(so, ident.1 = cl, only.pos = TRUE)
    de_cl$cluster = cl
    de_cl$gene = rownames(de_cl)
    de_cl = de_cl %>% dplyr::select(c('gene',everything()))
    
    
    if (is.null(results)){
      results = de_cl
    }else{
      results = rbind(results, de_cl)
    }
  }
  write.csv(markers, paste0("analysis/markers/",save_string_subset,"_",res,"_defaultparam_allmarkers_noann.csv"))
}


