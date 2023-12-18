# Create an integrated Seurat Object starting with either cellranger outputs or individual seurat objects  
# The workflow uses Seurat SCT integration, with the option of specifying reference samples

# HPC R library
LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(enrichplot)
  library(ggvenn)
  library(ggsankey)
  library(future)
})

library(speckle)
library(sc2marker)

options(future.globals.maxSize = 100000 * 1024^2)

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')


# Parameters ----
doublet_dir="analysis/doublet_removal/"
metadata_file="/data/metadata/mandf_combined_no128.csv"

print(doublet_dir)

# OBS/CTR PL 
input_sample_string = "mandf_pl_no128-4each_mf_obsctr_no115ref_nop6b2"

# OBS/CTR BR 
# input_sample_string = "mandf_br_no128-4each_mf_obsctr_no115ref_w35"

input_sample_file=paste0("data/all_cellranger/run_parameters/",input_sample_string,".csv")

integration_method="cca"
nfeatures=3000
npcadim=50
resolutions=c(0.2,0.4,0.6)

# Output file name and directory
save_string=paste(input_sample_string, integration_method, "nfeatures", nfeatures, "npcadim", npcadim, sep="_")


# Create merged seurat from doublet removed seurat  ----

#Sample input list
sample_list=read.csv(input_sample_file)
sample_list=sample_list$name

# Read in all samples
for (sample_s in sample_list) {
  print(sample_s)
  so = readRDS(paste0(doublet_dir, sample_s,"_doublet_removed.rds"))
  so@meta.data$orig.ident = sample_s
  so = RenameCells(so, new.names=paste0(sample_s,"_",colnames(so)))
  
  # trim metadata columns
  col_to_keep = c("orig.ident","nCount_RNA","nFeature_RNA")
  all_col = colnames(so@meta.data)
  col_to_remove = all_col[!(all_col %in% col_to_keep)]
  for (col in col_to_remove){
    so[[col]] <- NULL
  }
  sample_s = gsub("-","_",sample_s)
  assign(sample_s, so)
}

underscore_list=gsub("-","_",sample_list)
merged_seurat <- merge(x = get(underscore_list[1]), y = mget(underscore_list[-1]))
DefaultAssay(merged_seurat) = "RNA"
merged_seurat <- DietSeurat(merged_seurat, assays = "RNA")

# Add extra metadata
extra_metadata = read_csv(metadata_file)
sample_info = as.data.frame(colnames(merged_seurat))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = merged_seurat@meta.data$orig.ident
sample_info = dplyr::left_join(sample_info,extra_metadata)
rownames(sample_info) = sample_info$ID
merged_seurat = AddMetaData(object = merged_seurat, metadata = sample_info)
merged_seurat[['ID']] = NULL

# mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100


# Filter cells and genes ----

# Cell level filtering
filtered_seurat <- subset(x = merged_seurat,
                          subset= (nCount_RNA >= 500) &
                            (nFeature_RNA >= 250) &
                            (mitoRatio < 0.20))

rm(merged_seurat)

# Gene level filtering
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

rm(filtered_counts)

# SCT without cell cycle scoring ----

split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
rm(filtered_seurat)


for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  #split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  #split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio","S.Score","G2M.Score"))
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
  
}


# Integrate using reference samples ----
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = nfeatures)

## RunPCA
split_seurat <- lapply(X = split_seurat,
                       FUN = RunPCA,
                       verbose = FALSE,
                       features = integ_features)


# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)

print("LOG find anchors")

# Read in input samples
sample_list=read.csv(input_sample_file)
reference_sample_list = sample_list %>% dplyr::filter(ref == "Y")

print("using reference dataset")
print(reference_sample_list$name)

if (length(reference_sample_list) == 0){
  print("ERROR - no reference samples specified")
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


# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors,
                                   normalization.method = "SCT")
rm(split_seurat)

# PCA, UMAP, cluster ----

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


saveRDS(seurat_integrated, paste("analysis/final_rds/full/",save_string,"_integrated_seurat_noccreg.rds",sep=""))

# Find Cluster Markers ----

# Assign identity of clusters
Idents(object = seurat_integrated) <- res

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for DEG purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

clusters=unique(so@meta.data$sub_cluster_names)
colnames(so@meta.data)



# calc markers at various resolutions 
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
  write.csv(markers, paste0("analysis/markers/",save_string,"_",r,"_defaultparam_allmarkers_noann.csv"))
}