LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(tidyverse)
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')

# DEG - combined cluster sex stratified  ------

outdir="analysis/deg_seurat/obsctr/mergesome/mfseparate/"
c1 = read.csv(paste0(outdir, "de.mfseparate.cc.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_maleandfemale.csv"))
c1$X = NULL

c1$cluster_sex = paste0(c1$cluster,"_",c1$sex)

c1 = c1 %>% dplyr::select(c('p_val_adj','avg_log2FC','cluster_sex','gene'))

full = NULL

for (i in unique(c1$cluster_sex)){
  f = c1 %>% dplyr::filter(cluster_sex == i)
  
  names(f)[1] <- paste0(i,"_", names(f)[1])
  names(f)[2] <- paste0(i,"_", names(f)[2])
  
  f$cluster_sex = NULL
  if (is.null(full)){
    full = f
  }else{
    full = full_join(full,f, by = "gene")
  }
  
  full = full %>% dplyr::select(c("gene",everything()))
}


write.csv(full, paste0(outdir, "de.mfseparate.cc.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_maleandfemale_ipa.csv"),row.names=FALSE)




#  DEG - Male and female combined analysis ------
#This is what was used to generated the MF1 non-sex-stratified lists for IPA that Andrea analyzed
#c1 = read.csv(paste0(outdir, "de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_final.csv"))

outdir="analysis/deg_seurat/obsctr/mergesome/mfcombined/"
c1 = read.csv("analysis/deg_seurat/obsctr/mergesome/de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_logfcthresh0_3apr23.csv")
c1 = read.csv("analysis/deg_seurat/obsctr/mergesome/de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_sex_mast_4apr23.csv")
table(c1$cluster)

c1$cluster_sex = paste0(c1$cluster,"_",c1$sex)

c1 = c1 %>% dplyr::select(c('p_val_adj','avg_log2FC','cluster_sex','gene'))

full = NULL

for (i in unique(c1$cluster_sex)){
  f = c1 %>% dplyr::filter(cluster_sex == i)
  
  names(f)[1] <- paste0(i,"_", names(f)[1])
  names(f)[2] <- paste0(i,"_", names(f)[2])
  
  f$cluster_sex = NULL
  if (is.null(full)){
    full = f
  }else{
    full = full_join(full,f, by = "gene")
  }
  
  full = full %>% dplyr::select(c("gene",everything()))
}


write.csv(full, 
          "analysis/deg_seurat/obsctr/mergesome/de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_logfcthresh0_3apr23_ipa.csv",
          row.names=FALSE)


# MARKERS -----

outdir="analysis/markers/obsctr/"
br_marker=read.csv(paste0(outdir,"mandf_br_no128-4each_mf_obsctr_no115ref_w35_res_0.4_subset_cca_nfeatures_3000_npcadim_50_res_0.4_remove_13_15_17_rename_padj_0.05_pct.1_0.25_lfc_0.25.csv"),)
pl_marker=read.csv(paste0(outdir,"mandf_pl_no128-4each_mf_obsctr_no115ref_w35_res_0.4_subset_cca_nfeatures_3000_npcadim_50_res_0.4_remove_11_21_rename_padj_0.05_pct.1_0.25_lfc_0.25.csv"),)

head(br_marker)
c1 = br_marker %>% dplyr::select(c('p_val_adj','avg_log2FC','cluster','gene'))
c2 = pl_marker %>% dplyr::select(c('p_val_adj','avg_log2FC','cluster','gene'))

full = NULL
for (i in unique(c1$cluster)){
  f = c1 %>% dplyr::filter(cluster == i)
  head(f)
  names(f)[1] <- paste0(i,"_", names(f)[1])
  names(f)[2] <- paste0(i,"_", names(f)[2])
  
  f$cluster = NULL
  if (is.null(full)){
    full = f
  }else{
    full = full_join(full,f, by = "gene")
  }
  
  full = full %>% dplyr::select(c("gene",everything()))
}

for (i in unique(c2$cluster)){
  f = c2 %>% dplyr::filter(cluster == i)
  head(f)
  names(f)[1] <- paste0(i,"_", names(f)[1])
  names(f)[2] <- paste0(i,"_", names(f)[2])
  
  f$cluster = NULL
  full = full_join(full,f, by = "gene")
  full = full %>% dplyr::select(c("gene",everything()))
}

write.csv(full, paste0(outdir, "markers.obsctr.brpl.fix35_padj_0.05_pct.1_0.25_lfc_0.25_ipa.csv"),row.names=FALSE)
