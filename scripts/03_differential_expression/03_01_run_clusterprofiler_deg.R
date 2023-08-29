LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
REPO="http://cran.us.r-project.org"
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggvenn)
  library(pheatmap)
  library(ggpubr)
  library(ggrepel)
  library(openxlsx)
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')

# Read in features in order to use ensembl IDs ----
convert_ens =read.csv("data/all_cellranger/AEb2-10_B41_OBS_BR/outs/filtered_feature_bc_matrix/features.tsv.gz",
                      sep="\t",
                      col.names=c("ens","gene","na"),
                      header=F)
convert_ens = convert_ens %>% dplyr::select(-c('na'))
colnames(convert_ens)
convert_ens %>% dplyr::filter(grepl('Ccl19',gene))
convert_ens$gene = make.unique(convert_ens$gene)


convert_ens %>% dplyr::filter(grepl('Ccl19',gene))
## all genes are the genes that could have been measured
all_genes=convert_ens$ens


# MALE AND FEMALE together, obs vs. control ----

indir="analysis/deg_seurat/obsctr/mergesome/mfcombined/"
outdir="analysis/deg_seurat/obsctr/mergesome/mfcombined/plots/"

# Male and Female together, obs/ctr
merged_seurat = read.csv(paste0(indir, "de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_final.csv"))
merged_seurat$sex = "MF"

merged_seurat_ens= merged_seurat %>% dplyr::left_join(convert_ens,by="gene")
merged_seurat_ens %>% dplyr::filter(is.na(ens))

## SIG by LFC > 0.25 = is the default
# merged_seurat_sig = merged_seurat_ens %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) >0.25)
# 
# merged_seurat_sig$direction = ifelse( merged_seurat_sig$avg_log2FC > 0,"up","dn")
# 
# t1 = as.data.frame.array(table(merged_seurat_sig$cluster)) 
# t1$cluster = rownames(t1)
# write.xlsx(t1, paste0(indir, "de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_final_table.csv"))
# write.csv(merged_seurat_sig , paste0(indir, "de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.25_final.csv"), row.names = F)
# 
# 
# merged_seurat_sig  = read.csv(paste0(indir, "de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.25_final.csv"))
# 
# # Run the clusterprofiler, BP,  no direction
# 
# ck_indiv_nodir<- compareCluster(geneCluster = ens ~ cluster,
#                                 data = merged_seurat_sig,
#                                 OrgDb = org.Mm.eg.db,
#                                 keyType="ENSEMBL",
#                                 fun = "enrichGO",
#                                 ont="BP",
#                                 universe=all_genes,
#                                 readable=TRUE
# )
# 
# saveRDS(ck_indiv_nodir, file = paste0(outdir, "ck_de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.25_nodir_finalfig3.rds"))
# ck_df = data.frame(ck_indiv_nodir)
# ck_df$ratio = DOSE::parse_ratio(ck_indiv_nodir@compareClusterResult$GeneRatio)
# write.csv(ck_df, paste0(outdir, "ck_de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.2_nodir_finalfig3.csv"))

# try a higher fold change threshold -----


merged_seurat_sig = merged_seurat_ens %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) >0.58)

merged_seurat_sig$direction = ifelse( merged_seurat_sig$avg_log2FC > 0,"up","dn")

t1 = as.data.frame.array(table(merged_seurat_sig$cluster))
t1$cluster = rownames(t1)
write.xlsx(t1, paste0(indir, "de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_padj_0.05_abslfc_0.58_final_table.csv"))
write.csv(merged_seurat_sig , paste0(indir, "de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.58_final.csv"), row.names = F)


# Run the clusterprofiler, BP,  no direction

ck_indiv_nodir<- compareCluster(geneCluster = ens ~ cluster,
                                data = merged_seurat_sig,
                                OrgDb = org.Mm.eg.db,
                                keyType="ENSEMBL",
                                fun = "enrichGO",
                                ont="BP",
                                universe=all_genes,
                                readable=TRUE
)

saveRDS(ck_indiv_nodir, file = paste0(indir, "ck_de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.58_nodir.rds"))
ck_df = data.frame(ck_indiv_nodir)
ck_df$ratio = DOSE::parse_ratio(ck_indiv_nodir@compareClusterResult$GeneRatio)
write.csv(ck_df, paste0(indir, "ck_de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.58_nodir.csv"))




#M and F separate  -----
indir="~/analysis/deg_seurat/obsctr/mergesome/mfseparate/"
outdir="~/analysis/deg_seurat/obsctr/mergesome/mfseparate/plots/"

merged_seurat_nolb2 = read.csv(paste0(indir, "de.obsctr_fix35.obs_vs_ctr.latent_pair_mast_logfcthresh0_23Jun22.csv"))
unique(merged_seurat_nolb2$sex)
merged_seurat_m = merged_seurat_nolb2 %>% dplyr::filter(sex == "M")
merged_seurat_m$X = NULL

merged_seurat_lb2_f = read.csv(paste0(indir, "de.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_23Jun22.csv"))
## check that only has "F"
unique(merged_seurat_lb2_f$sex)
merged_seurat_lb2_f$sex = "F"
merged_seurat_lb2_f$X = NULL

# write combined
merged_seurat = rbind(merged_seurat_m, merged_seurat_lb2_f)
#write.csv(merged_seurat, paste0(indir, "de.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_maleandfemale_23Jun22.csv"))


merged_seurat_sig = merged_seurat %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC)>0.25)
write.csv(merged_seurat_sig, paste0(indir, "de.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_maleandfemale_padj_0.05_abslfc_0.25_23Jun22.csv"))

pl_hbc=c("HBC_Pf4",
         "HBC_Cd72",
         "HBC_cellcycle")

br_clusters=c("Mg_Ccl5",
              "Mg_cellcycle",
              "Mg_Hspb1",
              "Mg_Sparc",
              "Mg_Spp1",
              "Mg_YSI_Pf4",
              "Mg_YSI_cellcycle")

pamm_clusters=c("PAMM_Ccl8",
                "PAMM_Spp1",
                "PAMM_Chil3",
                "PAMM_MHCII",
                "PAMM_S100a9")

mono_clusters=c("Mono_FPl",
                "Mono_FBr")

## Label cell type
merged_seurat_sig$cc1 = ifelse(grepl('^Mg_YSI',merged_seurat_sig$cluster),'MgYSI',
                               ifelse(grepl('^Mg_',merged_seurat_sig$cluster),'Mg',
                                      ifelse(grepl('^HBC',merged_seurat_sig$cluster),'HBC',
                                             ifelse(grepl('^PAMM',merged_seurat_sig$cluster),'PAMM',
                                                    ifelse(grepl('^Mono', merged_seurat_sig$cluster),'Mono','NA')))))


# MF separate, combine cluster, no direction ------------------------------

indir="~/analysis/deg_seurat/obsctr/mergesome/mfseparate/"
outdir="~/analysis/deg_seurat/obsctr/mergesome/mfseparate/plots/"

merged_seurat_nolb2 = read.csv(paste0(indir, "de.mfseperate.cc.obsctr_fix35.obs_vs_ctr.latent_pair_mast_logfcthresh0.csv"))
unique(merged_seurat_nolb2$sex)
merged_seurat_m = merged_seurat_nolb2 %>% dplyr::filter(sex == "M")
merged_seurat_m$X = NULL

merged_seurat_lb2_f = read.csv(paste0(indir, "de.mfseparate.cc.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0.csv"))
## check that only has "F"
unique(merged_seurat_lb2_f$sex)
merged_seurat_lb2_f$sex = "F"
merged_seurat_lb2_f$X = NULL

# write combined
merged_seurat = rbind(merged_seurat_m, merged_seurat_lb2_f)
write.csv(merged_seurat, paste0(indir, "de.mfseparate.cc.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_maleandfemale.csv"))


merged_seurat_sig = merged_seurat %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC)>0.25)


merged_seurat_sig$dir = ifelse(merged_seurat_sig$avg_log2FC > 0, "up", "dn")
merged_seurat_sig$sex_dir = paste0(merged_seurat_sig$sex,"_",merged_seurat_sig$dir)
write.csv(merged_seurat_sig, paste0(indir, "de.mfseparate.cc.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_maleandfemale_padj_0.05_lfc_0.25.csv"))

merged_seurat_sig_ens= merged_seurat_sig %>% dplyr::left_join(convert_ens,by="gene")
merged_seurat_sig_ens %>% dplyr::filter(is.na(ens))

table(merged_seurat_sig$cluster, merged_seurat_sig$sex)

ck_nodir <- compareCluster(geneCluster = ens ~ cluster + sex ,
                           data = merged_seurat_sig_ens,
                           OrgDb = org.Mm.eg.db,
                           keyType="ENSEMBL",
                           fun = "enrichGO",
                           ont="BP",
                           universe=all_genes,
                           readable=TRUE
)

saveRDS(ck_nodir , file = paste0(indir, "ck_combinecluster_fix35_latentpair_latentb2_padj_0.05_lfc_0.25_nodir_23Jun22.rds"))
ck_nodir_df = data.frame(ck_nodir)
ck_nodir_df$ratio = DOSE::parse_ratio(ck_nodir@compareClusterResult$GeneRatio)
write.csv(ck_nodir_df, paste0(indir, "ck_combinecluster_fix35_latentpair_latentb2_padj_0.05_lfc_0.25_nodir_23Jun22.csv"))