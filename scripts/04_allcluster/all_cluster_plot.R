# LIB='/cluster/tufts/patralab/rbator01/R_libs/4.3.0'
# .libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')

#  brain data ----
seurat_integrated_br = readRDS("analysis/final_rds/full/mandf_br_no128-4each_mf_obsctr_no115ref_w35_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg_singler_rename.rds")
DefaultAssay(seurat_integrated_br) = "RNA"
Idents(object = seurat_integrated_br) <- "cell_type"

# pl data ----
seurat_integrated_pl = readRDS("analysis/final_rds/full/mandf_pl_no128-4each_mf_obsctr_no115ref_cca_nfeatures_3000_npcadim_50_integrated_seurat_noccreg_singler_rename.rds")
Idents(object = seurat_integrated_pl) <- "RNA"
DefaultAssay(seurat_integrated_pl) = "cell_type"


# deg ------
br_deg = read.csv("analysis/deg_seurat/obsctr/mergesome/mfcombined/mandf_br_alldata_named_cluster_deg.csv")
pl_deg = read.csv("analysis/deg_seurat/obsctr/mergesome/mfcombined/mandf_pl_alldata_named_cluster_deg.csv")

# low stringency threshold, like in rest of paper, are all cell types the same?!?
br_deg_sig = br_deg %>%
  dplyr::filter(cluster != "high-mt cluster") %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)  %>% 
  mutate(dir = ifelse(avg_log2FC > 0, "up","dn"))%>%
  mutate(cluster = ifelse( cluster == "Monocyte","Fetal Brain Monocyte",cluster)) %>%
  mutate(cluster = ifelse( cluster == "Granulocyte","Fetal Brain Granulocyte",cluster))

pl_deg_sig = pl_deg %>%
  dplyr::filter(cluster != "high-mt-cluster") %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>% 
  mutate(dir = ifelse(avg_log2FC > 0, "up","dn"))%>%
  mutate(cluster = ifelse(cluster == "Fetal Monocyte","Fetal Placental Monocyte",cluster)) %>%
  mutate(cluster = ifelse( cluster == "Granulocyte","Fetal Placental Granulocyte",cluster))

cl_to_select = c('Microglia','Microglia_YSI',"Fetal Brain Monocyte","Fetal Brain Granulocyte","HBC","PAMM", "Fetal Placental Monocyte","Fetal Placental Granulocyte")

unique(pl_deg_sig$cluster)
merge_deg = rbind(br_deg_sig, pl_deg_sig) %>%
  dplyr::filter(cluster %in% cl_to_select)

t = as.data.frame(table(merge_deg$cluster,merge_deg$dir))
t$cluster = t$Var1
t$direction_in_obese= t$Var2
t$Freq = ifelse(t$direction_in_obese == "dn",-t$Freq,t$Freq)

t = t %>%
  mutate(cluster = factor(cluster, levels = rev(cl_to_select))) %>%
  mutate(type = ifelse(cluster %in% c('Microglia','Microglia_YSI',"Fetal Brain Monocyte","Fetal Brain Granulocyte"),"Brain","Placenta"))

p = ggplot(data=t, aes(x=cluster, y=Freq, fill=type)) +
  geom_bar(stat="identity") +
  geom_bar(stat="identity", fill="black", aes(alpha=direction_in_obese)) +
  scale_alpha_manual(values=c(0,0.2),
                     name="dir",
                     breaks=c("up","dn"),         # can use to set order
                     labels=c("up","dn")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Number of Obese vs. Control DEG") +
  xlab("")+ coord_flip() + theme_classic()+ theme(axis.text=element_text(size=12),
                                                  axis.title=element_text(size=14))

show(p)
ggsave(p, filename = "analysis/plots/allcluster_subsetnames_deg_barplot_0.25lfcfilter.pdf",
       width = 6, height = 4, units = "in")


table(merge_deg$cluster)

ck <- compareCluster(geneCluster = gene ~ cluster,
                     data = merge_deg,
                     OrgDb = org.Mm.eg.db,
                     keyType="SYMBOL",
                     fun = "enrichGO",
                     ont="BP")

saveRDS(ck, file = "analysis/deg_seurat/obsctr/mergesome/mfcombined/ck_mandf_brpl_select_merge_deg_0.25lfcfilter.rds")
ck_df = data.frame(ck)
ck_df$ratio = DOSE::parse_ratio(ck@compareClusterResult$GeneRatio)
write.csv(ck_df, "analysis/deg_seurat/obsctr/mergesome/mfcombined/ck_mandf_brpl_select_merge_deg_0.25lfcfilter.csv")

ck = readRDS(file = "analysis/deg_seurat/obsctr/mergesome/mfcombined/ck_mandf_brpl_select_merge_deg_0.25lfcfilter.rds")
ck@compareClusterResult$Cluster = factor(ck@compareClusterResult$Cluster, levels = cl_to_select)

p = dotplot(ck, showCategory=3) + 
  FontSize(x.text = 8, y.text= 8) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
show(p)
ggsave(p, filename = "analysis/plots/allcluster_subsetnames_deg_goenrich_0.25lfcfilter.pdf",
       width = 6, height = 7, units = "in")





br_deg_sig = br_deg %>%
  dplyr::filter(cluster != "high-mt cluster") %>%
  mutate(exp_10 = ifelse(pct.1 > 0.1, 1,0)) %>%
  mutate(exp_10  = ifelse(exp_10 ==1, 1,
                                 ifelse(pct.2 > 0.1,1,0))) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.3 & exp_10  == 1) %>% 
  mutate(dir = ifelse(avg_log2FC > 0, "up","dn")) %>%
  mutate(cluster = ifelse( cluster == "Monocyte","Fetal Brain Monocyte",cluster)) %>%
  mutate(cluster = ifelse( cluster == "Granulocyte","Fetal Brain Granulocyte",cluster))
  
pl_deg_sig = pl_deg %>%
  dplyr::filter(cluster != "high-mt-cluster") %>%
  mutate(exp_10 = ifelse(pct.1 > 0.1, 1,0)) %>%
  mutate(exp_10  = ifelse(exp_10 ==1, 1,
                          ifelse(pct.2 > 0.1,1,0))) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.3 & exp_10  == 1) %>% 
  mutate(dir = ifelse(avg_log2FC > 0, "up","dn"))%>%
  mutate(cluster = ifelse(cluster == "Fetal Monocyte","Fetal Placental Monocyte",cluster)) %>%
  mutate(cluster = ifelse( cluster == "Granulocyte","Fetal Placental Granulocyte",cluster))

cl_to_select = c('Microglia','Microglia_YSI',"Fetal Brain Monocyte","Fetal Brain Granulocyte","HBC","PAMM", "Fetal Placental Monocyte","Fetal Placental Granulocyte")

unique(pl_deg_sig$cluster)
merge_deg = rbind(br_deg_sig, pl_deg_sig) %>%
  dplyr::filter(cluster %in% cl_to_select)

t = as.data.frame(table(merge_deg$cluster,merge_deg$dir))
t$cluster = t$Var1
t$direction_in_obese= t$Var2
t$Freq = ifelse(t$direction_in_obese == "dn",-t$Freq,t$Freq)

t = t %>%
  mutate(cluster = factor(cluster, levels = rev(cl_to_select))) %>%
  mutate(type = ifelse(cluster %in% c('Microglia','Microglia_YSI',"Fetal Brain Monocyte","Fetal Brain Granulocyte"),"Brain","Placenta"))

p = ggplot(data=t, aes(x=cluster, y=Freq, fill=type)) +
  geom_bar(stat="identity") +
  geom_bar(stat="identity", fill="black", aes(alpha=direction_in_obese)) +
  scale_alpha_manual(values=c(0,0.2),
                     name="dir",
                     breaks=c("up","dn"),         # can use to set order
                     labels=c("up","dn")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Number of Obese vs. Control DEG") +
  xlab("")+ coord_flip() + theme_classic()+ theme(axis.text=element_text(size=12),
                                                  axis.title=element_text(size=14))

show(p)
ggsave(p, filename = "analysis/plots/allcluster_subsetnames_deg_barplot.pdf",
       width = 6, height = 4, units = "in")


table(merge_deg$cluster)

ck <- compareCluster(geneCluster = gene ~ cluster,
                     data = merge_deg,
                     OrgDb = org.Mm.eg.db,
                     keyType="SYMBOL",
                     fun = "enrichGO",
                     ont="BP")

saveRDS(ck, file = "analysis/deg_seurat/obsctr/mergesome/mfcombined/ck_mandf_brpl_select_merge_deg_exp10.rds")
ck_df = data.frame(ck)
ck_df$ratio = DOSE::parse_ratio(ck@compareClusterResult$GeneRatio)
write.csv(ck_df, "analysis/deg_seurat/obsctr/mergesome/mfcombined/ck_mandf_brpl_select_merge_deg_exp10.csv")

ck = readRDS(file = "analysis/deg_seurat/obsctr/mergesome/mfcombined/ck_mandf_brpl_select_merge_deg_exp10.rds")
ck@compareClusterResult$Cluster = factor(ck@compareClusterResult$Cluster, levels = cl_to_select)

p = dotplot(ck, showCategory=3) + 
  FontSize(x.text = 8, y.text= 8) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
show(p)
ggsave(p, filename = "analysis/plots/allcluster_subsetnames_deg_goenrich.pdf",
       width = 6, height = 7, units = "in")




# find differences in mg and hbc ----

ck = readRDS(file = "analysis/deg_seurat/obsctr/mergesome/mfcombined/ck_mandf_brpl_select_merge_deg_0.25lfcfilter.rds")

ck_mg_hbc = ck %>%
  dplyr::filter(cluster %in% c('HBC','Microglia','Microglia_YSI'))

dotplot(ck_mg_hbc, showCategory = 30)+ 
  FontSize(x.text = 8, y.text= 8) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



# find differences in mg and hbc ----
ck_mg = ck %>%
  dplyr::filter(cluster %in% c('Microglia','Microglia_YSI'))
ck_hbc = ck %>%
  dplyr::filter(cluster %in% c('HBC'))

mg_cat = unique(ck_mg@compareClusterResult$Description)
hbc_cat = unique(ck_hbc@compareClusterResult$Description)

ck_mg_hbc_diff = ck %>%
  dplyr::filter(cluster %in% c('HBC','Microglia','Microglia_YSI')) %>%
  dplyr::filter(Description %in% c(setdiff(mg_cat, hbc_cat), setdiff(hbc_cat, mg_cat)))

p = dotplot(ck_mg_hbc_diff, showCategory = 5)+ 
  FontSize(x.text = 8, y.text= 8) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

show(p)

ggsave(p, filename = "analysis/plots/mg_hbc_diff_subsetnames_deg_goenrich_0.25lfcfilter.pdf",
       width = 5, height = 5.5, units = "in")



# markers, each in own dataset -----
br_markers = read.csv("analysis/markers/obsctr/mandf_br_alldata_named_cluster_markers.csv")
pl_markers = read.csv("analysis/markers/obsctr/mandf_pl_alldata_named_cluster_markers.csv")

Idents(seurat_integrated_br_rename) = "cell_type"
DefaultAssay(seurat_integrated_br_rename)
seurat_integrated_br_rename_sub = subset(seurat_integrated_br_rename, ident = "high-mt cluster", invert= T)

br_top_genes = c()
for (cl in levels(seurat_integrated_br_rename_sub@meta.data$cell_type)){
  c1 = br_markers_sig %>% dplyr::filter(cluster == cl)
  c1 = c1 %>% dplyr::arrange(-avg_log2FC) %>% head(3)
  genes=c1$gene
  print(paste0(cl,":",paste(genes,collapse = ", ")))
  br_top_genes = c(br_top_genes, c1$gene)
}

br_top_genes = unique(br_top_genes)

DotPlot(object = seurat_integrated_br_rename_sub, 
        features = br_top_genes) + 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,hjust=1, size=20))


br_markers_sig = br_markers %>%
  dplyr::filter(cluster != "high-mt cluster") %>%
  dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 1 & pct.1 > 0.25) %>%
  mutate(cluster = ifelse(cluster == "Monocyte", "Monocyte_brain",cluster))%>%
  mutate(cluster = ifelse(cluster == "Granulocyte", "Granulocyte_brain",cluster))

pl_markers_sig = pl_markers %>%
  dplyr::filter(cluster != "high-mt-cluster") %>%
  dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 1 & pct.1 > 0.25)%>%
  mutate(cluster = ifelse(cluster == "Fetal Monocyte", "Monocyte_placenta",cluster))%>%
  mutate(cluster = ifelse(cluster == "Granulocyte", "Granulocyte_placenta",cluster))

unique(br_markers_sig$cluster)
unique(pl_markers_sig$cluster)

cl_to_select = c('Microglia','Microglia_YSI',"Monocyte_brain","Monocyte_placenta","HBC","PAMM")
combine_sig = rbind(br_markers_sig, pl_markers_sig) %>%
  dplyr::filter(cluster %in% cl_to_select)
table(combine_sig$cluster)
ck <- compareCluster(geneCluster = gene ~ cluster,
                     data = combine_sig,
                     OrgDb = org.Mm.eg.db,
                     keyType="SYMBOL",
                     fun = "enrichGO",
                     ont="BP")

saveRDS(ck, file = "analysis/markers/obsctr/ck_mandf_plbr_separatecalc_subsetcluster_markers.rds")
ck_df = data.frame(ck)
ck_df$ratio = DOSE::parse_ratio(ck@compareClusterResult$GeneRatio)
write.csv(ck_df, "analysis/markers/obsctr/ck_mandf_plbr_separatecalc_subsetcluster_markers.csv")

# ggsave(br_dot, filename = paste0(outdir,"br_obsctr_dotplot_top3.pdf"),
#        width = 13, height = 5, units = "in")

ck = readRDS(file = "analysis/markers/obsctr/ck_mandf_plbr_separatecalc_subsetcluster_markers.rds")

dotplot(ck, showCategory=3) + 
  FontSize(x.title = 6, y.title = 6) 

# merged markers ----
# merge br and pl select clusters and calculate markers
merged_markers = read.csv("analysis/markers/obsctr/mandf_plbrmerge_select_markers.csv")
