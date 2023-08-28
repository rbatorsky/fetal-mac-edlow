## This script does GO BP enrichment for top markers in brain and placenta
LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

## Libraries
suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(fgsea)
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(enrichplot)})

outdir="/cluster/tufts/slonimlab/rbator01/mouse_scrna_edlow_2020/data/analysis/markers/obsctr/"


## Read in features in order to use ensembl IDs
convert_ens =read.csv("/cluster/tufts/slonimlab/rbator01/mouse_scrna_edlow_2020/data/all_cellranger/AEb2-10_B41_OBS_BR/outs/filtered_feature_bc_matrix/features.tsv.gz",
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

## first run up until line 96 in cell identification

padj_cut = 0.01
#pct_cut=0.1
lfc_cut=0.4
pct1_cut=0.2
nshow=3
namestring=paste0("padj_", padj_cut, "_filter_lfc_",lfc_cut,"_pct1_",pct1_cut)

brain_m_filter = br_markers
brain_m_filter = brain_m_filter %>% dplyr::filter(., p_val_adj < padj_cut )
brain_m_filter = brain_m_filter %>% dplyr::filter(avg_log2FC > lfc_cut)
brain_m_filter = brain_m_filter %>% dplyr::filter(pct.1 > pct1_cut)

## output this to file
brain_m_filter$X = NULL
brain_m_filter = brain_m_filter %>% dplyr::select(c("gene",everything()))
unique(brain_m_filter$cluster)
write.csv(brain_m_filter, gsub(".csv",paste0(namestring,".csv"),br_marker_file), row.names=F)

## give type column
brain_m_filter$type = 'Brain'

## find top markers
br_top_genes = c()
for (cl in levels(seurat_integrated_br_rename@meta.data$cluster_name)){
  c1 = brain_m_filter %>% dplyr::filter(cluster == cl)
  #c1 = c1 %>% dplyr::arrange(p_val_adj) %>% head(nshow)
  c1 = c1 %>% dplyr::arrange(-avg_log2FC) %>% head(nshow)
  genes=c1$gene
  print(paste0(cl,":",paste(genes,collapse = ", ")))
  br_top_genes = c(br_top_genes, c1$gene)
}

br_top_genes = unique(br_top_genes)

br_dot=DotPlot(object = seurat_integrated_br_rename, 
               features = br_top_genes) + 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,hjust=1, size=20))& theme(axis.title.x = element_blank())& theme(axis.title.y = element_blank())

show(br_dot)

ggsave(br_dot, filename = paste0(outdir,"br_obsctr_dotplot_top3.pdf"),
       width = 13, height = 5, units = "in")

pl_dot_br_genes=DotPlot(object = seurat_integrated_pl_rename, 
               features = br_top_genes,
               cols=c("grey","red")) + 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,hjust=1, size=20))& theme(axis.title.x = element_blank())& theme(axis.title.y = element_blank())

show(pl_dot_br_genes)

## PLACENTA

## df input
placenta_m_filter = pl_markers
placenta_m_filter = placenta_m_filter %>% dplyr::filter(., p_val_adj < padj_cut ) 
placenta_m_filter = placenta_m_filter %>% dplyr::filter(avg_log2FC > lfc_cut)
placenta_m_filter = placenta_m_filter %>% dplyr::filter(pct.1 > pct1_cut)

## output this to file
placenta_m_filter$X = NULL
placenta_m_filter = placenta_m_filter %>% dplyr::select(c("gene",everything()))
write.csv(placenta_m_filter, gsub(".csv",paste0(namestring,".csv"),pl_marker_file), row.names=F)

## give type column
placenta_m_filter$type = 'Placenta'

## top genes
pl_top_genes = c()

seurat_integrated_pl_rename@meta.data$cluster_name = relevel(seurat_integrated_pl_rename@meta.data$cluster_name, ref = "HBC_Pf4")

for (cl in levels(seurat_integrated_pl_rename@meta.data$cluster_name)){
  c1 = placenta_m_filter %>% dplyr::filter(cluster == cl)
  #c1 = c1 %>% dplyr::arrange(p_val_adj) %>% head(nshow)
  c1 = c1 %>% dplyr::arrange(-avg_log2FC) %>% head(nshow)
  
  genes=c1$gene
  print(paste0(cl, ":", paste(genes,collapse = ", ")))
  pl_top_genes = c(pl_top_genes, c1$gene)
}

pl_top_genes = unique(pl_top_genes)


pl_dot=DotPlot(object = seurat_integrated_pl_rename, 
               features = pl_top_genes, 
               cols=c("grey","red"))+ 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,hjust=1, size=20)) & theme(axis.title.x = element_blank())& theme(axis.title.y = element_blank())

show(pl_dot)

ggsave(pl_dot, filename = paste0(outdir,"pl_obsctr_dotplot_top3.pdf"),
       device = cairo_pdf,
       width = 14.5, height = 5.5, units = "in")

br_dot_pl_genes=DotPlot(object = seurat_integrated_br_rename, 
               features = pl_top_genes)+ 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,hjust=1, size=20)) & theme(axis.title.x = element_blank())& theme(axis.title.y = element_blank())

show(br_dot_pl_genes)

# dpi=300
# tiff(file=paste0(plot_markers,save_string_br,"dotplot_top",nshow,"_PLACENTA_padj_",padj_cut,"_",res_label,".tiff"), width = dpi*12, height = dpi*7, units = "px",res = dpi,type='cairo')
# DotPlot(object = seurat_integrated_br_rename, features = pl_top_genes)+ 
#   FontSize(x.title = 6, y.title = 15)  + 
#   theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 90,size=20))
# dev.off()


### Clusterprofiler

## put them together and see which terms are common
# save_string = "compare_br_cl"
# 
# cp_df = rbind(brain_m_filter, placenta_m_filter)
# cp_df_ens= cp_df %>% dplyr::left_join(convert_ens,by="gene")
# cp_df_ens %>% dplyr::filter(is.na(ens))
# 
# ck <- compareCluster(geneCluster = ens~cluster,
#                     data = cp_df_ens,
#                      OrgDb = org.Mm.eg.db,
#                      keyType= 'ENSEMBL',
#                      fun = "enrichGO",
#                      ont="BP",
#                     universe=all_genes,
#                     readable=TRUE)
# 
# saveRDS(ck, file = paste0(outdir, "ck_br_pl_markers_compare_br_cl_",namestring,"_ens_universe.rds"))
# ck_df = data.frame(ck)
# ck_df$ratio = DOSE::parse_ratio(ck@compareClusterResult$GeneRatio)
# write.csv(ck_df, file = paste0(outdir, "ck_br_pl_markers_compare_br_cl_",namestring,"_ens_universe.csv"))


## load ck_lfc

#ck_df = data.frame(ck)
#ck_df$ratio = DOSE::parse_ratio(ck@compareClusterResult$GeneRatio)
#write.csv(ck_df,paste0(outdir, "ck_br_pl_markers_compare_br_cl_",namestring,".csv") )

# cp_df_ens_lfc = cp_df_ens %>% dplyr::filter(abs(avg_log2FC)> 0.58)
# table(cp_df_ens_lfc$cluster)
# ck_lfc <- compareCluster(geneCluster = ens~cluster,
#                      data = cp_df_ens_lfc,
#                      universe = all_genes,
#                      OrgDb = org.Mm.eg.db,
#                      keyType= 'ENSEMBL',
#                      fun = "enrichGO",
#                      ont="BP",
#                      readable = TRUE)
# save(ck_lfc, file = paste0(outdir, "ck_br_pl_markers_lfc_0.58.Rdata"))
# ck_lfc_df = data.frame(ck_lfc)
# ck_lfc_df$ratio = DOSE::parse_ratio(ck_lfc@compareClusterResult$GeneRatio)
# write.csv(ck_lfc_df, file = paste0(outdir, "ck_br_pl_markers__lfc_0.58.csv"))

# 
# outdir="/cluster/tufts/slonimlab/rbator01/mouse_scrna_edlow_2020/data/analysis/markers/obsctr_marker_analysis/"
# rm(ck_lfc)
# load(file = paste0(outdir, "ck_br_pl_markers_lfc_0.58.Rdata"))
# 
# nshow=3
# p = clusterProfiler::dotplot(ck_lfc, showCategory =nshow, font.size = 12) + 
#   theme(axis.text.x=element_text(angle=45, hjust=1))+ 
#   scale_y_discrete(labels=function(x) str_wrap(x, width=1000))
# print(p)


######
## Load clusterprofiler and make plots
####
outdir="/cluster/tufts/slonimlab/rbator01/mouse_scrna_edlow_2020/data/analysis/markers/obsctr/marker_analysis/"
padj_cut = 0.01
#pct_cut=0.1
lfc_cut=0.4
pct1_cut=0.2
nshow=3
namestring=paste0("padj_", padj_cut, "_filter_lfc_",lfc_cut,"_pct1_",pct1_cut)
namestring
ck = readRDS(file = paste0(outdir, "ck_br_pl_markers_compare_br_cl_",namestring,"_ens_universe.rds"))
#write.csv(ck, file = paste0(outdir, "ck_br_pl_markers_compare_br_cl_",namestring,"_ens_universe.csv"))
#unique(ck@compareClusterResult$cluster)
### New terms to select

immune=c(
  "positive regulation of response to external stimulus",
  "leukocyte mediated immunity",
  "leukocyte migration",
  "positive regulation of cytokine production",
  "regulation of immune effector process",
  "response to lipopolysaccharide",
  "response to virus",
  "response to molecule of bacterial origin",
  "regulation of inflammatory response",
  "antigen processing and presentation"
)

cellcell=c("positive regulation of cell adhesion",
  "regulation of cell-cell adhesion",
  "cellular response to hormone stimulus",
  "cell chemotaxis"
)

meta=c("ATP metabolic process",
       "regulation of vasculature development",
       "generation of precursor metabolites and energy",
       "carbohydrate catabolic process",
       "purine ribonucleotide metabolic process")

development=c("gliogenesis",
              "regulation of hemopoiesis",
              "regulation of angiogenesis")


stress=c("response to wounding",
         "ERK1 and ERK2 cascade",
         "neuron death",
        "neuron apoptotic process")


lipid=c("lipid localization",
        "lipid transport",
        "cellular response to lipid")

protein=c("positive regulation of protein kinase activity",
        "protein localization to extracellular region",
        "protein secretion"
)

select_cat=c(immune, cellcell, meta, development, stress, lipid, protein)


## Select clusters only

pl_hbc=c("HBC_Pf4",
         "HBC_Cd72")

br_select_clusters=c("Mg_Ccl5",
                     "Mg_Hspb1",
                     "Mg_Sparc",
                     "Mg_Spp1",
                     "Mg_YSI_Pf4")

pamm_select_clusters=c("PAMM_Spp1","PAMM_S100a9")


# plot for fig 2 ------
 
ck_select = data.frame(ck) %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::filter(Description %in% select_cat)  %>% 
  dplyr::filter(cluster %in% c(br_select_clusters,pl_hbc)) %>%
  mutate(cluster = factor(cluster, levels=c(br_select_clusters,pl_hbc))) %>% 
  mutate(Description = factor(Description, levels=rev(select_cat))) 

# ck_select_notdf = ck %>% 
#   dplyr::filter(p.adjust < 0.05) %>% 
#   dplyr::filter(Description %in% select_cat)  %>% 
#   dplyr::filter(cluster %in% c(br_select_clusters,pl_hbc, pamm_select_clusters)) 
# 
# ## try emap
# ck_select_notdf = pairwise_termsim(ck_select_notdf)
# emapplot(ck_select_notdf)
# ##


## The below is to order the factors
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

object = ck_select
x= "cluster"
colorBy="p.adjust"
showCategory=10

by="count"
size="count"

#by="geneRatio"
#size="geneRatio"


split=NULL
includeAll=TRUE
font.size=12
title=""
label_format = 30
group = FALSE
shape = FALSE
color <- NULL

df <- fortify(object, showCategory=showCategory, by=by,
              includeAll=includeAll, split=split) 

label_func <- default_labeller(label_format)

by2 <- switch(size, rowPercentage = "Percentage", 
              count         = "Count", 
              geneRatio     = "GeneRatio")    

df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)

p <- ggplot(df, aes_string(x = 'cluster', y = "Description", size = by2))      

p = p +  geom_point(aes_string(color = colorBy)) +
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) + ggtitle(title) + 
  DOSE::theme_dose(font.size) + 
  scale_size_continuous(range=c(3, 8)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 

show(p)

ggsave(p, filename = paste0(outdir, "marker_select_cat_goenrich.pdf"),
       device = cairo_pdf,
       width = 8, height = 11,
       units = "in")


