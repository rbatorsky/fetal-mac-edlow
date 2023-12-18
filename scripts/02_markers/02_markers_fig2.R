# Script for creating Fig. 2 in manuscript

LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(Seurat)
  library(tidyverse)
  library(pheatmap)
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')
outdir ="analysis/plots/"
# read in the data ----
br_markers = read.csv("analysis/markers/obsctr/mandf_br_named_cluster_markers.csv")
pl_markers = read.csv("analysis/markers/obsctr/mandf_pl_named_cluster_markers.csv")

br_so = readRDS("analysis/final_rds/subset/br_obsctr_final.rds")
pl_so = readRDS("analysis/final_rds/subset/pl_obsctr_final.rds")

# Fig 2 A,B UMAP ----
umap_br=DimPlot(br_so, 
                label = F,
                cols=br_col, 
                raster=T,
                raster.dpi = c(900, 900))

ggsave(umap_br, filename = paste0(outdir, "umap/br_dimplot_rename_noraster.pdf"),
       device = cairo_pdf,
       width = 7, height = 5, units = "in")

umap_pl=DimPlot(pl_so, 
                label = F,
                cols=pl_col, 
                raster=T,
                raster.dpi = c(900, 900))

ggsave(umap_pl, filename = paste0(outdir, "umap/pl_dimplot_rename_noraster.pdf"),
       device = cairo_pdf,
       width = 7, height = 5, units = "in")

# Fig 2 A,B Markers  ----
padj_cut = 0.01
lfc_cut=0.4
pct1_cut=0.2
namestring=paste0("padj_", padj_cut, "_filter_lfc_",lfc_cut,"_pct1_",pct1_cut)

# filter
br_filter = br_markers  %>% 
  dplyr::filter(., p_val_adj < padj_cut ) %>% 
  dplyr::filter(avg_log2FC > lfc_cut) %>% 
  dplyr::filter(pct.1 > pct1_cut) %>% 
  dplyr::select(c("gene",everything())) %>%
  mutate(type = 'Brain')

pl_filter = pl_markers %>% 
  dplyr::filter(., p_val_adj < padj_cut )  %>% 
  dplyr::filter(avg_log2FC > lfc_cut) %>% 
  dplyr::filter(pct.1 > pct1_cut) %>% 
  dplyr::select(c("gene",everything())) %>%
  mutate(type = 'Placenta')

# top markers
br_top_genes = c()
for (cl in levels(br_so@meta.data$cluster_name)){
  c1 = br_filter %>% dplyr::filter(cluster == cl)
  c1 = c1 %>% dplyr::arrange(-avg_log2FC) %>% head(3)
  genes=c1$gene
  print(paste0(cl,":",paste(genes,collapse = ", ")))
  br_top_genes = c(br_top_genes, c1$gene)
}

br_top_genes = unique(br_top_genes)

pl_top_genes = c()
pl_so@meta.data$cluster_name = relevel(pl_so@meta.data$cluster_name, ref = "HBC_Pf4")

for (cl in levels(pl_so@meta.data$cluster_name)){
  c1 = pl_filter %>% dplyr::filter(cluster == cl)
  c1 = c1 %>% dplyr::arrange(-avg_log2FC) %>% head(3)
  genes=c1$gene
  print(paste0(cl, ":", paste(genes,collapse = ", ")))
  pl_top_genes = c(pl_top_genes, c1$gene)
}

pl_top_genes = unique(pl_top_genes)

# dotplot 
br_dot=DotPlot(object = br_so, 
               features = br_top_genes) + 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,hjust=1, size=20))& theme(axis.title.x = element_blank())& theme(axis.title.y = element_blank())

ggsave(br_dot, filename = paste0(outdir,"br_obsctr_dotplot_top3.pdf"),
       width = 13, height = 5, units = "in")

pl_dot=DotPlot(object = pl_so, 
               features = pl_top_genes, 
               cols=c("grey","red"))+ 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,hjust=1, size=20)) & theme(axis.title.x = element_blank())& theme(axis.title.y = element_blank())

ggsave(pl_dot, filename = paste0(outdir,"pl_obsctr_dotplot_top3.pdf"),
       device = cairo_pdf,
       width = 14.5, height = 5.5, units = "in")

# Fig 2 C, Module Scores ----
outdir="analysis/plots/module/"
markers.bao <- read.csv("~/slonimlab/rbator01/reference_data/thomas_2020/jem_20200891_Sdat_combined_convert.csv")

## add the module score
clusters=unique(markers.bao$cluster)

for (cl in clusters){
  
  clusters=unique(markers.bao$cluster)
  cluster_filt = markers.bao %>% dplyr::filter(cluster == cl)
  feature=cluster_filt$gene_name
  
  all_genes=rownames(br_so@assays$RNA@counts)
  
  feature=intersect(feature,all_genes)
  feature=list(c(feature))
  br_so <- AddModuleScore(
    object = br_so,
    features = feature,
    name = cl
  )
}

all_genes=rownames(pl_so@assays$RNA@counts)
clusters=unique(markers.bao$cluster)

for (cl in clusters){
  cluster_filt = markers.bao %>% dplyr::filter(cluster == cl)
  feature=cluster_filt$gene_name
  
  feature=intersect(feature,all_genes)
  feature=list(c(feature))
  pl_so <- AddModuleScore(
    object = pl_so,
    features = feature,
    name = cl
  )
}
features = c("YS1","Mono1","CS101","SAMac1","Kup_cell1" )

pl_so[['module']]= CreateAssayObject(data = t(x = FetchData(object = pl_so, vars = features)))
br_so[['module']]= CreateAssayObject(data = t(x = FetchData(object = br_so, vars = features)))

pl_avg_nos <- AverageExpression(pl_so)
br_avg_nos <- AverageExpression(br_so)

full=t(pl_avg_nos$module)
full_br=t(br_avg_nos$module)

sub_pl=t(pl_avg_nos$module) %>% as.data.frame() %>% dplyr::select('YS1','Mono1')
sub_br=t(br_avg_nos$module) %>% as.data.frame() %>% dplyr::select('YS1','Mono1')

both = rbind(sub_br, sub_pl)
both_subset = both[c('Mg_YSI_Pf4','Mg_Spp1', 'Mg_Ccl5','Mg_Hspb1','Mg_Sparc','HBC_Cd72','HBC_Pf4','PAMM_Spp1','PAMM_Chil3', 'PAMM_MHCII','PAMM_Ccl8'),]

both_long <- data.frame(both_subset) %>%
  rownames_to_column("id") %>%
  mutate(color = YS1) %>%
  pivot_longer(cols='YS1':'Mono1') %>% 
  mutate(id = factor(id, levels=unique(both_long$id))) 

p <- ggplot(data = both_long) +
  geom_tile(aes(x=reorder(id, -color), y=name, fill = value), colour = "white") +
  scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white") + 
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 10,
                                   hjust=1, vjust = 0.5),
        axis.text.y = element_text(size=10),
        axis.ticks.length = unit(0.15, units = c('lines')),
        legend.title = element_text(size = 10),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.8, units = c('lines'))) +
  ylab("Module Score") +
  xlab("Placenta Cluster")


show(p)
ggsave(p, filename = paste0(outdir, "both_module_heatmap_sub_all_ordered_nocc.pdf"),
       device = cairo_pdf,
       width = 6.5, height = 2, 
       units = "in")

# Fig 2 D, Compare Br and Pl datasets ----
# see 02_compare_datasets_fig2.R

# Fig 2E, Canonical Markers ----
## markers to plot

cannonical=c("C1qa", "Fcrls","Cx3cr1", "Csf1r","Cd163")
cannonical_no163=c("C1qa", "Fcrls","Cx3cr1", "Csf1r")
cannonical_rev = rev(cannonical)
cannonical_no163_rev = rev(cannonical_no163)

cannonical_si=c("Aif1","Trem2","Cd68","Cd14","Itgam","Tmem119","P2ry12")


## umap
can_br_u = FeaturePlot(br_so,
                       reduction = "umap",
                       features = cannonical,
                       raster=T,
                       sort.cell = T,
                       label = F,
                       cols=c("lightgrey","blue"),
                       combine = F)

can_br_u  <- lapply(X = can_br_u , FUN = function(x) x + theme(plot.title = element_text(size = 20)))

for(i in 1:length(can_br_u )) {
  can_br_u [[i]] <- can_br_u [[i]] + NoAxes()
}

can_br_u_comb = cowplot::plot_grid(plotlist = can_br_u, ncol=5 )
show(can_br_u_comb)

ggsave(can_br_u_comb, filename = paste0(outdir,"br_cannonical_umap.pdf"),
       device = cairo_pdf,
       width = 12, height = 3, units = "in")

dpi=300
png(file=paste0(outdir,"br_cannonical_umap.png"), 
    width = dpi*15, height = dpi*3, units = "px",res = dpi,type='cairo')
print(can_br_u_comb)
dev.off()


can_pl_u = FeaturePlot(pl_so,
                       reduction = "umap",
                       features = cannonical,
                       label = F, 
                       raster = T,
                       raster.dpi=c(900,900),
                       sort.cell = T,
                       cols=c("lightgrey","red"),
                       combine = F)

can_pl_u  <- lapply(X = can_pl_u , FUN = function(x) x + theme(plot.title = element_text(size = 20)))


for(i in 1:length(can_pl_u )) {
  can_pl_u [[i]] <- can_pl_u [[i]] + NoAxes()
}


can_pl_u_comb = cowplot::plot_grid(plotlist = can_pl_u, ncol=5 )

show(can_pl_u_comb)
ggsave(can_pl_u_comb, filename = paste(outdir,"pl_cannonical_umap.pdf",sep=""),
       device = cairo_pdf,
       width = 15, height = 3, units = "in")

dpi=300
png(file=paste0(outdir,"pl_cannonical_umap.png"), 
    width = dpi*12, height = dpi*3, units = "px",res = dpi,type='cairo')
print(can_pl_u_comb)
dev.off()

## dot

can_br_dot_si = DotPlot(object = br_so, 
                        features = cannonical_si,
                        scale=T) + 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,size=20, hjust=1))& 
  theme(axis.title.x = element_blank())& theme(axis.title.y = element_blank())

show(can_br_dot_si)
ggsave(can_br_dot_si, filename = paste(outdir,"br_cannonical_si_dot.pdf",sep=""),
       device = cairo_pdf,
       width = 8, height = 4.5, units = "in")

can_br_dot = DotPlot(object = br_so, 
                     features = cannonical,
                     scale=T) + 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,size=20, hjust=1))& 
  theme(axis.title.x = element_blank())& theme(axis.title.y = element_blank())

show(can_br_dot)
ggsave(can_br_dot, filename = paste(outdir,"br_cannonical_dot.pdf",sep=""),
       device = cairo_pdf,
       width = 8, height = 4.5, units = "in")



can_pl_dot = DotPlot(object = pl_so, 
                     features = cannonical,
                     cols=c("grey","red"),
                     scale=T) + 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,size=20, hjust=1))& theme(axis.title.x = element_blank())& 
  theme(axis.title.y = element_blank())

show(can_pl_dot)
ggsave(can_pl_dot, filename = paste(outdir,"pl_cannonical_dot.pdf",sep=""),
       device = cairo_pdf,
       width = 8, height = 5.25, units = "in")


can_pl_dot_si = DotPlot(object = pl_so, 
                        features = cannonical_si,
                        cols=c("grey","red"),
                        scale=T) + 
  FontSize(x.title = 6, y.title = 15) + 
  theme(axis.text.y=element_text(size=20), axis.text.x = element_text(angle = 45,size=20, hjust=1))& theme(axis.title.x = element_blank())& 
  theme(axis.title.y = element_blank())

show(can_pl_dot_si)
ggsave(can_pl_dot_si, filename = paste(outdir,"pl_cannonical_si_dot.pdf",sep=""),
       device = cairo_pdf,
       width = 8, height = 5.25, units = "in")




# violin
can_br_v = VlnPlot(br_so, cannonical_rev, 
                   stack = TRUE, 
                   flip = TRUE,
                   sort = F)+
  theme(legend.position = "none") + 
  xlab("") 

ggsave(can_br_v, filename = paste(outdir,"br_cannonical_violin.pdf",sep=""),
       device = cairo_pdf,
       width = 5, height = 3, units = "in")


can_pl_v = VlnPlot(pl_so, cannonical_rev, 
                   stack = TRUE, 
                   flip = TRUE,
                   sort = F)+
  theme(legend.position = "none") + 
  xlab("") + 
  ylab("")

ggsave(can_pl_v, filename = paste(outdir,"pl_cannonical_violin.pdf",sep=""),
       device = cairo_pdf,
       width = 5, height = 3, units = "in")


# Fig 2G, GO BP enrichment -----

# Read in features in order to use ensembl ID
convert_ens =read.csv("data/all_cellranger/AEb2-10_B41_OBS_BR/outs/filtered_feature_bc_matrix/features.tsv.gz",
                      sep="\t",
                      col.names=c("ens","gene","na"),
                      header=F)
convert_ens = convert_ens %>% dplyr::select(-c('na'))
colnames(convert_ens)

# make unique
convert_ens %>% dplyr::filter(grepl('Ccl19',gene))
convert_ens$gene = make.unique(convert_ens$gene)
convert_ens %>% dplyr::filter(grepl('Ccl19',gene))
all_genes=convert_ens$ens

# run clusterprofiler
cp_df = rbind(br_filter, pl_filter)

cp_df_ens= cp_df %>% 
  dplyr::left_join(convert_ens,by="gene")

cp_df_ens %>% dplyr::filter(is.na(ens))

ck <- compareCluster(geneCluster = ens~cluster,
                    data = cp_df_ens,
                     OrgDb = org.Mm.eg.db,
                     keyType= 'ENSEMBL',
                     fun = "enrichGO",
                     ont="BP",
                    universe=all_genes,
                    readable=TRUE)

saveRDS(ck, file = "analysis/markers/marker_analysis/ck_br_pl_markers_compare_br_cl_ens_universe.rds")

# organize into conceptual groups

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


# Select clusters to plot

pl_hbc=c("HBC_Pf4",
         "HBC_Cd72")

br_select_clusters=c("Mg_Ccl5",
                     "Mg_Hspb1",
                     "Mg_Sparc",
                     "Mg_Spp1",
                     "Mg_YSI_Pf4")

pamm_select_clusters=c("PAMM_Spp1","PAMM_S100a9")
 
ck_select = data.frame(ck) %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::filter(Description %in% select_cat)  %>% 
  dplyr::filter(cluster %in% c(br_select_clusters,pl_hbc)) %>%
  mutate(cluster = factor(cluster, levels=c(br_select_clusters,pl_hbc))) %>% 
  mutate(Description = factor(Description, levels=rev(select_cat))) 

#The below is to order the factors
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

ggsave(p, filename = "markers/marker_analysis/marker_select_cat_goenrich.pdf",
       device = cairo_pdf,
       width = 8, height = 11, units = "in")
