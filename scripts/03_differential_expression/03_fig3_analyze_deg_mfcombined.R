# Code for Figure 3, Male and Female combined analysis
LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(scales)
  library(pheatmap)
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')

indir="analysis/deg_seurat/obsctr/mergesome/mfcombined/"
outdir="analysis/deg_seurat/obsctr/mergesome/mfcombined/plots/"

## Types of cells
pl_hbc=c("HBC_Pf4",
         "HBC_Cd72")

br_clusters=c("Mg_Ccl5",
              "Mg_Hspb1",
              "Mg_Sparc",
              "Mg_Spp1",
              "Mg_YSI_Pf4")

pamm_clusters=c("PAMM_Ccl8",
                "PAMM_Spp1",
                "PAMM_Chil3",
                "PAMM_MHCII",
                "PAMM_S100a9")

mono_clusters=c("Mono_FPl",
                "Mono_FBr")


non_br=c(pl_hbc, pamm_clusters, "Mono_FPl")

select_mg_hbc = c("Mg_Ccl5","Mg_Hspb1", "Mg_Sparc", "Mg_Spp1", "Mg_YSI_Pf4","HBC_Pf4","HBC_Cd72")

# Read in data to plot ----

# GO results
padj_cut = 0.05

ck=readRDS(file = paste0(indir, "ck_de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.2_nodir_finalfig3.rds"))
ck_filter = ck %>% dplyr::filter(p.adjust < padj_cut)
ck_filter@compareClusterResult$cluster = factor(ck_filter@compareClusterResult$cluster, levels=c(br_clusters, pl_hbc, pamm_clusters, mono_clusters))
ck_hbc_mg = ck_filter %>% dplyr::filter(cluster %in% c(pl_hbc, br_clusters))

# DEG results with zero lfc threshold (for heatmap plotting)
deg = read.csv(paste0(indir, "de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_logfcthresh0_final.csv"))
table(deg$cluster)

# make an example plot
p = clusterProfiler::dotplot(ck_hbc_mg, x =~cluster,showCategory =5, font.size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
show(p)

# 3B Look at the similarities between MgYSI_Pf4 and HBC_Pf4 -----

library(enrichplot)
ck_pf4 = ck_filter %>%
  dplyr::filter(cluster %in% c("HBC_Cd72","Mg_YSI_Pf4"))

ck_pf4_df = data.frame(ck_pf4)

ck_pf4_df = ck_pf4_df %>%
  mutate(GeneRatio =  DOSE::parse_ratio(ck_pf4_df$GeneRatio)) %>%
  dplyr::filter(cluster == "Mg_YSI_Pf4") %>%
  arrange(-GeneRatio) %>%
  head(15)

ck_pf4 = ck_pf4 %>%
  dplyr::filter(Description %in% ck_pf4_df$Description)

ck_pf4 = pairwise_termsim(ck_pf4)
p1 = emapplot(ck_pf4,
              show=30,
              legend_n = 2)+ theme(text=element_text(size=6, family="Arial"))

print(p1)
ggsave(p1, filename = paste0(outdir,"3b_mgysipf4_hbccd72_15_emmapplot_deg_mfcombined_latent_pair_lb2_sex.pdf"),
       device = cairo_pdf,
       width = 5, height = 4.5,
       units = "in")

 
# Fig C. Separate the categories to show functionally similar groups -----

stress=c("neuron death",
         "response to metal ion",
         "response to oxidative stress",
         "stress-activated protein kinase signaling cascade",
         "ERK1 and ERK2 cascade"
)

immune=c("positive regulation of cytokine production",
         "response to molecule of bacterial origin",
         "response to lipopolysaccharide",
         "myeloid leukocyte migration",
         "regulation of tumor necrosis factor production",
         "positive regulation of endocytosis",
         "tumor necrosis factor production",
         "regulation of inflammatory response",
         "humoral immune response")

development=c("myeloid cell differentiation")

metabolic=c("ATP metabolic process",
            "regulation of mRNA metabolic process",
            "regulation of cellular amide metabolic process",
            "generation of precursor metabolites and energy",
            "carbohydrate catabolic process",
            "glycolytic process")

protein=c("protein folding",
          "response to unfolded protein"
)
ribosome=c("cytoplasmic translation",
           "ribonucleoprotein complex assembly",
           "posttranscriptional regulation of gene expression",
           "ribosome assembly",
           "negative regulation of translation",
           "viral genome replication"
)

lipid=c("cellular response to lipid")


select_terms=c(immune,  metabolic, development, stress, lipid, protein, ribosome)

ck_hbc_mg_select = ck_hbc_mg %>% 
  dplyr::filter(Description %in% select_terms)

ck_hbc_mg_select@compareClusterResult$Description = factor(ck_hbc_mg_select@compareClusterResult$Description, levels=rev(select_terms))


## The below is to order the factors
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

object = ck_hbc_mg_select@compareClusterResult
x= "Cluster"
colorBy="p.adjust"
showCategory=15
by="Count"
size="Count"
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

p <- ggplot(df, aes_string(x = 'cluster', y = "Description", size = 'Count'))      

p = p +  geom_point(aes_string(color = colorBy)) + 
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) + ggtitle(title) + 
  DOSE::theme_dose(font.size) + 
  scale_size_continuous(range=c(3, 8)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
show(p)

ggsave(p, filename = paste0(outdir,"select_go_deg_mfcombined_latent_pair_lb2_sex.pdf"),
       device = cairo_pdf,
       width = 8, height = 9, 
       units = "in")



# Fig 3D single category heatmap  -----
select_mg_hbc = c("Mg_Ccl5","Mg_Hspb1","Mg_Sparc","Mg_YSI_Pf4","HBC_Pf4","HBC_Cd72")


# Take the terms we will plot and choose a single scale so they will be comparable
terms=c("regulation of inflammatory response",
                  "response to oxidative stress",
                  "ATP metabolic process")

ck_select_all = ck %>% dplyr::filter(Description %in% terms)
ck_select_all_genes = data.frame(genes = ck_select_all@compareClusterResult$geneID)
ck_select_all_genes$split = gsub("/",",",ck_select_all_genes$genes)
ck_select_all_genes_list = unique(unlist(strsplit(ck_select_all_genes$split,",")))

deg_select_all = deg %>% 
  dplyr::filter(gene %in% ck_select_all_genes_list) %>%
  dplyr::filter(cluster %in% select_mg_hbc) %>%
  arrange(factor(cluster, levels=select_mg_hbc)) %>%
  dplyr::select("gene","avg_log2FC","p_val_adj", "cluster")

# now choose the limits

heatmap_max = max(deg_select_all$avg_log2FC)
heatmap_min = min(deg_select_all$avg_log2FC)


# do individual plots

plot_single <-function(ck, term, outname, heatmap_min, heatmap_max){
  ck_select = ck %>% dplyr::filter(Description == term)
  ck_select_genes = data.frame(genes = ck_select@compareClusterResult$geneID)
  ck_select_genes$split = gsub("/",",",ck_select_genes$genes)
  ck_select_genes_list = unique(unlist(strsplit(ck_select_genes$split,",")))
  
  deg_select = deg %>% 
    dplyr::filter(gene %in% ck_select_genes_list) %>%
    dplyr::filter(cluster %in% select_mg_hbc) %>%
    arrange(factor(cluster, levels=select_mg_hbc)) %>%
    dplyr::select("gene","avg_log2FC","p_val_adj", "cluster")
  
  deg_select$cluster = factor(deg_select$cluster, levels=select_mg_hbc)
  
  lfc_select = deg_select %>% 
    dplyr::select("gene","avg_log2FC", "cluster") %>% 
    spread("cluster","avg_log2FC")
  
  rownames(lfc_select) = lfc_select$gene
  lfc_select$gene = NULL
  lfc_select_mat = as.matrix(lfc_select)
  lfc_select_mat[is.na(lfc_select_mat)] <- 0
  
  p_select = deg_select %>% 
    dplyr::select("gene","p_val_adj", "cluster") %>% 
    spread("cluster","p_val_adj")
  rownames(p_select) = p_select$gene
  p_select$gene = NULL
  p_select_mat = as.matrix(p_select)
  p_select_mat[is.na(p_select_mat)] <- 1
  
  
  
  library(circlize)
  library(ComplexHeatmap)
  ph = Heatmap(lfc_select_mat, 
               col = circlize::colorRamp2(c(heatmap_min, 0, heatmap_max), c("Darkblue", "white", "red")),
               cluster_columns=FALSE, 
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(p_select_mat[i, j] < 0.05 & abs(lfc_select_mat[i,j])>0.2) {
                   grid.text("*", x, y)
                 }
               },
               heatmap_legend_param = list(title = "Scaled log\nfold change")) 
  
  
  
  return(ph)
}

term="ATP metabolic process"
ph = plot_single(ck, term, outname, heatmap_min, heatmap_max)
pdf(paste0(outdir, 'atp_meta.pdf'), height=5.5, width=3)
print(ph)
dev.off()


term="response to oxidative stress"
ph = plot_single(ck, term, outname, heatmap_min, heatmap_max)
pdf(paste0(outdir, 'oxi_stress_heat2.pdf'), height=7.5, width=3)
print(ph)
dev.off()


term="regulation of inflammatory response"
ph = plot_single(ck, term, outname, heatmap_min, heatmap_max)
pdf(paste0(outdir, 'reg_of_inflam_heat2.pdf'), height=7.5, width=3)
print(ph)
dev.off()


# Fig 3E: IPA read in combined ipa data -----
overall_scores = read_excel("analysis/ipa/mf_together/Canonical Pathways/export_cp_mftogether_abslfc_0.25.xlsx", na="NA")
overall_scores$count = sapply(strsplit(as.character(overall_scores$Molecules), ","), length)
#overall_scores[is.na(overall_scores)] = 0

unique(overall_scores$cluster)
select_mg_hbc = c("Mg_Ccl5","Mg_Hspb1", "Mg_Sparc", "Mg_Spp1", "Mg_YSI_Pf4","HBC_Pf4","HBC_Cd72")

select_terms = c(
  "CLEAR Signaling Pathway",  
  "Glycolysis I", 
  "BAG2 Signaling Pathway", 
  "LXR/RXR Activation",
  "NRF2-mediated Oxidative Stress Response", 
  "HIF1Œ± Signaling",
  "FcŒ≥ Receptor-mediated Phagocytosis in Macrophages and Monocytes",
  "EIF2 Signaling",
  "IL-8 Signaling", 
  "Role of PKR in Interferon Induction and Antiviral Response",
  "ERK5 Signaling",  
  "Oxytocin Signaling Pathway", 
  "Ephrin Receptor Signaling",
  "Role of NFAT in Regulation of the Immune Response", 
  "Neuroinflammation Signaling Pathway",
  "Phagosome Formation")

select_terms = gsub("Œ±", "alpha", select_terms)
select_terms = gsub("Œ≥", "gamma", select_terms)

overall_scores_select = data.frame(overall_scores) %>% 
  mutate(Canonical.Pathway = gsub("Œ±", "alpha", `Ingenuity.Canonical.Pathways`)) %>%
  mutate(Canonical.Pathway = gsub("Œ≥", "gamma", Canonical.Pathway)) %>%
  mutate(cluster = gsub("Hsbp1", "Hspb1", cluster)) %>%
  dplyr::filter(Canonical.Pathway %in% select_terms) %>%
  mutate(z.score = as.numeric(z.score)) %>%
  mutate(cluster = factor(cluster, levels=select_mg_hbc)) %>%
  mutate(Canonical.Pathway = factor(Canonical.Pathway, levels = rev(select_terms)))

head(overall_scores_select)
overall_scores_select$`Ingenuity.Canonical.Pathways` = NULL

overall_scores_select$count_na = ifelse(is.na(overall_scores_select$z.score), NA,  overall_scores_select$count)
unique(overall_scores_select$count_na)
p <- ggplot(overall_scores_select, aes_string(x = 'cluster', y = 'Canonical.Pathway', size = 'count_na'))      

p = p +  geom_point(aes_string(color = 'z.score')) + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
  ylab(NULL) +  
  DOSE::theme_dose(12) + 
  scale_size_continuous(range=c(3, 8)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  xlab("") + 
  ylab("")+ 
  theme(rect = element_rect(fill = "transparent"))

show(p)

ggsave(p, filename = paste0(outdir,"ipa_dotplot_can_path_abslfc0.25.pdf"),
       device = cairo_pdf,
       width = 10, height = 5.5, 
       units = "in")


# Fig 3F: Single pathway plots -----
plot_ipa_gene_heatmap <-function(pathway_name, deg, p_thresh, lfc_thresh, select_clusters){
  
  deg_filter = deg %>% 
    dplyr::filter(p_val_adj< p_thresh & abs(avg_log2FC) > lfc_thresh) %>%
    dplyr::filter(cluster %in% select_clusters)
  
  ## IPA pathway genes
  ipa_pathway = read.table(paste0("analysis/ipa/genes_in_cp/",pathway_name,".txt"),
                           header=T, 
                           sep="\t",
                           fill = TRUE )
  
  mouse_convert = read.table("analysis/ipa/mg_ysi_pf4_ipa_export.txt",
                             skip = 2, 
                             header=T, 
                             sep="\t",
                             fill = TRUE )
  
  mouse_convert = mouse_convert %>% dplyr::select(c("ID","Symbol"))
  ipa_pathway = ipa_pathway %>% dplyr::select(c("Symbol"))
  ipa_pathway = ipa_pathway %>% dplyr::inner_join(mouse_convert, by="Symbol")
  
  deg_select = ipa_pathway %>% inner_join(deg, by = c("ID" = "gene"))%>% 
    dplyr::filter(cluster %in% select_clusters) %>%
    mutate(gene = ID) %>%
    dplyr::filter(gene %in% unique(deg_filter$gene))
  
  deg_select$cluster = factor(deg_select$cluster, levels=select_clusters)
  
  lfc_select = deg_select %>% 
    dplyr::select("gene","avg_log2FC", "cluster") %>% 
    spread("cluster","avg_log2FC")
  
  rownames(lfc_select) = lfc_select$gene
  lfc_select$gene = NULL
  lfc_select_mat = as.matrix(lfc_select)
  lfc_select_mat[is.na(lfc_select_mat)] <- 0
  #lfc_select_mat = t(lfc_select_mat)
  
  p_select = deg_select %>% 
    dplyr::select("gene","p_val_adj", "cluster") %>% 
    spread("cluster","p_val_adj")
  
  rownames(p_select) = p_select$gene
  
  p_select$gene = NULL
  p_select_mat = as.matrix(p_select)
  p_select_mat[is.na(p_select_mat)] <- 1
  #p_select_mat=t(p_select_mat)
  
  library(circlize)
  ph = Heatmap(lfc_select_mat, 
               col = circlize::colorRamp2(c(-1, 0, 1), c("Darkblue", "white", "red")),
               column_title =pathway_name,
               cluster_rows=FALSE, 
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(p_select_mat[i, j] < 0.05 & abs(lfc_select_mat[i,j])>0.2) {
                   grid.text("*", x, y)
                 }
               },
               heatmap_legend_param = list(title = "Scaled log\nfold change")) 
  return(ph)
}  

# for male and female together, cluster-level deg

p_thresh=0.05
lfc_thresh=0.25
select_mg_hbc = c("Mg_Ccl5","Mg_Hspb1", "Mg_Sparc", "Mg_Spp1", "Mg_YSI_Pf4","HBC_Pf4","HBC_Cd72")
outdir="analysis/ipa/mf_together/plots/"

deg = read.csv("analysis/deg_seurat/obsctr/mergesome/mfcombined/de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_logfcthresh0_final.csv")

p = plot_ipa_gene_heatmap("CLEAR Signaling Pathway", deg, p_thresh, lfc_thresh,
                          select_clusters = select_mg_hbc)
pdf(paste0(outdir, "clear.pdf"), height=6, width=3.5)
print(p)
dev.off()

p = plot_ipa_gene_heatmap("Fcg Receptor-mediated Phagocytosis in Macrophages and Monocytes", deg, p_thresh, lfc_thresh,
                          select_clusters = select_mg_hbc)
pdf(paste0(outdir, "fcg.pdf"), height=4, width=3.5)
print(p)
dev.off()

p = plot_ipa_gene_heatmap("Neuroinflammation Signaling Pathway", deg, p_thresh, lfc_thresh,
                          select_clusters = select_mg_hbc)
pdf(paste0(outdir, "neuro_inflam.pdf"), height=7, width=3.5)
print(p)
dev.off()



### START HERE IN UPDATING
#####
## Fig 3 A
## Overall correlation of expression changes: use file correlation_of_deg.R
######

indir="~/analysis/deg_seurat/obsctr/mergesome/mfcombined/"
outdir="~/analysis/deg_seurat/obsctr/mergesome/mfcombined/plots/"

library(ggpubr)
library(ggrepel)

merged_seurat = read.csv(paste0(indir, "de.mfcombined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_sex_mast_logfcthresh0_final.csv"))
merged_seurat$X = NULL

pl_hbc=c("HBC_Pf4",
         "HBC_Cd72")

br_clusters=c("Mg_Ccl5",
              "Mg_Hspb1",
              "Mg_Sparc",
              "Mg_Spp1",
              "Mg_YSI_Pf4")


### Need to modify this to plot the black points on top
calc_corr <- function(deg_table,cl1,cl2, padj_thresh, lfc_thresh){
  
  deg_table = merged_seurat
  deg_table$sig = ifelse(deg_table$p_val_adj < padj_thresh & abs(deg_table$avg_log2FC)>lfc_thresh,1,0)
  cl1_table = deg_table %>% dplyr::filter(cluster %in% c(cl1)) %>% dplyr::select(c('gene','avg_log2FC','sig'))
  cl2_table = deg_table %>% dplyr::filter(cluster %in% c(cl2)) %>% dplyr::select(c('gene','avg_log2FC','sig'))
  
  suff1=paste0(".",paste(cl1, collapse=","))
  suff2=paste0(".",paste(cl2, collapse=","))
  name1 = paste0("avg_log2FC", suff1)
  name2=paste0("avg_log2FC", suff2)
  
  
  joined = cl1_table %>% dplyr::full_join(cl2_table, by=c("gene"), suffix=c(suff1,suff2))
  joined$col = ifelse(joined[[paste0('sig.',cl1)]] == 1 & joined[[paste0('sig.',cl2)]] == 1, "black",
                                       ifelse(joined[[paste0('sig.',cl1)]] == 1 & joined[[paste0('sig.',cl2)]] == 0, "blue",
                                              ifelse(joined[[paste0('sig.',cl2)]] == 1 & joined[[paste0('sig.',cl1)]] == 0, "red", "grey")))

  joined = joined %>% arrange(factor(col, levels = c("grey","blue","red","black")))
  corr_value = cor(joined[[name1]], joined[[name2]], method = "spearman",use = "complete.obs")
  
  m <- lm(get(name1) ~ get(name2), joined)
  
  p = ggplot(data=joined,
             aes_string(name2, name1, label = "gene"))+
    geom_point(size=0.1,color = joined$col) + theme_minimal() +
    geom_text_repel(data = . %>% dplyr::filter(col == "black" & 
                                                      abs(get(name1))>0.5 | 
                                                      abs(get(name2))>0.5), max.overlaps = 50,  size=1.75) +
    geom_smooth(method = "lm", se=FALSE, color="red", linetype="dashed",size=0.5, formula = y ~ x)  +
    #geom_text(x = -0.75, y = 1.25, label = paste0("Slope = ", format(unname(coef(m)[2]), digits = 2)),size=1.5) +
    theme(text = element_text(size = 6)) + 
    xlab("") + 
    ylab("")

  
  #show(p)
  ggsave(p, filename = paste0(outdir,cl1,"_",cl2,"_cor.pdf"),
         device = cairo_pdf,
         width = 2.8, height = 2.8, 
         units = "in")
  
  return(corr_value)
}

cor_hbcpf4_mgysi_pf4 = calc_corr(merged_seurat, cl1 = 'Mg_YSI_Pf4',   cl2 = 'HBC_Pf4', padj_thresh=0.05, lfc_thresh = 0)

merged_seurat_sig = merged_seurat %>% 
  dplyr::filter(p_val_adj<0.05) %>%
  dplyr::filter(abs(avg_log2FC)>0.25)

pl_hbc=c("HBC_Pf4",
         "HBC_Cd72")

br_clusters=c("Mg_Ccl5",
              "Mg_Hspb1",
              "Mg_Sparc",
              "Mg_Spp1",
              "Mg_YSI_Pf4")

pamm_clusters=c("PAMM_Ccl8",
                "PAMM_Spp1",
                "PAMM_Chil3",
                "PAMM_MHCII",
                "PAMM_S100a9")

pl_clusters = c(pamm_clusters,pl_hbc)

cor_df = data.frame(matrix(nrow = length(br_clusters), ncol = length(pl_clusters)))
overlap_df = data.frame(matrix(nrow = length(br_clusters), ncol = length(pl_clusters)))
mgcount = data.frame(matrix(nrow = length(br_clusters), ncol = length(pl_clusters)))

rownames(cor_df) = br_clusters
colnames(cor_df) = pl_clusters

rownames(overlap_df) = br_clusters
colnames(overlap_df) = pl_clusters

rownames(mgcount) = br_clusters
colnames(mgcount) = pl_clusters

for (i in 1:length(br_clusters)){
  for (j in 1:length(pl_clusters)){
    brcl = br_clusters[i]
    plcl = pl_clusters[j]
    cor_df[i,j] = calc_corr(merged_seurat, cl1 = brcl,   cl2 = plcl, padj_thresh=0.05, lfc_thresh = 0)
    
    br_gene = merged_seurat_sig %>% 
      dplyr::filter(cluster == brcl)
    
    pl_gene = merged_seurat_sig %>% 
      dplyr::filter(cluster == plcl)
    
    int = intersect(br_gene$gene, pl_gene$gene)
    overlap_df[i,j] = length(int)
    
    mgcount[i,] = length(br_gene$gene)
  }
}

overlap_df_mgfrac = overlap_df/mgcount

#Not the final plot
#pheatmap(overlap_df_mgfrac, cluster_cols=F, cluster_rows=F)
#pheatmap(overlap_df, cluster_cols=F, cluster_rows=F)
#pheatmap(cor_df, cluster_cols=F, cluster_rows=F)
# p = pheatmap(cor_df, cluster_cols=F)
# show(p)
# ggsave(p, filename = paste0(outdir,"mg_hbc_heatmap_cor.pdf"),
#        device = cairo_pdf,
#        width = 3, height = 4, 
#        units = "in")

##########
### heatmap of DEG overlap and correlations
###########
pl_clusters = c(pamm_clusters,pl_hbc)

cor_df = data.frame(matrix(nrow = length(br_clusters), ncol = length(pl_clusters)))
overlap_df = data.frame(matrix(nrow = length(br_clusters), ncol = length(pl_clusters)))

rownames(cor_df) = br_clusters
colnames(cor_df) = pl_clusters

rownames(overlap_df) = br_clusters
colnames(overlap_df) = pl_clusters

for (i in 1:length(br_clusters)){
  for (j in 1:length(pl_clusters)){
    brcl = br_clusters[i]
    plcl = pl_clusters[j]
    cor_df[i,j] = calc_corr(merged_seurat, cl1 = brcl,   cl2 = plcl, padj_thresh=0.05, lfc_thresh = 0)
    
    br_gene = merged_seurat_sig %>% 
      dplyr::filter(cluster == brcl)
    
    pl_gene = merged_seurat_sig %>% 
      dplyr::filter(cluster == plcl)
    
    int = intersect(br_gene$gene, pl_gene$gene)
    overlap_df[i,j] = length(int)
  }
}

#One line per gene 

cl1=br_clusters
cl2=pl_hbc
deg_table=merged_seurat
cl1_table = deg_table %>% dplyr::filter(cluster %in% c(cl1)) %>% dplyr::select(c('gene','avg_log2FC'))
cl2_table = deg_table %>% dplyr::filter(cluster %in% c(cl2)) %>% dplyr::select(c('gene','avg_log2FC'))

cl1_dedup = cl1_table %>%
  group_by(gene) %>%
  dplyr::filter(abs(avg_log2FC) == max(abs(avg_log2FC)))

cl2_dedup = cl2_table %>%
  group_by(gene) %>%
  dplyr::filter(abs(avg_log2FC) == max(abs(avg_log2FC)))


suff1=paste0(".","mg")
suff2=paste0(".","hbc")
name1 = paste0("avg_log2FC", suff1)
name2=paste0("avg_log2FC", suff2)

joined = cl1_dedup %>% dplyr::full_join(cl2_dedup, by=c("gene"), suffix=c(suff1,suff2))
view(joined)
p = ggplot(data=joined, 
           aes_string(name1, name2, label = "gene"))+ 
  geom_point(size=1) +
  stat_cor(method = "spearman", label.x=-1, label.y=1)+ 
  ylab(paste0(cl2,"\nLog2FC")) + 
  xlab(paste0(cl1,"\nLog2FC")) +
  geom_text_repel(max.overlaps=50) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)

print(p)
ggsave(p, filename = paste0(outdir,"mg_hbc_dedup_cor.pdf"),
       device = cairo_pdf,
       width = 5, height = 5, 
       units = "in")






## Don't use this, use the more updated version below
# #########
# ## EGO for a few clusters 
# ## This is the starting point for the single category heatplot
# #########
# run_ego <- function(deg_df, c, all_genes, outdir) {
#   select = deg_df %>% dplyr::filter(cluster == c)
#   genes = select$ens
#   lfc =  select$avg_log2FC
#   names(lfc) <- genes
#   
#   ego <- enrichGO(genes,
#                   OrgDb = org.Mm.eg.db,
#                   ont = "BP",
#                   universe=all_genes,
#                   keyType= 'ENSEMBL',
#                   readable=TRUE)
#   saveRDS(ego, file = paste0(outdir, paste0(c, "_ego.rds")))
#   
#   return(list(ego, lfc))
# }
# ###
# convert_ens =read.csv("/cluster/tufts/slonimlab/rbator01/mouse_scrna_edlow_2020/data/all_cellranger/AEb2-10_B41_OBS_BR/outs/filtered_feature_bc_matrix/features.tsv.gz",
#                       sep="\t",
#                       col.names=c("ens","gene","na"),
#                       header=F)
# convert_ens = convert_ens %>% dplyr::select(-c('na'))
# colnames(convert_ens)
# convert_ens %>% dplyr::filter(grepl('Ccl19',gene))
# convert_ens$gene = make.unique(convert_ens$gene)
# convert_ens %>% dplyr::filter(grepl('Ccl19',gene))
# ## all genes are the genes that could have been measured
# all_genes=convert_ens$ens
# merged_seurat_ens = read.csv(paste0(indir, "de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.2_final.csv"))
# 
# ## SIG by LFC
# merged_seurat_sig = merged_seurat_ens %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) >0.2) 
# mg_ccl5_list = run_ego(merged_seurat_sig, "Mg_Ccl5", all_genes, indir)
# mg_sparc_list = run_ego(merged_seurat_sig, "Mg_Sparc", all_genes, indir)
# mg_hsp_list = run_ego(merged_seurat_sig, "Mg_Hspb1", all_genes, indir)
# mgysi_pf4_list = run_ego(merged_seurat_sig, "Mg_YSI_Pf4", all_genes, indir)
# mg_spp1_list = run_ego(merged_seurat_sig, "Mg_Spp1", all_genes, indir)
# hbc_cd72_list = run_ego(merged_seurat_sig, "HBC_Cd72", all_genes, indir)
# hbc_pf4_list = run_ego(merged_seurat_sig, "HBC_Pf4", all_genes, indir)
# 
# 
# #####
# ## Generate pheat for important clusters
# #####
# 
# select_heat_cat=c("positive regulation of cytokine production",
#                   "response to lipopolysaccharide",
#                   "regulation of inflammatory response",
#                   "response to oxidative stress",
#                   "ATP metabolic process")
# 
# 
# mg_sparc_select=mg_sparc_list[[1]] %>% 
#   dplyr::filter(Description %in% select_heat_cat)
# mg_sparc_pheat=heatplot(mg_sparc_select,
#                          showCategory = 5,
#                          foldChange = mg_sparc_list[[2]]) + 
#   ggplot2::scale_fill_gradient2(low = "blue",
#                                 mid = "white",
#                                 high = "red",
#                                 midpoint = 0)+ 
#   coord_equal() + theme(axis.text=element_text(size=12))+ 
#   ylab("Mg_Sparc")
# ggsave(mg_sparc_pheat, filename = paste0(outdir,"mg_sparc_mfcombined_heatplot.pdf"),
#        device = cairo_pdf,
#        width = 15, height = 2, 
#        units = "in")
# 
# mg_spp1_select=mg_spp1_list[[1]] %>%dplyr::filter(Description %in% select_heat_cat)
# mg_spp1_pheat=heatplot(mg_spp1_select,
#                         showCategory = 5,
#                         foldChange = mg_spp1_list[[2]]) + 
#   ggplot2::scale_fill_gradient2(low = "blue",
#                                 mid = "white",
#                                 high = "red",
#                                 midpoint = 0)+ 
#   coord_equal() + theme(axis.text=element_text(size=12))+ 
#   ylab("Mg_Spp1")
# ggsave(mg_spp1_pheat, filename = paste0(outdir,"mg_spp1_mfcombined_heatplot.pdf"),
#        device = cairo_pdf,
#        width = 15, height = 2, 
#        units = "in")
# 
# 
# mg_ccl5_select=mg_ccl5_list[[1]] %>%dplyr::filter(Description %in% select_heat_cat)
# mg_ccl5_pheat=heatplot(mg_ccl5_select,
#                        showCategory = 5,
#                        foldChange = mg_ccl5_list[[2]])+ 
#   ggplot2::scale_fill_gradient2(low = "blue",
#                                 mid = "white",
#                                 high ="red",
#                                 midpoint = 0)+ 
#   coord_equal() + theme(axis.text=element_text(size=10)) + 
#   ylab("Mg_CCl5")
# show(mg_ccl5_pheat)
# ggsave(mg_ccl5_pheat, filename = paste0(outdir,"mg_ccl5_mfcombined_heatplot.pdf"),
#        device = cairo_pdf,
#        width = 15, height = 2, 
#        units = "in")
# 
# mgysi_pf4_select=mgysi_pf4_list[[1]] %>%dplyr::filter(Description %in% select_heat_cat)
# mgysi_pf4_pheat=heatplot(mgysi_pf4_select,
#                        showCategory = 5,
#                        foldChange = mgysi_pf4_list[[2]]) + 
#   ggplot2::scale_fill_gradient2(low = "blue",
#                                 mid = "white",
#                                 high = "red",
#                                 midpoint = 0)+ 
#   coord_equal() + theme(axis.text=element_text(size=12))+ 
#   ylab("MgYSI_Pf4")
# ggsave(mgysi_pf4_pheat, filename = paste0(outdir,"mgysi_pf4_mfcombined_heatplot.pdf"),
#        device = cairo_pdf,
#        width = 15, height = 2, 
#        units = "in")
# 
# 
# hbc_cd72_select=hbc_cd72_list[[1]] %>%dplyr::filter(Description %in% select_heat_cat)
# hbc_cd72_pheat=heatplot(hbc_cd72_select,
#                        showCategory = 5,
#                        foldChange = hbc_cd72_list[[2]])+ 
#   ggplot2::scale_fill_gradient2(low = "blue",
#                                 mid = "white",
#                                 high = "red",
#                                 midpoint = 0)+ 
#   coord_equal() + theme(axis.text=element_text(size=12))+ 
#   ylab("HBC_Cd72")
# ggsave(hbc_cd72_pheat, filename = paste0(outdir,"hbc_cd72_mfcombined_heatplot.pdf"),
#        device = cairo_pdf,
#        width = 15, height = 2, 
#        units = "in")
# 
# hbc_pf4_select=hbc_pf4_list[[1]] %>%dplyr::filter(Description %in% select_heat_cat)
# hbc_pf4_pheat=heatplot(hbc_pf4_select,
#                         showCategory = 5,
#                         foldChange = hbc_pf4_list[[2]])+ 
#   ggplot2::scale_fill_gradient2(low = "blue",
#                                 mid = "white",
#                                 high = "red",
#                                 midpoint = 0)+ 
#   coord_equal() + theme(axis.text=element_text(size=12))+ 
#   ylab("HBC_Pf4")
# ggsave(hbc_pf4_pheat, filename = paste0(outdir,"hbc_pf4_mfcombined_heatplot.pdf"),
#        device = cairo_pdf,
#        width = 21, height = 2, 
#        units = "in")
# 
# ### Combine heatplot data
# mg_ccl5_heatplot_data = mg_ccl5_pheat$data
# mg_ccl5_heatplot_data$categoryID = paste0("Mg_Ccl5: ",mg_ccl5_heatplot_data$categoryID)
# 
# mg_sparc_heatplot_data = mg_sparc_pheat$data
# mg_sparc_heatplot_data$categoryID = paste0("Mg_Sparc: ",mg_sparc_heatplot_data$categoryID)
# 
# mg_spp1_heatplot_data = mg_spp1_pheat$data
# mg_spp1_heatplot_data$categoryID = paste0("Mg_Spp1: ",mg_spp1_heatplot_data$categoryID)
# 
# 
# mg_hsp1b_heatplot_data = mg_ccl5_pheat$data
# mg_hsp1b_heatplot_data$categoryID = paste0("Mg_Hsp1b: ",mg_hsp1b_heatplot_data$categoryID)
# 
# 
# mgysi_pf4_heatplot_data = mgysi_pf4_pheat$data
# mgysi_pf4_heatplot_data$categoryID = paste0("MgYSI_Pf4: ", mgysi_pf4_heatplot_data$categoryID)
# 
# 
# hbc_cd72_heatplot_data = hbc_cd72_pheat$data
# hbc_cd72_heatplot_data$categoryID = paste0("HBC_Cd72: ", hbc_cd72_heatplot_data$categoryID)
# 
# 
# hbc_pf4_heatplot_data = hbc_pf4_pheat$data
# hbc_pf4_heatplot_data$categoryID = paste0("HBC_Pf4: ", hbc_pf4_heatplot_data$categoryID)
# 
# 
# combined_heatplot = rbind(mg_ccl5_heatplot_data, mg_sparc_heatplot_data, mg_spp1_heatplot_data, mg_hsp1b_heatplot_data, mgysi_pf4_heatplot_data, hbc_pf4_heatplot_data, hbc_cd72_heatplot_data)
# write.csv(combined_heatplot, paste0(outdir, "combined_heatplot_data.csv"))
# 
# 
# plot_single_term <- function(df, term){
#   rng = range(df$foldChange)
#   df_filter = df %>% dplyr::filter(grepl(term, categoryID))
#   df_filter$categoryID = gsub(paste0(": ",term),"",df_filter$categoryID)
#   plot = ggplot(df_filter, aes(categoryID, Gene, fill= foldChange)) + 
#     geom_tile(color = "white",
#               lwd = 1.5,
#               linetype = 1) +
#     coord_fixed() + 
#     ggplot2::scale_fill_gradient2(low = "blue",
#                                   mid = "white",
#                                   high = "red",
#                                   midpoint=0,    #same midpoint for plots (mean of the range)
#                                   limits=c(floor(rng[1]), ceiling(rng[2]))) + #same limits for plots
#     ggtitle("X") + 
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
#     xlab("") + ylab("")+ ggtitle(term)+ coord_flip()+ theme(text = element_text(size = 20))
#   
#   return(plot)
# }  
# 
# select_heat_cat=c("positive regulation of cytokine production",
#                   "response to lipopolysaccharide",
#                   "regulation of inflammatory response",
#                   "response to oxidative stress",
#                   "ATP metabolic process")
# 
# 
# plot = plot_single_term(combined_heatplot, "positive regulation of cytokine production")
# 
# show(plot)
# ggsave(plot, filename = paste0(outdir,"positive_regulation_cytokine_production.mfcombined.heatplot.pdf"),
#        device = cairo_pdf,
#        width = 16, height = 3.4, 
#        units = "in")
# 
# plot = plot_single_term(combined_heatplot, "response to lipopolysaccharide")
# show(plot)
# ggsave(plot, filename = paste0(outdir,"response_to_lps.mfcombined.heatplot.pdf"),
#        device = cairo_pdf,
#        width = 15, height = 3.4, 
#        units = "in")
# 
# plot = plot_single_term(combined_heatplot, "regulation of inflammatory response")
# show(plot)
# ggsave(plot, filename = paste0(outdir,"regulation_inflammatory_response.mfcombined.heatplot.pdf"),
#        device = cairo_pdf,
#        width = 15, height = 3.4, 
#        units = "in")
# 
# plot = plot_single_term(combined_heatplot, "ATP metabolic process")
# show(plot)
# ggsave(plot, filename = paste0(outdir,"atp_metabolic_process.mfcombined.heatplot.pdf"),
#        device = cairo_pdf,
#        width = 10, height = 3.4, 
#        units = "in")
# 
# plot = plot_single_term(combined_heatplot, "response to oxidative stress")
# show(plot)
# ggsave(plot, filename = paste0(outdir,"response_oxidative_stress.mfcombined.heatplot.pdf"),
#        device = cairo_pdf,
#        width = 15, height = 3.4, 
#        units = "in")
# 
# 
### M and F separately
# 
# ck=readRDS(file = paste0(indir, "ck_cluster_fix35_latentpair_latentb2_padj_0.05_lfc_0.25_nodir.rds"))
# padj_cut = 0.05
# ck_filter = ck %>% dplyr::filter(p.adjust < padj_cut)
# 
# ck_filter@compareClusterResult$Cluster = gsub('Mg_IEG','Mg_Hspb1',ck_filter@compareClusterResult$Cluster)
# ck_filter@compareClusterResult$cluster = gsub('Mg_IEG','Mg_Hspb1',ck_filter@compareClusterResult$cluster)
# 
# ## There is a problem here with NA, probably Mg_YSI_CCl24->Pf4 not changed correctly
# unique(ck_filter@compareClusterResult$Cluster)
# 
# ## set the levels for plotting
# ck_filter@compareClusterResult$cluster = factor(ck_filter@compareClusterResult$cluster, levels=c(br_clusters, pl_hbc, pamm_clusters, mono_clusters))
# 
# ## Do Male First
# ck_filter_m = ck_filter %>% dplyr::filter(sex == "M")
# 
# ## Break into cell types
# ck_brain_m= ck_filter_m %>% dplyr::filter(cluster %in% c(br_clusters))
# ck_hbc_m =  ck_filter_m %>% dplyr::filter(cluster %in% c(pl_hbc))
# ck_hbc_mg_m = ck_filter_m %>% dplyr::filter(cluster %in% c(pl_hbc, br_clusters))
# ck_pamm_m = ck_filter_m %>% dplyr::filter(cluster %in% c(pamm_clusters))
# ck_mg_hbc_pamm = ck_filter_m %>% dplyr::filter(cluster %in% c(br_clusters, pl_hbc, pamm_clusters))
# ck_select_mg_hbc = ck_filter_m %>% dplyr::filter(cluster %in% select_mg_hbc)
# 
# # #####
# ## 0. Make a plot of a DEG for illustration
# ####
# 
# seurat_integrated_br_rename_m = subset(seurat_integrated_br_rename, sex == 'M')
# seurat_integrated_pl_rename_m = subset(seurat_integrated_pl_rename, sex == 'M')
# 
# heat_shock=list(c("Hspa1a","Hspa1b", "Hsph1", "Hspa8","Hsph1","Hspe1"))
# 
# seurat_integrated_br_rename_m <- AddModuleScore(
#   object = seurat_integrated_br_rename_m,
#   features = heat_shock,
#   name = 'heat_shock'
# )
# 
# seurat_integrated_pl_rename_m <- AddModuleScore(
#   object = seurat_integrated_pl_rename_m,
#   features = heat_shock,
#   name = 'heat_shock'
# )
# 
# colnames(seurat_integrated_br_rename_m@meta.data)
# FeaturePlot(seurat_integrated_br_rename_m, features = "heat_shock1", split.by = "group")
# FeaturePlot(seurat_integrated_pl_rename_m, features = "heat_shock1", split.by = "group")
# 
# colnames(seurat_integrated_br_rename@meta.data)
# 
# 
# 
# 
# #####
# ## 0. Plot a single category with a heatmap
# ####
# 
# mg_ccl5_ego = readRDS(file = paste0(outdir, paste0("Mg_Ccl5", "_ego.rds")))
# hbc_pf4_ego = readRDS(file = paste0(outdir, paste0("HBC_Pf4", "_ego.rds")))
# hbc_cd72_ego = readRDS(file = paste0(outdir, paste0("HBC_Cd72", "_ego.rds")))
# mg_sparc_ego = readRDS(file = paste0(outdir, paste0("Mg_Sparc", "_ego.rds")))
# mg_ysl_ego = readRDS(file = paste0(outdir, paste0("Mg_YSI_Ccl24", "_ego.rds")))
# 
# 
# 
# 
# 
# rwb <- colorRampPalette(colors = c("blue", "white", "red"))
# 
# mg_ccl5_pheat=heatplot(mg_ccl5_list[[1]],
#                     showCategory = 10,
#                     foldChange = mg_ccl5_list[[2]])+ ggplot2::scale_fill_gradientn(colours = rwb(100))
# 
# show(mg_ccl5_pheat)
# 
# ## Plot the genes in "positive regulation of cytokine production"
# 
# 
# 
# 
# 
# ## Find uniquely overlaping genes Mg/HBC
# merged_seurat_sig_mg_hbc_not_pamm_m = merged_seurat_sig_m %>%
#   dplyr::filter(gene %in% mg_gene) %>%
#   dplyr::filter(gene %in% hbc_gene) %>%
#   dplyr::filter(!(gene %in% pamm_gene))
# 
# 
# 
# ##############
# # 1. all categories with overlap
# ############
# p = clusterProfiler::dotplot(ck_hbc_mg_m, x =~cluster,showCategory =20, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# ##############
# # 1. just mg, hbc in male for selection
# ############
# p = clusterProfiler::dotplot(ck_hbc_mg_m, x =~cluster,showCategory =20, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# 
# ##############
# ## 2. the brain clusters in male
# ###############
# 
# p = clusterProfiler::dotplot(ck_brain_m, x =~cluster,showCategory =8, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# show(p)
# 
# # emmap plot
# termsim = enrichplot::pairwise_termsim(ck_brain_m)
# emapplot(termsim, showCategory = 30)
# 
# ##########
# ## 3. Then look at all three categories
# ##########
# ## set the levels
# ck_mg_hbc_pamm@compareClusterResult$cluster = factor(ck_mg_hbc_pamm@compareClusterResult$cluster, levels=c(br_clusters, pl_hbc, pamm_clusters))
# 
# 
# p = clusterProfiler::dotplot(ck_mg_hbc_pamm, x =~cluster,showCategory =3, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# termsim2 = enrichplot::pairwise_termsim(ck_mg_hbc_pamm)
# 
# library(RColorBrewer)
# require(viridis)
# pal = brewer.pal(n = 8, name = "Set2")
# 
# emapplot(termsim2)
# 
# ##########
# ## 3. Look at the Top 5 categories in all Brain and plot them in all clusters
# ##########
# top_5_each_brain_m <- fortify(ck_brain_m, showCategory = 5)
# 
# 
# ck_filter_m_top5_brain_m = ck_mg_hbc_pamm %>% dplyr::filter(ID %in% top_5_each_brain_m$ID)
# 
# p = clusterProfiler::dotplot(ck_filter_m_top5_brain_m, x =~cluster,showCategory =5, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# 
# ##########
# ## 3. Look at the Top 5 categories in all Brain and plot them in Mg and HBC only
# ##########
# 
# ck_filter_m_top5_brain_m = ck_hbc_mg_m %>% dplyr::filter(ID %in% top_5_each_brain_m$ID)
# 
# p = clusterProfiler::dotplot(ck_filter_m_top5_brain_m, x =~cluster,showCategory =5, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# termsim2 = enrichplot::pairwise_termsim(ck_hbc_mg_m)
# 
# library(RColorBrewer)
# require(viridis)
# pal = brewer.pal(n = 8, name = "Set2")
# 
# emapplot(termsim2) + scale_color_gradient(low = "#132B43", high = "#132B43")
# 
# 
# ##########
# ## 3. Separate the categories to show functionally similar groups
# ##########
# ## These cat are selected in Brain clusters, male, taking top 5 UNIQUE per cluster 
# ###(e.g. if top 5 in a cluster are already selected due to another cluster, we go further down the list until we have 5 unique)
# 
# rna=c("RNA splicing",
#               "cytoplasmic translation")
# protein_folding=c("protein folding",
#                   "response to heat")
# immune=c("leukocyte migration",
#          "tumor necrosis factor production",
#          "leukocyte chemotaxis",
#          "viral process",
#          "symbiotic process",
#          "positive regulation of cytokine production",
#          "cellular response to lipoprotein particle stimulus",
#          "cell chemotaxis")
# 
# metabolic=c("regulation of mRNA metabolic process",
#             "carbohydrate catabolic process",
#             "regulation of cellular amide metabolic process",
#             "positive regulation of proteolysis",
#             "positive regulation of protein kinase activity",
#             "peptidyl-tyrosine phosphorylation",
#             "reactive oxygen species biosynthetic process")
# 
# ribosome=c("ribonucleoprotein complex biogenesis",
#            "ribosome biogenesis",
#            "ribosome assembly"
#            )
# 
# homeostasis=c("myeloid cell homeostasis",
#               "homeostasis of number of cells"
#               )
# 
# development=c("regulation of angiogenesis",
#               "regulation of vasculature development",
#               "myeloid cell differentiation",
#               "positive regulation of neuron differentiation",
#               "developmental cell growth"              )
# 
# cell_death=c("neuron death")
# 
# adhesion=c("positive regulation of cell adhesion",
#            "protein localization to extracellular region")
# 
# localization=c("establishment of protein localization to organelle")
# 
# signal_transduction=c("ERK1 and ERK2 cascade")
# 
# stress_response=c("response to wounding",
#                   "response to oxidative stress")
# 
# 
# ck_select_brain_m_immune = ck_hbc_mg_m %>% dplyr::filter(Description %in% c(immune, rna))
# ck_select_brain_m_rna = ck_hbc_mg_m %>% dplyr::filter(Description %in% c(rna))
# 
# ck_select_brain_m_immune@compareClusterResult$Description = factor(ck_select_brain_m_immune@compareClusterResult$Description, levels=c(immune, rna))
# 
# p = clusterProfiler::dotplot(ck_select_brain_m_immune, x =~cluster,showCategory =5, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# 
# 
# ##########
# ## 3. Look at a specific cluster up and down genes
# ##########
# 
# 
# 
# 
# ##########
# ## 3. Show just Mg_Sparc, Mg_YSI_Ccl5
# ##########
# 
# p = clusterProfiler::dotplot(ck_select_mg_hbc, x =~cluster,showCategory =5, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# 
# 
# ###
# ## 4. Filter down the ID to focus on the categories that are more important to Mg and HBC and less to PAMM
# ###
# ck_m_df = ck_filter_m@compareClusterResult
# id_ck_brain_m = unique(ck_brain_m@compareClusterResult$ID)
# id_ck_hbc_m = unique(ck_hbc_m@compareClusterResult$ID)
# id_ck_pamm_m = unique(ck_pamm_m@compareClusterResult$ID)
# 
# ## rank by rank
# net_df = NULL
# for (i in 1:length(id_ck_brain_m)){
#   id = id_ck_brain_m[i]
#   if (id %in% id_ck_hbc_m ){
#     hbc_sub = ck_hbc_m@compareClusterResult %>% dplyr::filter(ID == id)
#     pamm_sub = ck_pamm_m@compareClusterResult %>% dplyr::filter(ID == id)
#     mg_sub = ck_brain_m@compareClusterResult %>% dplyr::filter(ID == id)
#     
#     
#     unique_mg = length(unique(mg_sub$cluster))
#     unique_hbc = length(unique(hbc_sub$cluster))
#     unique_pamm = length(unique(pamm_sub$cluster))
#     
#     count_mg = sum(mg_sub$Count)
#     count_hbc = sum(hbc_sub$Count)
#     count_pamm = sum(pamm_sub$Count)
#     
#     net_i = data.frame("id" = id, 
#                        "net" = unique_mg + unique_hbc - 2*unique_pamm, 
#                        "count" = count_mg + count_hbc - 2*count_pamm,
#                        "rank" = sum(hbc_sub$rank) + sum(mg_sub$rank) - 4*sum(pamm_sub$rank))
#     
#     if (is.null(net_df)){
#       net_df = net_i
#     }else{
#       net_df = rbind(net_df, net_i)
#     }
#   }
# }
# 
# 
# hist(net_df$count)
# 
# net_df_net = (net_df %>% dplyr::filter(net >0))$id
# net_df_count = (net_df %>% dplyr::filter(count >30))$id
# net_df_rank = (net_df %>% dplyr::filter( rank>0))$id
# 
# 
# ck_m_net = ck_filter_m %>% dplyr::filter(ID %in% net_df_count)
# 
# p = clusterProfiler::dotplot(ck_m_count, x =~cluster,showCategory =5, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# termsim2 = enrichplot::pairwise_termsim(ck_m_net)
# 
# emapplot(termsim2) 
# 
# ##########
# ## 5. See which HBC cluster shows the most Mg cat enrichment terms
# ########
# ## I don't think this is a good idea, too many terms with 1 or 2 of the same genes, counted many times.
# 
# non_br=c(pl_hbc, pamm_clusters, 'Mono_FPl')
# mat = matrix(0, length(non_br), length(br_clusters))
# rownames(mat) = non_br
# colnames(mat) = br_clusters
# 
# for (i in 1:length(non_br)){
#   ck_pl = ck_filter_m %>% dplyr::filter(cluster == non_br[i])
#   go_pl = unique(ck_pl@compareClusterResult$ID)
#   
#   for (j in 1:length(br_clusters)){
#     ck_br = ck_filter_m %>% dplyr::filter(cluster == br_clusters[j])
#     go_br = unique(ck_br@compareClusterResult$ID)
#     overlap = intersect(go_br,go_pl)
#     fraction = length(overlap)/length(go_br)
#     mat[i,j] = fraction
#   }
# }
# library(pheatmap)
# pheatmap(mat)
# view(mat)
# 
# ##########
# ## 5. See which HBC cluster shows the most Mg DEG
# ########
# 
# ## Read in the DEG for the individual clusters
# merged_seurat_nolb2 = read.csv(paste0(outdir, "de.combined.obsctr_fix35.obs_vs_ctr.latent_pair_mast.csv"))
# merged_seurat_m = merged_seurat_nolb2 %>% dplyr::filter(sex == "M")
# merged_seurat_lb2_f = read.csv(paste0(outdir, "de.combined.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast.csv"))
# merged_seurat_lb2_f$sex = "F"
# merged_seurat = rbind(merged_seurat_m,merged_seurat_lb2_f)
# 
# ## LOG FOLD CHANGE ADDED
# merged_seurat_sig = merged_seurat %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC) >0.25) 
# merged_seurat_sig_m = merged_seurat_sig %>% dplyr::filter(sex == "M")
# 
# ## total genes for phyper calc
# #total= length(merged_seurat_nolb2$gene) # All genes in genome
# total= length(merged_seurat_sig_m$gene) # All genes in genome
# 
# 
# ## make an empty matrix to store results
# mat = matrix(0, length(br_clusters), length(non_br))
# rownames(mat) = br_clusters
# colnames(mat) = non_br
# 
# phyper_mat = matrix(0, length(br_clusters), length(non_br))
# rownames(phyper_mat) = br_clusters
# colnames(phyper_mat) = non_br
# 
# for (i in 1:length(non_br)){
#   deg_pl = merged_seurat_sig_m %>% dplyr::filter(cluster == non_br[i])
#   deg_pl = unique(deg_pl$gene)
#   group2 = length(deg_pl)
#   for (j in 1:length(br_clusters)){
#     deg_br = merged_seurat_sig_m %>% dplyr::filter(cluster == br_clusters[j])
#     deg_br = unique(deg_br$gene)
#     group1 = length(deg_br)
#     overlap = length(intersect(deg_br,deg_pl))
#     fraction = overlap/group1
#     phyper_mat[j,i] = phyper(overlap-1, group2, total-group2, group1,lower.tail= FALSE)
#     mat[j,i] = fraction
#   }
# }
# library(pheatmap)
# 
# phyper_mat = formatC(phyper_mat, format = "e", digits = 0)
# 
# png(file=paste0(home, "/male_brain_deg_overlap.png"), width = 300*4, height = 300*3, units = "px",res = 300,type='cairo')
# p = pheatmap(mat,display_numbers =phyper_mat)
# print(p)
# dev.off()
# 
# ##########
# ## 5. Calculate the hypergeometric p-value for the two gene sets
# ## https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
# ########
# 
# group1 = 3+2 =5 # Mg genes 
# group2 = 11+2 = 13 # HBC genes
# Overlap = 2. # Overlap
# total= length(merged_seurat_nolb2$gene) # All genes in genome
# 
# 
# fisher.test(matrix(c(Overlap, group2-Overlap, group1-Overlap, Total-group2-group1 +Overlap), 2, 2), alternative='greater')$p.value
# 
# 
# ## make an empty matrix to store results
# phyper_mat = matrix(0, length(non_br), length(br_clusters))
# rownames(mat) = non_br
# colnames(mat) = br_clusters
# 
# for (i in 1:length(non_br)){
#   deg_pl = merged_seurat_sig_m %>% dplyr::filter(cluster == non_br[i])
#   deg_pl = unique(deg_pl$gene)
#   
#   for (j in 1:length(br_clusters)){
#     deg_br = merged_seurat_sig_m %>% dplyr::filter(cluster == br_clusters[j])
#     deg_br = unique(deg_br$gene)
#     overlap = intersect(deg_br,deg_pl)
#     
#     group1 = length(deg_br) # Mg genes 
#     group2 = length(deg_pl) # HBC genes
#     Overlap = length(overlap(deg_br))
#     total= length(merged_seurat_nolb2$gene) # All genes in genome
#     
#     
#     fisher.test(matrix(c(Overlap, group2-Overlap, group1-Overlap, Total-group2-group1 +Overlap), 2, 2), alternative='greater')$p.value
#     
#     mat[i,j] = fraction
#   }
# }
# library(pheatmap)
# pheatmap(mat)
# view(mat)
# 
# 
# 
# ##########
# ## 5. Look at the overlap in genes between Mg_YSI_Spp1 and 2 HBC populations
# ########
# 
# hbc_check=c('HBC_Pf4','HBC_Cd72')
# mg_check=c('Mg_YSI_Spp1')
# 
# list <- vector(mode="list", length=length(hbc_check))
# names(list) <- hbc_check
# 
# for (i in 1:length(hbc_check)){
#   deg_pl = merged_seurat_sig_m %>% dplyr::filter(cluster == hbc_check[i])
#   deg_pl = unique(deg_pl$gene)
#   
#   for (j in 1:length(mg_check)){
#     deg_br = merged_seurat_sig_m %>% dplyr::filter(cluster == mg_check[j])
#     deg_br = unique(deg_br$gene)
#     list[[hbc_check[i]]] <- intersect(deg_br,deg_pl)
# 
#   }
# }
# 
# print(list)
# df <- stack(list)
# view(df)
# ck<- compareCluster(geneCluster = values ~ ind,
#                                 data = df,
#                                 OrgDb = org.Mm.eg.db,
#                                 keyType="SYMBOL",
#                                 fun = "enrichGO",
#                                 ont="BP")
# 
# 
# p = clusterProfiler::dotplot(ck, showCategory =10, 
#                              font.size = 12, 
#                              includeAll = F) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# 
# ###
# ## 6. Look for terms that are only sig in Mg and HBC, not in PAMM
# ###
# 
# 
# id_ck_brain_m = unique(ck_brain_m@compareClusterResult$ID)
# id_ck_hbc_m = unique(ck_hbc_m@compareClusterResult$ID)
# id_ck_pamm_m = unique(ck_pamm_m@compareClusterResult$ID)
# 
# 
# ck_brain_hbc_not_pamm_m = ck_brain_m %>% dplyr::filter(ID %in% id_ck_hbc_m) %>% dplyr::filter(!(ID %in% id_ck_pamm_m))
# 
# id_to_plot_m = ck_brain_hbc_not_pamm_m@compareClusterResult$ID
# ck_to_plot_m = ck_m_filter %>% dplyr::filter(ID %in% id_to_plot_m) %>% dplyr::filter(!(cluster == "Mono_FBr"))
# head(ck_to_plot_m)
# 
# p = clusterProfiler::dotplot(ck_to_plot_m, x =~cluster,showCategory =10, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# 
# ##
# ## Instead only look at Genes that dont overla
# ##
# # ck_overlap = readRDS(file = paste0(outdir, "ck_cluster_fix35_latentpair_latentb2_padj_0.05_lfc_0.25_male_mg_and_hbc_not_pamm_no_ccl8.rds"))
# # 
# # p = clusterProfiler::dotplot(ck_overlap, x =~cluster,showCategory =3, font.size = 12) +
# #   theme(axis.text.x=element_text(angle=45, hjust=1))
# # 
# # show(p)
# 
# 
# ## Brain
# 
# ck_m_brain = ck_m %>% dplyr::filter(cluster %in% br_clusters)
# 
# p = clusterProfiler::dotplot(ck_m_brain, x =~cluster,showCategory =5, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# 
# 
# ck_cc=readRDS(file = paste0(outdir, "ck_cc_indiv_cluster_fix35_latentpair_latentb2_padj_0.05_lfc_0.25_nodir.rds"))
# p = clusterProfiler::dotplot(ck_cc, x =~cluster,showCategory =2, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# write.csv(ck, paste0(outdir, "ck_cluster_fix35_latentpair_latentb2_padj_0.05_lfc_0.25_nodir.csv"))
# p = clusterProfiler::dotplot(ck, x =~cluster,showCategory =1, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# ck_dir=readRDS(file = paste0(outdir, "ck_cluster_fix35_latentpair_latentb2_padj_0.05_lfc_0.25_withdir.rds"))
# p = clusterProfiler::dotplot(ck_dir, x =~cluster,showCategory =1, font.size = 12) +
#   facet_grid(~type_direction) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# 
# ck_select = ck %>% dplyr::filter(cluster %in% c(pl_hbc, pamm_select_clusters, br_select_clusters))
# ck_dir_select = ck_dir %>% dplyr::filter(cluster %in% c(pl_hbc, pamm_select_clusters, br_select_clusters))
# 
# p = clusterProfiler::dotplot(ck_select, x =~cluster,showCategory =2, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# show(p)
# 
# 
# ## male brain only
# 
# unique(ck@compareClusterResult$sex)
# p = clusterProfiler::dotplot(ck_male_brain, x =~cluster,showCategory =5, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# term_ad = c("adhesion")
# ck_ad = ck_select %>% dplyr::filter(grepl(paste(term_ad, collapse='|' ),Description))
# 
# p = clusterProfiler::dotplot(ck_ad, x =~cluster,showCategory =3, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# term_mig = c("chemotaxis","migration")
# 
# 
# term_immune=c("defense", "inflammatory", "interferon", "immune", "cytokine", "interleukin", "tumor", "viral", "antigen")
# 
# 
# ck_immune = ck_select %>% dplyr::filter(grepl(paste(term_immune, collapse='|' ),Description))
# 
# write.csv(ck_immune, paste0(outdir, "ck_cluster_fix35_latentpair_latentb2_padj_0.05_lfc_0.25_nodir_immune_terms.csv"))
# 
# p = clusterProfiler::dotplot(ck_immune, x =~cluster,showCategory =3, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# ck_dir_immune = ck_dir_select %>% dplyr::filter(grepl(paste(term_immune, collapse='|' ),Description))
# p = clusterProfiler::dotplot(ck_dir_immune, x =~cluster,showCategory =2, font.size = 12) +
#   facet_grid(~type_direction) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# 
# term_mig = c("chemotaxis","migration")
# 
# 
# ck_mig = ck_select %>% dplyr::filter(grepl(paste(term_mig, collapse='|' ),Description))
# 
# p = clusterProfiler::dotplot(ck_mig, x =~cluster,showCategory =3, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# 
# term_pha = c("phagocytosis")
# ck_pha = ck_select %>% dplyr::filter(grepl(paste(term_pha, collapse='|' ),Description))
# 
# p = clusterProfiler::dotplot(ck_pha, x =~cluster,showCategory =3, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# 
# term_glu = c("glucose","insulin","glycolytic")
# term_glu = ck_select %>% dplyr::filter(grepl(paste(term_glu, collapse='|' ),Description))
# 
# p = clusterProfiler::dotplot(term_glu, x =~cluster,showCategory =3, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# term_meta = c("metabol")
# term_meta = ck_select %>% dplyr::filter(grepl(paste(term_meta, collapse='|' ),Description))
# 
# p = clusterProfiler::dotplot(term_meta, x =~cluster,showCategory =3, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# 
# term_heat = c("heat","folding","shock")
# term_heat = ck_select %>% dplyr::filter(grepl(paste(term_heat, collapse='|' ),Description)) %>% dplyr::filter(!grepl("ensheathment",Description))
# 
# p = clusterProfiler::dotplot(term_heat, x =~cluster,showCategory =3, font.size = 12) +
#   facet_grid(~sex) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# print(p)
# 
# pl_select_clusters=c("HBC_Pf4",
#                      "HBC_Cd72")
# 
# br_select_clusters=c("Mg_YSI_Ccl24",
#                      "Mg_Sparc",
#                      "Mg_Spp1",
#                      "Mg_Ccl5")
# 
# pamm_select_clusters=c("PAMM_Spp1",
#                        "PAMM_Ly6c2",
#                        "PAMM_Chil3",
#                        "PAMM_Ccl8")
# 
# 
# ck_large = ck %>% dplyr::filter(cluster %in% c(br_select_clusters, pl_select_clusters, pamm_select_clusters))
# 
# p = clusterProfiler::dotplot(ck_large, x =~cluster,showCategory =1, font.size = 12) +
#   facet_grid(~type_direction) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# 
# print(p)
# 

## n


# 
# ## Add the direction of the fold change and make a composite "type_direction" so that up/down and male/female can be plotted together
# merged_seurat_brain = merged_seurat_sig %>% dplyr::filter(sample == "Brain")
# 
# ## Is the gene a ribosomal gene
# # merged_seurat_brain$rib = ifelse(grepl("^Rps",merged_seurat_brain$gene) | grepl("^Rpl",merged_seurat_brain$gene), 1,0)
# #table(merged_seurat_brain$rib, merged_seurat_brain$type_direction)
# #
# # ## First analyze the ribosomal things and understand what is going on
# # merged_seurat_brain_rib = merged_seurat_brain %>% dplyr::filter(grepl("^Rps",gene)| grepl("^Rpl",gene))
# # 
# # ## clusterprofiler
# # ck_brain_rib_dir<- compareCluster(geneCluster = gene~cluster + type_direction,
# #                           data = merged_seurat_brain_rib,
# #                           OrgDb = org.Mm.eg.db,
# #                           keyType="SYMBOL",
# #                           fun = "enrichGO",
# #                           ont="BP")
# # 
# # save(ck_brain_rib_dir, file = paste0(outdir,"ck_brain_rib_dir.Rdata"))
# # 
# # write.csv(ck_brain_rib_dir, paste0(outdir, "ck_brain_rib_dir_latentpair_padj_0.05.csv"), row.names = F)
# # 
# # p = clusterProfiler::dotplot(ck_brain_rib_dir, x =~cluster,showCategory =2, font.size = 12) + 
# #   facet_grid(~type_direction) + 
# #   theme(axis.text.x=element_text(angle=45, hjust=1))
# # 
# # dpi = 300
# # png(file=paste0(outdir, "ck_brain_rib_dir_latentpair_padj_0.05.png"), width = dpi*16, height = dpi*10, units = "px",res = dpi,type='cairo')
# # print(p)
# # dev.off()
# 
# 
# ## remove ribosomal proteins
# merged_seurat_brain_norib = merged_seurat_brain %>% dplyr::filter(!grepl("^Rps",gene)) %>% dplyr::filter(!grepl("^Rpl",gene))
# 
# # ck_brain_norib_dir<- compareCluster(geneCluster = gene~ cluster + type_direction,
# #                           data = merged_seurat_brain_norib,
# #                           OrgDb = org.Mm.eg.db,
# #                           keyType="SYMBOL",
# #                           fun = "enrichGO",
# #                           ont="BP")
# # save(ck_brain_norib_dir, file = paste0(outdir,"ck_brain_norib_dir.Rdata"))
# # 
# # write.csv(ck_brain_norib_dir, paste0(outdir, "ck_brain_norib_dir_padj_0.05.csv"), row.names = F)
# # 
# # p = clusterProfiler::dotplot(ck_brain_norib_dir, x =~cluster,showCategory =1, font.size = 12) + 
# #   facet_grid(~type_direction) + 
# #   theme(axis.text.x=element_text(angle=45, hjust=1))
# # 
# # show(p)
# # dpi = 300
# # png(file=paste0(outdir, "ck_brain_norib_dir_padj_0.05.png"), width = dpi*20, height = dpi*9, units = "px",res = dpi,type='cairo')
# # print(p)
# # dev.off()
# # 
# # ## PLACENTA
# # 
# merged_seurat_placenta = merged_seurat_sig %>% dplyr::filter(sample == "Placenta")
# 
# 
# # ## First analyze the ribosomal things and understand what is going on
# # merged_seurat_placenta_rib = merged_seurat_placenta %>% dplyr::filter(grepl("^Rps",gene)| grepl("^Rpl",gene))
# # merged_seurat_placenta_rib$type_direction = paste0(merged_seurat_placenta_rib$type, "_", merged_seurat_placenta_rib$direction)
# # 
# # ## clusterprofiler
# # ck_placenta_rib_dir<- compareCluster(geneCluster = gene~cluster + type_direction,
# #                                   data = merged_seurat_placenta_rib,
# #                                   OrgDb = org.Mm.eg.db,
# #                                   keyType="SYMBOL",
# #                                   fun = "enrichGO",
# #                                   ont="BP")
# # 
# # save(ck_placenta_rib_dir, file = paste0(outdir,"ck_placenta_rib_dir.Rdata"))
# # 
# # write.csv(ck_placenta_rib_dir, paste0(outdir, "ck_placenta_rib_dir_latentpair_padj_0.05.csv"), row.names = F)
# # 
# # p = clusterProfiler::dotplot(ck_placenta_rib_dir, x =~cluster,showCategory =2, font.size = 12) + 
# #   facet_grid(~type_direction) + 
# #   theme(axis.text.x=element_text(angle=45, hjust=1))
# # 
# # dpi = 300
# # png(file=paste0(outdir, "ck_placenta_rib_dir_latentpair_padj_0.05.png"), width = dpi*16, height = dpi*10, units = "px",res = dpi,type='cairo')
# # print(p)
# # dev.off()
# 
# 
# ## remove ribosomal proteins
# merged_seurat_placenta_norib = merged_seurat_placenta %>% dplyr::filter(!grepl("^Rps",gene)) %>% dplyr::filter(!grepl("^Rpl",gene))
# 
# # ck_placenta_norib_dir<- compareCluster(geneCluster = gene~ cluster + type_direction,
# #                                     data = merged_seurat_placenta_norib,
# #                                     OrgDb = org.Mm.eg.db,
# #                                     keyType="SYMBOL",
# #                                     fun = "enrichGO",
# #                                     ont="BP")
# # 
# # save(ck_placenta_norib_dir, file = paste0(outdir,"ck_placenta_norib_dir.Rdata"))
# # 
# # rm(ck_placenta_norib_dir,)
# # test = load(paste0(outdir,"ck_placenta_norib_dir.Rdata"))
# #   
# #   write.csv(ck_placenta_norib_dir, paste0(outdir, "ck_placenta_norib_dir_padj_0.05.csv"), row.names = F)
# # 
# # p = clusterProfiler::dotplot(ck_placenta_norib_dir, x =~cluster,showCategory =1, font.size = 12) + 
# #   facet_grid(~type_direction) + 
# #   theme(axis.text.x=element_text(angle=45, hjust=1))
# # 
# # show(p)
# # dpi = 300
# # png(file=paste0(outdir, "ck_placenta_norib_dir_padj_0.05.png"), width = dpi*20, height = dpi*9, units = "px",res = dpi,type='cairo')
# # print(p)
# # dev.off()
# 
# 
# ### make some individual plots
# combined_norib = rbind(merged_seurat_brain_norib, merged_seurat_placenta_norib)
# combined_norib = combined_norib %>% dplyr::filter(grepl("^Mg" , cluster) | grepl("^HBC", cluster))
# unique(combined_norib$cluster)
# 
# ck_combined_norib_dir<- compareCluster(geneCluster = gene~ cluster + type_direction,
#                                        data = combined_norib,
#                                        OrgDb = org.Mm.eg.db,
#                                        keyType="SYMBOL",
#                                        fun = "enrichGO",
#                                        ont="BP")
# 
# save(ck_combined_norib_dir, file = paste0(outdir,"ck_combined_norib_dir.Rdata"))
# 
# write.csv(ck_combined_norib_dir, paste0(outdir, "ck_combined_norib_dir_padj_0.05.csv"), row.names = F)
# 
# load(file = paste0(outdir,"ck_combined_norib_dir.Rdata"))
# 
# p = clusterProfiler::dotplot(ck_combined_norib_dir, x =~cluster,showCategory =1, font.size = 12) +
#   facet_grid(~type_direction) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# dpi = 300
# png(file=paste0(outdir, "ck_combined_norib_dir_padj_0.05.png"), width = dpi*20, height = dpi*9, units = "px",res = dpi,type='cairo')
# print(p)
# dev.off()


# Fig 3C alternative: for comparison, check abs(lfc)>0.58
# 
# ck=readRDS(file = paste0(indir, "ck_de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.58_nodir.rds"))
# padj_cut = 0.05
# ck_filter = ck %>% dplyr::filter(p.adjust < padj_cut)
# ck_filter@compareClusterResult$cluster = factor(ck_filter@compareClusterResult$cluster, levels=c(br_clusters, pl_hbc, pamm_clusters, mono_clusters))
# 
# ## Break into cell types
# ck_hbc_mg = ck_filter %>% dplyr::filter(cluster %in% c(pl_hbc, br_clusters))
# 
# # response to oxidative stress GO:0006979
# # response to lipopolysaccharide GO:0032496
# # regulation of inflammatory response GO:0050727
# # positive regulation of cytokine production GO:0001819
# # ATP metabolic process GO:0046034
# 
# ## All categories with overlap, select clusters
# 
# p = clusterProfiler::dotplot(ck_hbc_mg, x =~cluster,showCategory =5, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 3. Look at the Top 5 categories in all Brain and plot them in Mg and HBC only -> not use
# 
# ck_filter_m_top5_brain_m = ck_hbc_mg_m %>% dplyr::filter(ID %in% top_5_each_brain_m$ID)
# 
# p = clusterProfiler::dotplot(ck_filter_m_top5_brain_m, x =~cluster,showCategory =5, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# 
# show(p)
# 
# termsim2 = enrichplot::pairwise_termsim(ck_hbc_mg_m)
# 
# library(RColorBrewer)
# require(viridis)
# pal = brewer.pal(n = 8, name = "Set2")
# 
# emapplot(termsim2) + scale_color_gradient(low = "#132B43", high = "#132B43")
#####
## BP with direction
######
# 
# 
# ck_bp_dir=readRDS(file = paste0(indir, "ck_bp_de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.2_dir.rds"))
# padj_cut = 0.05
# ck_bp_dir_filter = ck_bp_dir %>% dplyr::filter(p.adjust < padj_cut)
# ck_bp_dir_filter@compareClusterResult$cluster = factor(ck_bp_dir_filter@compareClusterResult$cluster, levels=c(br_clusters, pl_hbc, pamm_clusters, mono_clusters))
# 
# ## Break into cell types
# ck_hbc_mg = ck_bp_dir_filter %>% dplyr::filter(cluster %in% c(pl_hbc, br_clusters))
# 
# ## All categories with overlap, select clusters
# 
# p = clusterProfiler::dotplot(ck_hbc_mg, x =~cluster,showCategory =5, font.size = 12) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))+ facet_grid(~direction)
# 
# show(p)
# 
# ggsave(p, filename = paste0(outdir,"top5_ck_bp_de.mfcombined.obsctr.latent_pair_sex_mast_padj_0.05_abslfc_0.2_dir.pdf"),
#        device = cairo_pdf,
#        width = 20, height = 10, 
#        units = "in")

