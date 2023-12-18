# Calculate the correlation of average gene expression between two integrated seurat objects
# Much of this code has been adapted from this https://github.com/archavan/covid-placenta/blob/master/code/02_reference-cell-atlas.Rmd

# HPC R library
LIB='/cluster/tufts/patralab/rbator01/R_libs/4.0.0'
.libPaths(c("",LIB))

suppressPackageStartupMessages({
  library(tidyverse)
  library(remotes)
  library(celldex)
  library(scRNAseq)
  library(SummarizedExperiment)
  library(SingleR)
  library(scuttle)
  library(Seurat)
  library(RColorBrewer)
  library(viridis)
})

get_avg_so <-function(so, res, string){
  Idents(object = so) <- res
  so.avg <- AverageExpression(so)
  so.avg.rna<- so.avg$RNA
  so.avg.rna = data.frame(so.avg.rna)
  colnames(so.avg.rna) = gsub("^",string,colnames(so.avg.rna))
  so.avg.rna$gene_name <- rownames(so.avg.rna)
  return(so.avg.rna)
}

corr_dat <-function(avg1, avg2, id1, id2){
  
  combined=dplyr::inner_join(avg1, avg2, by = "gene_name")
  rownames(combined) = combined$gene_name
  combined$gene_name=NULL
  cor.matrix <- cor(combined, method = 'spearman')
  
  # reorder correlation matrix based on clustering
  dd <- as.dist((1 - cor.matrix)) 
  hc <- hclust(dd, method = 'complete')
  cor.matrix <- cor.matrix[hc$order, hc$order]
  
  # reorder correlation matrix based on clustering
  cormat <- reshape2::melt(cor.matrix, na.rm = T)
  cormat$Var1_source <- sapply(strsplit(as.character(cormat$Var1), split = "_"), "[[", 1)
  cormat$Var2_source <- sapply(strsplit(as.character(cormat$Var2), split = "_"), "[[", 1)
  
  # subset
  cormat <- cormat[which(cormat$Var1_source == id1  & cormat$Var2_source == id2), ]
  
  # melt correlation matrix
  dat <- reshape2::melt(cormat, na.rm = T)
  
  return(list("hc" = hc, "dat" = dat, "cormat" = cormat))
}

# plotting function
clustAnnoPlot <- function(dat, query, reference, plot_title) {
  
  dat <- dat[which(dat$Var2_source %in% query & dat$Var1_source == reference), ]
  
  p <- ggplot(data = dat, 
              aes(Var1, Var2, fill = value)) +
    geom_tile(colour = "white") +
    scale_fill_viridis(name = "Spearman\nCorrelation") +
    geom_point(aes(Var1, Var2, alpha = top3),
               size = 1.5, shape = 19, stroke  = 0) +
    scale_alpha_manual(values = c(1, 0.5, 0.25), 
                       breaks = c(1, 2, 3),
                       name = "top3", na.value = 0) +
    coord_fixed(ratio = 1) +
    xlab("reference") +
    ylab("query") +
    
    labs(title = plot_title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 8, 
                                     hjust=1, vjust = 0.5),
          axis.text.y = element_text(size=8),
          axis.ticks.length = unit(0.15, units = c('lines')),
          legend.title = element_text(size = 10),
          axis.title = element_text(size = 8),
          plot.caption = element_text(size = 7),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.key.size = unit(0.8, units = c('lines')))
  
  return(p)
}


# Load our placenta data ---------------
pl_so = readRDS("analysis/final_rds/subset/pl_obsctr_final.rds")

avg_pl=get_avg_so(so=pl_so, res = "cluster_name", string="p_")

#Suryawanshi - they provide the average data -------------
avg_surya = read.csv("~/slonimlab/rbator01/reference_data/suryawanshi/aau4788_data_file_s1_mousenames.csv")
rownames(avg_surya) = avg_surya$X
avg_surya$gene_name = avg_surya$X
avg_surya$X = NULL

l = corr_dat(avg1 = avg_pl, avg2 = avg_surya, id1 = "surya", id2 = "p")
hc = l$hc
dat = l$dat
cormat = l$cormat
ll = c("vil.HC", "dec.MAC","dec.DC2",
       "dec.DC1","dec.TC","dec.NK1","dec.NK2",
       "dec.LEC","dec.VEC","vil.VEC","vil.FB2","vil.FB1",
       "vil.FB3","dec.SMC","dec.FB2","dec.DSC",  
       "dec.FB1","dec.EEC","vil.EB","vil.EVT","vil.VCT","vil.SCT")


anno.surya <- clustAnnoPlot_notop(dat = cormat, 
                                  query = c("p"), 
                                  reference = "surya", 
                                  plot_title="", 
                                  order = ll)

cowplot::ggsave2(anno.surya, device = "pdf", width = 6, height = 4, units = "in",
                 filename = "analysis/plots/compare_reference/surya_subset_heatmap.pdf")

# Vento-Tormo et al ---------

mouse_so = readRDS("~/slonimlab/rbator01/reference_data/vento-tormo_2018/clustered_so_mousenames.rds")

avg_vt=get_avg_so(so=mouse_so, res = "annotation",  string="vt_")

l = corr_dat(avg1 = avg_pl, avg2 = avg_vt, id1 = "vt", id2 = "p")
hc = l$hc
dat = l$dat
cormat = l$cormat

ll =c("HB","dM1","dM2","dM3" , "MO" ,"DC1","DC2", "dNK1","dNK2","dNK3","dNK.p" ,"Granulocytes", "Tcells","NK.CD16neg", "NK.CD16pos","ILC3"  ,"Plasma",
      "Endo.f" ,"Endo.m", "Endo.L" ,"fFB1"   ,"fFB2"   ,"dP1"    ,"dP2"    ,"dS3"    ,"dS1" ,"dS2","EVT"    ,"VCT"    ,"SCT","Epi1","Epi2" )

anno.vento <- clustAnnoPlot_notop(dat = cormat, query = "p", reference = "vt", plot_title="", order=ll)
show(anno.vento)
cowplot::ggsave2(anno.vento, device = "pdf", width = 8, height = 4, units = "in",
                 filename = paste0("analysis/plots/compare_reference/", "vento_subset_heatmap.pdf"))


# Lu-Culligan et al ------------

mouse_so = readRDS("~/slonimlab/rbator01/reference_data/lu_culligan/GSE171381_so_mousenames.rds")

mouse_so = readRDS("~/slonimlab/rbator01/reference_data/lu_culligan/GSE171381_so_mousenames.rds")
avg_lc=get_avg_so(so=mouse_so, res = "annotation_merged",  string="lc_")

l = corr_dat(avg1 = avg_pl, avg2 = avg_lc, id1 = "lc", id2 = "p")
hc = l$hc
dat = l$dat
cormat = l$cormat

# for annotation merged column
ll=c("vil.Hofb","Mono_1","Mono_2","APC","vil.Ery","Gran","NK_1","NK_2","NK_3","Tcell_1","Tcell_2","Tcell_3","Bcell",
     "dec.DSC","dec.FB","dec.SMC","dec.Endo","vil.EVT ","vil.VCT ","vil.SCT","vil.FB")

anno.covid <- clustAnnoPlot_notop(dat = cormat, query = c("p"), reference = "lc", plot_title="", order = ll)

show(anno.covid)
cowplot::ggsave2(anno.covid, device = "pdf", width = 6, height = 4, units = "in",
                 filename = paste0("analysis/plots/compare_reference/", "lu_c_subset_heatmap.pdf"))

