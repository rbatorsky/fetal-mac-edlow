# Fig 2D - This script calculates the correlation of average gene expression between two integrated seurat objects
# Much of this code has been adapted from this https://github.com/archavan/covid-placenta/blob/master/code/02_reference-cell-atlas.Rmd

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
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')
outdir ="analysis/plots/"


# read in the data ----
br_so = readRDS("analysis/final_rds/subset/br_obsctr_final.rds")
pl_so = readRDS("analysis/final_rds/subset/pl_obsctr_final.rds")

# functions ----
get_avg <-function(infile, res, string){
  so = readRDS(infile)
  Idents(object = so) <- res
  so.avg <- AverageExpression(so)
  so.avg.rna<- so.avg$RNA
  so.avg.rna = data.frame(so.avg.rna)
  colnames(so.avg.rna) = gsub("^",string,colnames(so.avg.rna))
  so.avg.rna$gene_name <- rownames(so.avg.rna)
  return(so.avg.rna)
}

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

# analysis ----
avg_br=get_avg_so(so=seurat_integrated_br_rename, res = "cluster_name",  string="b_")
avg_pl=get_avg_so(so=seurat_integrated_pl_rename, res = "cluster_name", string="p_")

colnames(avg_pl)

rm(dat)
rm(l)
l = corr_dat(avg1 = avg_pl, avg2 = avg_br, id1 = "p", id2 = "b")
hc = l$hc
dat = l$dat
cormat = l$cormat

# top 3 matches with highest correlation coefficients
## max correlation coefficients
dat$top3 <- NA
for(i in unique(dat$Var2)){
  ind <- which(dat$Var2 == i & dat$Var1_source != strsplit(i, "_")[[1]][1])
  val <- dat$value[ind]
  top3ind <- ind[order(val, decreasing = TRUE)[1:3]]
  dat$top3[top3ind[1]] <- "1"
  # dat$top3[top3ind[2]] <- "2"
  # dat$top3[top3ind[3]] <- "3"
}


dat$Var1 = gsub("pl_", "", dat$Var1)
dat$Var2 = gsub("br_", "", dat$Var2)

## heatmap plot 
dpi = 300

dat$Var1 = gsub('p_','',dat$Var1)
dat$Var2 = gsub('b_','',dat$Var2)

#png(file=paste0(home, "/plots/br_pl_compare/heatmap_compare.png"), width = dpi*4, height = dpi*3, units = "px",res = dpi,type='cairo')
p <- ggplot(data = dat, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white") +
  scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white") + 
  geom_point(aes(Var1, Var2, alpha = top3),
             size = 2, shape = 19, stroke  = 0) +
  scale_alpha_manual(values = c(1, 0.5, 0.25),
                     breaks = c(1, 2, 3),
                     name = "", na.value = 0) +
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
        ylab("Brain Cluster") +
        xlab("Placenta Cluster")


ggsave(p, filename = paste0(home, "/plots/br_pl_compare/heatmap_compare.pdf"),
       device = cairo_pdf,
       width = 5, height = 4, 
       units = "in")

dev.off()


# Heatmap, not used ----
dpi = 300
png(file=paste0(home, "/plots//br_subset_heatmap_0.4_0.2.png"), width = dpi*4, height = dpi*3, units = "px",res = dpi,type='cairo')
p <- ggplot(data = dat, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = 'white', high = 'blue',
                      name = "Spearman\nCorrelation") +
  geom_point(aes(Var1, Var2, alpha = top3),
             size = 1.5, shape = 19, stroke  = 0) +
  scale_alpha_manual(values = c(1, 0.5, 0.25),
                     breaks = c(1, 2, 3),
                     name = "top3", na.value = 0) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 8,
                                   hjust=1, vjust = 0.5),
        axis.text.y = element_text(size=8),
        axis.ticks.length = unit(0.15, units = c('lines')),
        legend.title = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.size = unit(0.8, units = c('lines')))


# Dendogram ----
dpi = 300
png(file="br_subset_hc.png", width = dpi*6, height = dpi*12, units = "px",res = dpi,type='cairo')
hcd=as.dendrogram(hc)
par(mar = c(2,2,2,10)) 
p=plot(hcd,horiz=T)
print(p)
dev.off()
