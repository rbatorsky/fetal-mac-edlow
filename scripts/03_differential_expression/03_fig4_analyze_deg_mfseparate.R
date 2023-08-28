## This is MF separate

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
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')

indir="analysis/deg_seurat/obsctr/mergesome/mfseparate/"
outdir="analysis/deg_seurat/obsctr/mergesome/mfseparate/plots/"


# fig 4A: DEG barplot -------------------------------------------------------------

indir="analysis/deg_seurat/obsctr/mergesome/mfseparate/"
outdir="analysis/deg_seurat/obsctr/mergesome/mfseparate/plots/"

merged_seurat_nolb2 = read.csv(paste0(indir, "de.obsctr_fix35.obs_vs_ctr.latent_pair_mast_logfcthresh0_23Jun22.csv"))
merged_seurat_m = merged_seurat_nolb2 %>% dplyr::filter(sex == "M")
merged_seurat_m$X = NULL

merged_seurat_lb2_f = read.csv(paste0(indir, "de.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_23Jun22.csv"))
merged_seurat_lb2_f$sex = "F"
merged_seurat_lb2_f$X = NULL

merged_seurat = rbind(merged_seurat_m, merged_seurat_lb2_f)
merged_seurat_sig = merged_seurat %>% dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC)>0.25)

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

t = as.data.frame(table(merged_seurat_sig$cluster,merged_seurat_sig$sex,merged_seurat_sig$cc1))
t$cluster = t$Var1
t$sex = t$Var2
t$type = t$Var3
t$Freq = ifelse(t$sex == "F",-t$Freq,t$Freq)

# Mg colors
mg_col=c("blue","cornflowerblue","turquoise4","deepskyblue","dodgerblue4","steelblue1","turquoise1","lightskyblue")
# BAM colors
bam_col=c("limegreen","darkgreen","palegreen","yellowgreen","green","olivedrab")
# mono colors
mono_col=c("grey10","grey25","grey40","grey55","grey70")


brain_col = tidyr::tribble(
  ~name, ~color,
  "Mg_YSI_Pf4", bam_col[1],
  "Mg_YSI_cellcycle",bam_col[2],
  "Mg_Sparc",mg_col[1],
  "Mg_Ccl5",mg_col[5],
  "Mg_Spp1",mg_col[4],
  "Mg_Hspb1",mg_col[3],
  "Mg_cellcycle",mg_col[2],
  "Mono_FBr",mono_col[1],
)
# HBC colors
hbc_col=c("plum","purple1","mediumpurple4","orchid","mediumpurple1")
# dec mac colors
dm_col=c("red","lightcoral","firebrick","tomato","brown","darkred")

pal <- colorRampPalette(c("brown","red","orange"))
dm_col=pal(6)

# mono colors
mono_col=c("grey10","grey25","grey40","grey55","grey70")

pl_col = tidyr::tribble(
  ~name, ~color,
  "HBC_Cd72",hbc_col[1],
  "HBC_Pf4",hbc_col[2],
  "HBC_cellcycle",hbc_col[3],
  "PAMM_Chil3",dm_col[1],
  "PAMM_Spp1",dm_col[2],
  "PAMM_MHCII",dm_col[3],
  "PAMM_S100a9",dm_col[5],
  "PAMM_Ccl8",dm_col[6],
  "Mono_FPl",mono_col[1],
)

levels=rev(c(deframe(brain_col[,1]), deframe(pl_col[,1])))
t$cluster = factor(t$cluster, levels=levels)

p = ggplot(data=t, aes(x=cluster, y=Freq, fill=type)) +
  geom_bar(stat="identity") +
  geom_bar(stat="identity", fill="black", aes(alpha=sex)) +
  scale_alpha_manual(values=c(0,0.2),
                     name="sex",
                     breaks=c("M","F"),         # can use to set order
                     labels=c("Male","Female")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Number of DEG") +
  xlab("")+ coord_flip() + theme_classic()+ theme(axis.text=element_text(size=12),
                                                  axis.title=element_text(size=14))

show(p)

ggsave(p, filename = paste0(outdir,"obsvsctr_m_compare_f_latent_pair_lb2_deg_padj_0.05_abslfc_0.25_barplot.pdf"),
       device = cairo_pdf,
       width = 5, height = 7,
       units = "in")


# Statistical Test

merged_seurat_labelsig = merged_seurat %>% mutate(sig = ifelse(p_val_adj<0.05 & abs(avg_log2FC)>0.25, "sig","nosig"))

c = c()
p = numeric()
nf=numeric()
nm=numeric()

for (cl in unique(merged_seurat_labelsig$cluster)){
  test = merged_seurat_labelsig %>% dplyr::filter(cluster == cl)
  table = as.data.frame.matrix(table(test$sig, test$sex))
  print(cl)
  print(table)
  result = chisq.test(test$sig, test$sex)
  c = c(c,cl)
  p = c(p, result$p.value)
  nf=c(nf, table[2,1])
  nm=c(nm, table[2,2])
  
}

df = data.frame(cluster= c, pval=p, nf = nf, nm=nm)
view(df)


# Fig 4B sex different or consistent -----

ck = readRDS(paste0(outdir,"cp_hbc_mgysi_sexconsistency.rds"))

p = clusterProfiler::dotplot(ck,showCategory =5, by = "count", font.size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
show(p)
ggsave(p, filename = paste0(outdir,"go_sankey_hbc_mgysi.pdf"),
       device = cairo_pdf,
       width = 7, height = 6, 
       units = "in")

# do the ordering
select_terms = c('response to lipopolysaccharide',
                 'response to molecule of bacterial origin',
                 'response to interferon-alpha',
                 'myeloid cell differentiation',
                 'myeloid cell homeostasis',
                 'regulation of hemopoiesis',
                 'negative regulation of cell-cell adhesion',
                 "'de novo' protein folding",                                                                                                     
                 'chaperone-mediated protein folding',
                 'striated muscle tissue development',
                 'neuron death')


ck = ck %>% 
  dplyr::filter(Description %in% select_terms)

ck@compareClusterResult$Description = factor(ck@compareClusterResult$Description, levels=rev(select_terms))
## The below is to order the factors
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

object = ck@compareClusterResult
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

p <- ggplot(df, aes_string(x = 'Cluster', y = "Description", size = 'Count'))      

p = p +  geom_point(aes_string(color = colorBy)) + 
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) + ggtitle(title) + 
  DOSE::theme_dose(font.size) + 
  scale_size_continuous(range=c(3, 8)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
show(p)
dev.off()
ggsave(p, filename = paste0(outdir,"sex_consistent_go.pdf"),
       device = cairo_pdf,
       width = 7, height = 6, 
       units = "in")


# ck = pairwise_termsim(ck)
# p1 = emapplot(ck, 
#               show=30,
#               legend_n = 2)+ theme(text=element_text(size=6, family="Arial"))
# 
# show(p1)
# ggsave(p1, filename = paste0(outdir,"3b_pf4_10_emmapplot_deg_mfcombined_latent_pair_lb2_sex.pdf"),
#        device = cairo_pdf,
#        width = 5, height = 4.5, 
#        units = "in")


# Fig 4C - IPA sex different ----

# Fig 4D - Single category IPA plot ----

