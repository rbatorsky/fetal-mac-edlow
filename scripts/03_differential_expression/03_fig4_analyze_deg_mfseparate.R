# MF separate - this contains the complete code, but some functions can be written and moved to a utility startup script

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
  library(ComplexHeatmap)
})

setwd('/cluster/tufts/slonimlab/rbator01/fetal-mac-edlow/')

indir="analysis/deg_seurat/obsctr/mergesome/mfseparate/"
outdir="analysis/deg_seurat/obsctr/mergesome/mfseparate/plots/"


# Fig 4A: DEG barplot -------------------------------------------------------------

merged_seurat_nolb2 = read.csv(paste0(indir, "de.obsctr_fix35.obs_vs_ctr.latent_pair_mast_logfcthresh0_23Jun22.csv"))
merged_seurat_m = merged_seurat_nolb2 %>% dplyr::filter(sex == "M")
merged_seurat_m$X = NULL

merged_seurat_lb2_f = read.csv(paste0(indir, "de.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_23Jun22.csv"))
merged_seurat_lb2_f$sex = "F"
merged_seurat_lb2_f$X = NULL

merged_seurat = rbind(merged_seurat_m, merged_seurat_lb2_f)
deg_indiv = merged_seurat %>% 
  dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC)>0.25)

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
deg_indiv$cc1 = ifelse(grepl('^Mg_YSI',deg_indiv$cluster),'MgYSI',
                               ifelse(grepl('^Mg_',deg_indiv$cluster),'Mg',
                                      ifelse(grepl('^HBC',deg_indiv$cluster),'HBC',
                                             ifelse(grepl('^PAMM',deg_indiv$cluster),'PAMM',
                                                    ifelse(grepl('^Mono', deg_indiv$cluster),'Mono','NA')))))

t = as.data.frame(table(deg_indiv$cluster,deg_indiv$sex,deg_indiv$cc1))
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
  ylab("Number of deg_indiv") +
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

# Fig 4B sex different or consistent compare between two clusters -----

merged_seurat = read.csv(paste0(indir, "de.mfseparate.cc.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_maleandfemale.csv"))
merged_seurat$dir = ifelse(merged_seurat$avg_log2FC > 0, "up","dn")

cl1="MgYSI"
cl2="HBC"

deg_combined = merged_seurat %>% 
  dplyr::filter(p_val_adj<0.05 & abs(avg_log2FC)>0.3) %>%
  dplyr::filter(cluster %in% c(cl1,cl2)) %>%
  mutate(sex_dir = paste0(sex,"_",dir)) %>%
  dplyr::select(gene, cluster, sex_dir)

head(deg_combined)
# check for different or consistent directions in gene changes in all the DEG
deg_sexdir = merged_seurat %>% 
  dplyr::filter(gene %in% deg_combined$gene) %>%
  dplyr::filter(cluster %in% c(cl1,cl2)) %>%
  dplyr::select(gene, cluster, sex, dir ) %>%
  pivot_wider(names_from = sex, values_from = dir) %>%
  mutate(mf_agree = ifelse(M == F, 'sex_consistent','sex_different')) %>% 
  dplyr::filter(!(is.na(mf_agree)))

table(deg_sexdir$cluster,deg_sexdir$mf_agree)

# ck <- compareCluster(geneCluster = gene ~ mf_agree + cluster ,
#                      data = deg_sexdir,
#                      OrgDb = org.Mm.eg.db,
#                      keyType="SYMBOL",
#                      fun = "enrichGO",
#                      ont="BP" )
# 
# saveRDS(ck, paste0(outdir,"cp_hbc_mgysi_sexconsistency.rds"))
# write.csv(ck, paste0(outdir,"cp_hbc_mgysi_sexconsistency.csv"))

ck = readRDS(paste0(outdir,"cp_hbc_mgysi_sexconsistency.rds"))

p = clusterProfiler::dotplot(ck,showCategory =5, by = "count", font.size = 12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
show(p)
ggsave(p, filename = paste0(outdir,"go_sankey_hbc_mgysi.pdf"),
       device = cairo_pdf,
       width = 7, height = 6, 
       units = "in")

# order these terms
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

## Canonical Pathways
cp_scores = read.csv("analysis/ipa/mf_separate_combined_cluster/sex_stratified_combinecluster_ipa_cp.csv", 
                     header=T, na.strings ="N/A")

overall_scores_long = cp_scores %>% 
  gather(key=cluster, value="z.score", -Canonical.Pathways) %>%
  mutate(cluster_sex = cluster) %>%
  separate(cluster, c("cluster","sex"), sep="_") 

overall_scores_long %>% dplyr::filter(Canonical.Pathways == "EIF2 Signaling")


cl="MgYSI"
n_select=20
f_select_up = overall_scores_long %>% dplyr::filter(sex == "F" & z.score > 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(-z.score) %>% head(n_select) 
m_select_up = overall_scores_long %>% dplyr::filter(sex == "M" & z.score > 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(-z.score) %>% head(n_select)
f_select_dn = overall_scores_long %>% dplyr::filter(sex == "F" & z.score < 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(z.score) %>% head(n_select) 
m_select_dn = overall_scores_long %>% dplyr::filter(sex == "M" & z.score < 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(z.score) %>% head(n_select)

select = rbind(f_select_up, m_select_up, f_select_dn, m_select_dn)
des = unique(select$Canonical.Pathways)

'Neuroinflammation Signaling Pathway' %in% des

## Filter the main data frame for only the descriptions in the top n of either 
f_select = overall_scores_long  %>% dplyr::filter(sex == "F") %>% dplyr::filter(cluster == cl & Canonical.Pathways %in% des)
m_select = overall_scores_long  %>% dplyr::filter(sex == "M") %>% dplyr::filter(cluster == cl & Canonical.Pathways %in% des) 

select_all = rbind(f_select, m_select) 

z_score_sum = m_select %>% 
  dplyr::rename(m.z.score = z.score) %>% 
  dplyr::select(c(Canonical.Pathways, sex, m.z.score)) %>% 
  full_join(f_select %>% 
              dplyr::rename(f.z.score = z.score)  %>% 
              dplyr::select(c(Canonical.Pathways, sex, f.z.score)), by="Canonical.Pathways") %>%
  replace_na(replace = list(m.z.score=0,f.z.score=0)) %>%
  mutate(z.score.sum = m.z.score + f.z.score) %>%
  dplyr::select(c(Canonical.Pathways,z.score.sum)) %>% 
  dplyr::full_join(select_all, by="Canonical.Pathways")

mg_select_all = z_score_sum
p = ggplot(data=mg_select_all, aes(x=reorder(Canonical.Pathways,z.score.sum), 
                                   y=z.score, 
                                   fill=sex,
                                   color=sex)) +
  geom_bar(stat="identity",colour="grey",size=0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(aes(yintercept = 0)) +
  ylab("z.score") +
  xlab("")+ coord_flip() + 
  theme_minimal() + 
  ggtitle(cl) 

show(p)
ggsave(p, filename = paste0(outdir,paste0(cl,"_top10_mf_barplot_combinecluster_cp_df.pdf")),
       device = cairo_pdf,
       width = 8, height = 5,
       units = "in")

cl="HBC"
n_select=50
f_select_up = overall_scores_long %>% dplyr::filter(sex == "F" & z.score > 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(-z.score) %>% head(n_select) 
m_select_up = overall_scores_long %>% dplyr::filter(sex == "M" & z.score > 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(-z.score) %>% head(n_select)
f_select_dn = overall_scores_long %>% dplyr::filter(sex == "F" & z.score < 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(z.score) %>% head(n_select) 
m_select_dn = overall_scores_long %>% dplyr::filter(sex == "M" & z.score < 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(z.score) %>% head(n_select)

f_select_dn

select = rbind(f_select_up, m_select_up, f_select_dn, m_select_dn)

des = unique(select$Canonical.Pathways)
'Neuroinflammation Signaling Pathway' %in% des


## Filter the main data frame for only the descriptions in the top 10 of either 
f_select = overall_scores_long  %>% dplyr::filter(sex == "F") %>% dplyr::filter(cluster == cl & Canonical.Pathways %in% des)
m_select = overall_scores_long  %>% dplyr::filter(sex == "M") %>% dplyr::filter(cluster == cl & Canonical.Pathways %in% des) 

select_all = rbind(f_select, m_select) 

library(dplyr)
z_score_sum = m_select %>% 
  dplyr::rename(m.z.score = z.score) %>% 
  dplyr::select(c(Canonical.Pathways, sex, m.z.score)) %>% 
  full_join(f_select %>% 
              dplyr::rename(f.z.score = z.score)  %>% 
              dplyr::select(c(Canonical.Pathways, sex, f.z.score)), by="Canonical.Pathways") %>%
  replace_na(replace = list(m.z.score=0,f.z.score=0)) %>%
  mutate(z.score.sum = m.z.score + f.z.score) %>%
  dplyr::select(c(Canonical.Pathways,z.score.sum)) %>% 
  dplyr::full_join(select_all, by="Canonical.Pathways")

hbc_select_all = z_score_sum

p = ggplot(data=hbc_select_all, aes(x=reorder(Canonical.Pathways,z.score.sum), 
                                    y=z.score, 
                                    fill=sex,
                                    color=sex)) +
  geom_bar(stat="identity",colour="grey",size=0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(aes(yintercept = 0)) +
  ylab("z.score") +
  xlab("")+ coord_flip() + 
  theme_minimal() + 
  ggtitle(cl) 

show(p)
ggsave(p, filename = paste0(outdir,paste0(cl,"_top10_mf_barplot_combinecluster_cp_df.pdf")),
       device = cairo_pdf,
       width = 8, height = 5,
       units = "in")

#%>%
#  dplyr::filter(abs(z.score) > 0.75)

#select pathways that are the most different from male to female, then those that overlap Mg/HBC

library(dplyr)

cl="HBC"
select_all_diff_hbc = overall_scores_long %>%
  dplyr::filter(cluster == cl ) %>%
  dplyr::select(-c(cluster,sex)) %>%
  spread(cluster_sex, z.score, fill = 0) %>%
  mutate(diff = (abs(get(paste0(cl, "_M")) - get(paste0(cl, "_F"))))/get(paste0(cl, "_M"))) %>%
  dplyr::filter(abs(diff) > .5) %>%
  dplyr::rename(!!paste0(cl, "_diff") := diff)

cl="MgYSI"
select_all_diff_mg = overall_scores_long %>%
  dplyr::filter(cluster == cl ) %>%
  dplyr::select(-c(cluster,sex)) %>%
  spread(cluster_sex, z.score, fill = 0) %>%
  mutate(diff = (abs(get(paste0(cl, "_M")) - get(paste0(cl, "_F"))))/get(paste0(cl, "_M"))) %>%
  dplyr::filter(abs(diff) > .5) %>%
  dplyr::rename(!!paste0(cl, "_diff") := diff)

## select the ones in both HBC and Mg, removing the NA in each
merge_diff = select_all_diff_hbc %>%
  dplyr::full_join(select_all_diff_mg, by="Canonical.Pathways") %>% 
  dplyr::filter(!is.na(HBC_diff)) %>%
  dplyr::filter(!is.na(MgYSI_diff)) 

# select HBC
select_all_diff_hbc = merge_diff %>%
  dplyr::select(c('Canonical.Pathways','HBC_F','HBC_M'))%>% 
  gather(key=cluster_sex, value="z.score", -Canonical.Pathways) %>%
  separate(cluster_sex, c("cluster","sex"), sep="_") 

## select mgysi
select_all_diff_mg = merge_diff %>%
  dplyr::select(c('Canonical.Pathways','MgYSI_F','MgYSI_M'))%>% 
  gather(key=cluster_sex, value="z.score", -Canonical.Pathways) %>%
  separate(cluster_sex, c("cluster","sex"), sep="_")%>% 
  arrange(z.score)

## arrange levels
lvl = unique(select_all_diff_mg$Canonical.Pathways)

select_all_diff_mg = select_all_diff_mg %>%
  mutate(Canonical.Pathways = factor(Canonical.Pathways, levels=lvl))

select_all_diff_hbc = select_all_diff_hbc %>%
  mutate(Canonical.Pathways = factor(Canonical.Pathways, levels=lvl))


## HBC plot
p = ggplot(data=select_all_diff_hbc, aes(x=Canonical.Pathways,
                                         y=z.score,
                                         fill=sex,
                                         color=sex)) +
  geom_bar(stat="identity",colour="grey",size=0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(aes(yintercept = 0)) +
  ylab("z.score") +
  xlab("")+ coord_flip() +
  theme_minimal() +
  ggtitle('HBC')

show(p)

ggsave(p, filename = "",
       device = cairo_pdf,
       width = 8, height = 5,
       units = "in")



p = ggplot(data=select_all_diff_mg, aes(x=Canonical.Pathways,
                                        y=z.score,
                                        fill=sex,
                                        color=sex)) +
  geom_bar(stat="identity",colour="grey",size=0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(aes(yintercept = 0)) +
  ylab("z.score") +
  xlab("")+ coord_flip() +
  theme_minimal() +
  ggtitle('MgYSI')

show(p)


ggsave(p, filename = "",
       device = cairo_pdf,
       width = 8, height = 5,
       units = "in")




#Find the ones that are in both

mg_hbc_select_all = mg_select_all %>% 
  dplyr::inner_join(hbc_select_all, by = "Canonical.Pathways")

cat_to_show = unique(mg_hbc_select_all$Canonical.Pathways)

cat_to_show

#Now plot Mg

cl="MgYSI"
mg_select_all_both = mg_select_all %>%
  dplyr::filter(Canonical.Pathways %in% cat_to_show)

p = ggplot(data=mg_select_all_both, aes(x=reorder(Canonical.Pathways,z.score.sum), 
                                        y=z.score, 
                                        fill=sex,
                                        color=sex)) +
  geom_bar(stat="identity",colour="grey",size=0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(aes(yintercept = 0)) +
  ylab("z.score") +
  xlab("")+ coord_flip() + 
  theme_minimal() + 
  ggtitle(cl) 

show(p)

ggsave(p, filename = paste0(outdir,"MgYSI_ipa_cp_mg_hbc_overlap.pdf"),
       device = cairo_pdf,
       width = 7, height = 3.5,
       units = "in")



ll = mg_select_all_both %>%
  dplyr::select(c(Canonical.Pathways, z.score.sum)) %>%
  distinct() %>%
  arrange(-z.score.sum) %>% 
  dplyr::rename(z.score.sum.mg = z.score.sum)

ll

cl="HBC"

hbc_select_all_both = hbc_select_all %>%
  dplyr::filter(Canonical.Pathways %in% cat_to_show) %>%
  left_join(ll, by="Canonical.Pathways")

hbc_select_all_both
p = ggplot(data=hbc_select_all_both, aes(x=reorder(Canonical.Pathways,z.score.sum.mg), 
                                         y=z.score, 
                                         fill=sex,
                                         color=sex)) +
  geom_bar(stat="identity",colour="grey",size=0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(aes(yintercept = 0)) +
  ylab("z.score") +
  xlab("")+ coord_flip() + 
  theme_minimal() + 
  ggtitle(cl) 

show(p)

ggsave(p, filename = paste0(outdir,"HBC_ipa_cp_mg_hbc_overlap.pdf"),
       device = cairo_pdf,
       width = 7, height = 3.5,
       units = "in")


# select pathways that are the most different from male to female, then sort by zscore

for (cl in c('HBC','Mg','MgYSI','PAMM')){
  select_all_diff = overall_scores_long %>% 
    dplyr::filter(cluster == cl) %>% 
    mutate(cluster_sex = paste0(cluster,"_",sex)) %>%
    dplyr::select(-c(cluster,sex)) %>%
    spread(cluster_sex, z.score) %>%
    mutate(diff = (abs(get(paste0(cl, "_M")) - get(paste0(cl, "_F"))))/get(paste0(cl, "_M"))) %>%
    dplyr::filter(abs(diff) > .5)  %>% 
    dplyr::select(-c(diff)) %>%
    gather(key=cluster_sex, value="z.score", -Canonical.Pathways) %>%
    separate(cluster_sex, c("cluster","sex"), sep="_") 
  
  nhead=8
  f_select_up = select_all_diff %>% dplyr::filter(sex == "F" & z.score > 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(-z.score) %>% head(nhead) 
  m_select_up = select_all_diff %>% dplyr::filter(sex == "M"& z.score > 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(-z.score) %>% head(nhead)
  f_select_dn = select_all_diff %>% dplyr::filter(sex == "F"& z.score < 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(z.score) %>% head(nhead) 
  m_select_dn = select_all_diff %>% dplyr::filter(sex == "M"& z.score < 0)%>% dplyr::filter(cluster == cl) %>% dplyr::arrange(z.score) %>% head(nhead)
  
  select = rbind(f_select_up, m_select_up, f_select_dn, m_select_dn)
  
  des = unique(select$Canonical.Pathways)
  
  ## Filter the main data frame for only the descriptions in the top 10 of either 
  f_select = overall_scores_long  %>% dplyr::filter(sex == "F") %>% dplyr::filter(cluster == cl & Canonical.Pathways %in% des)
  m_select = overall_scores_long  %>% dplyr::filter(sex == "M") %>% dplyr::filter(cluster == cl & Canonical.Pathways %in% des) 
  select_all = rbind(f_select, m_select) 
  
  p = ggplot(data=select_all, aes(x=reorder(Canonical.Pathways,z.score), 
                                  y=z.score, 
                                  fill=sex,
                                  color=sex)) +
    geom_bar(stat="identity",colour="grey",size=0.25) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    geom_hline(aes(yintercept = 0)) +
    ylab("z.score") +
    xlab("")+ coord_flip() + 
    theme_minimal() + 
    ggtitle(cl) 
  
  show(p)
  
  
  ggsave(p, filename = paste0(outdir,paste0(cl,"_top8_50percentdiff_zscoresort_mf_barplot_combinecluster.pdf")),
         device = cairo_pdf,
         width = 8, height = 5,
         units = "in")
  
}

show(p)


# This is code for Lydia's select pathways

# Diseases and Functions
df_scores = read.csv("~/analysis/ipa/mf_separate_combined_cluster/Diseases and Biofunctions M and F sep cell type_rebeccaedit.csv",
                     header=T, na.strings ="N/A")

df_scores[is.na(df_scores)] <- 0
df_scores$Canonical.Pathways = df_scores$Diseases.and.Bio.Functions
df_scores$Diseases.and.Bio.Functions = NULL
df_scores = df_scores %>% dplyr::select(colnames(cp_scores))


overall_scores = rbind(df_scores,cp_scores)

overall_scores_long = overall_scores %>% 
  gather(key=cluster, value="z.score", -Canonical.Pathways) %>%
  separate(cluster, c("cluster","sex"), sep="_") 

head(overall_scores_long)



path_to_select=c("Glycolysis I",
                 "BAG2 Signaling Pathway",
                 "Immunogenic Cell Death Signaling Pathway",
                 "CLEAR Signaling Pathway",
                 "EIF2 Signaling",
                 "Cell movement of phagocytes",
                 "MSP-RON Signaling In Macrophages Pathway",
                 "Sirtuin Signaling Pathway",
                 "Oxytocin Signaling Pathway",
                 "ERK5 Signaling")


overall_scores_long_select = overall_scores_long %>% dplyr::filter(Canonical.Pathways %in% path_to_select)
for (cl in c('HBC','Mg','MgYSI','PAMM')){
  
  select_all_path= overall_scores_long_select %>% 
    dplyr::filter(cluster == cl ) %>%
    dplyr::filter(cluster == cl & abs(z.score)>0)
  
  des = unique(select$Canonical.Pathways)
  
  p = ggplot(data=select_all_path, aes(x=reorder(Canonical.Pathways,z.score), 
                                       y=z.score, 
                                       fill=sex,
                                       color=sex)) +
    geom_bar(stat="identity",colour="grey",size=0.25) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    geom_hline(aes(yintercept = 0)) +
    ylab("z.score") +
    xlab("")+ coord_flip() + 
    theme_minimal() + 
    ggtitle(cl)   
  p = ggplot(data=select_all_path, aes(x=reorder(Canonical.Pathways,z.score), 
                                       y=z.score, 
                                       fill=sex,
                                       color=sex)) +
    geom_bar(stat="identity",colour="grey",size=0.25) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    geom_hline(aes(yintercept = 0)) +
    ylab("z.score") +
    xlab("")+ coord_flip() + 
    theme_minimal() + 
    ggtitle(cl) 
  
  show(p)
  
  
  ggsave(p, filename = paste0(outdir,paste0(cl,"_lydiaselect_mf_barplot_combinecluster.pdf")),
         device = cairo_pdf,
         width = 7, height = 3,
         units = "in")
  
}


for (cl in c('HBC','Mg','MgYSI','PAMM')){
  
  select_all_path= overall_scores_long_select %>% 
    dplyr::filter(cluster == cl ) %>%
    mutate(Canonical.Pathways = factor(Canonical.Pathways, levels = rev(path_to_select))) 
  
  p = ggplot(data=select_all_path, aes(x=Canonical.Pathways, 
                                       y=z.score, 
                                       fill=sex,
                                       color=sex)) +
    geom_bar(stat="identity",colour="grey",size=0.25) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    geom_hline(aes(yintercept = 0)) +
    ylab("z.score") +
    xlab("")+ coord_flip() + 
    theme_minimal() + 
    ggtitle(cl)   
  
  
  ggsave(p, filename = paste0(outdir,paste0(cl,"_lydiaselect_customorder_mf_barplot_combinecluster.pdf")),
         device = cairo_pdf,
         width = 7, height = 3,
         units = "in")
  
}


# Fig 4D - Single category IPA plot ----

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
  
  
  # check max and min
  # now choose the limits
  
  heatmap_max = max(deg_select$avg_log2FC)
  heatmap_min = min(deg_select$avg_log2FC)
  
  print(paste0("heatmap max ", heatmap_max))
  print(paste0("heatmap min " , heatmap_min))
  print("using limits -1,1")
  
  
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
               cluster_columns = FALSE,
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(p_select_mat[i, j] < 0.05 & abs(lfc_select_mat[i,j])>0.2) {
                   grid.text("*", x, y)
                 }
               },
               heatmap_legend_param = list(title = "Scaled log\nfold change")) 
  return(ph)
}  

select_clusters_combined=c("Mg","MgYSI","HBC")

# make the combined cluster column the main cluster column for the plotting function

merged_seurat = read.csv(paste0(indir, "de.mfseparate.cc.obsctr_fix35.obs_vs_ctr.latent_pair_lb2_mast_logfcthresh0_maleandfemale.csv"))


deg_m = merged_seurat %>% dplyr::filter(sex == "M")
deg_f = merged_seurat %>% dplyr::filter(sex == "F")
p_thresh = 0.05
lfc_thresh = 0.1

# check max and min


p = plot_ipa_gene_heatmap("Neuroinflammation Signaling Pathway", deg_m, p_thresh, lfc_thresh, select_clusters=select_clusters_combined)
pdf(paste0(outdir, "neuro_inflam_m.pdf"), height=2.5, width=7.5)
print(p)
dev.off()

p = plot_ipa_gene_heatmap("Neuroinflammation Signaling Pathway", deg_f, p_thresh, lfc_thresh, select_clusters=select_clusters_combined)
pdf(paste0(outdir, "neuro_inflam_f.pdf"), height=2.5, width=7.5)
print(p)
dev.off()
