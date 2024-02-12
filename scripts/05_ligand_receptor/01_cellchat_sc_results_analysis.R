library(svglite)
library(circlize)
library(ComplexHeatmap)
library(CellChat)
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html


# cellchat analysis ------------------------------------------------------------

cellchat = readRDS("" )
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

groupSize <- as.numeric(table(cellchat@idents))
groupSize

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

p = netVisual_heatmap(cellchat, signaling = "GALECTIN", color.heatmap = "Reds") +
  xlab("") + ylab("")
print(p)

for (path in cellchat@netP$pathways){
  print(path)
  p = netVisual_heatmap(cellchat, signaling = path, color.heatmap = "Reds")
  pdf(paste0(path,'.pdf'),width = 4,height = 4)
  draw(p)
  dev.off()

}

for (path in cellchat@netP$pathways){
  print(path)
  p = netVisual_heatmap(cellchat, signaling = path, color.heatmap = "Reds")
  pdf(paste0(path,'.pdf'),width = 4,heigh = 4)
  draw(p)
  dev.off()
  
}

gg <- netAnalysis_contribution(cellchat, signaling = cellchat@netP$pathways[1])
show(gg)
ggsave(filename="",
       plot=p, width = 3, height = 2, units = 'in', dpi = 300)


# write out all  ----------------------------------------------
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "", row.names=F)

# what are the top pairs?
df.net_filter = df.net %>%
  dplyr::filter(prob > 0.1)

unique(df.net$pathway_name)
write.xlsx(df.net_filter, "", row.names=F)

df.net <- subsetCommunication(cellchat, sources.use = c("Ciliated Epithelial","Epithelial"), targets.use = c("neut_0","neut_1","neut_2"))
write.csv(df.net, "", row.names=F)

view(df.net)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(neut_epi_cellchat, 
                                         pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(neut_epi_cellchat, pattern = "all")
ht1 + ht2




netAnalysis_signalingRole_scatter(neut_epi_cellchat, signaling = c("CXCL", "CCL"))

netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)



netAnalysis_signalingRole_heatmap(cellchat, 
                                        signaling = cellchat@netP$pathway)
