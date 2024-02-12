library(svglite)
library(circlize)
library(ComplexHeatmap)
library(Seurat)
library(CellChat)

#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

# all cell types ----
so = readRDS("")
cellchat <- createCellChat(object = so, group.by = "rna_cell_type_sub")
cellchat@DB <- CellChatDB.human
cellchat = subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat, "" )