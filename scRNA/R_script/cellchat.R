library(CellChat)

cellranger_filter_count <- readRDS("~/Project/MultiOmics/data/skin/res/cellbender_anno_count.rds")
#anno <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/cellranger_filter_anno.tsv", row.names=1)

#seu = AddMetaData(cellranger_filter_count,anno)
seu = cellranger_filter_count %>% NormalizeData(.)

select_cellmajor = c('Fibroblast','Keratinocyte','T cell','Myeloid','Endothelial')
input_seu = seu %>% subset(., cl_major %in% select_cellmajor)

callcellchat = function(seu,type,cluster = 'cl_subset',ncore = 8){
  input_seu = subset(seu, SampleType == type)
  input_seu$cl_subset = droplevels(input_seu$cl_subset)
  data_input <- GetAssayData(input_seu,assay = "RNA", layer = "data")
  cellchat <- createCellChat(object = data_input, meta = input_seu@meta.data, group.by = cluster)
  cellchat@DB <- CellChatDB.human
  future::plan("multisession", workers = ncore) # do parallel
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')
  return(cellchat)
}

cellchat_H = callcellchat(seu = input_seu, type = 'H')
cellchat_N = callcellchat(seu = input_seu, timepoint = 'N')

#Calculate the aggregated cell-cell communication network
#cellchat <- mergeCellChat(list(cellchat_H,cellchat_N), add.names = c('ARD','Normal'))

object.list = list(ARD = cellchat_H, Normal = cellchat_N)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

netVisual_bubble(object.list[[1]])

pt_source =  c(4,7,9,10,13)
pt_target = c(14,15,23)

netVisual_bubble(object.list[[1]], sources.use = pt_source, targets.use = pt_target, remove.isolate = FALSE)+coord_flip()
#netVisual_chord_gene(object.list[[1]], sources.use = c(4,7,9,10,13), targets.use = c(14:15,23), lab.cex = 0.5,legend.pos.y = 30)
#netVisual_circle(object.list[[1]]@net$count,sources.use = c(4,7,9,10,13),targets.use = c(14:15,23), weight.scale = T)


pathways.show <- c("TGFb") 
netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "circle")

netVisual_aggregate(object.list[[1]], signaling = pathways.show,sources.use = pt_source,targets.use = pt_target, layout = "chord")
netVisual_heatmap(object.list[[1]], signaling = pathways.show, sources.use = pt_source,targets.use = pt_target,color.heatmap = "Reds")

pairLR.CXCL <- extractEnrichedLR(object.list[[1]], signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,]
netVisual_individual(object.list[[1]], signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

netAnalysis_signalingRole_network(object.list[[1]],signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_scatter(object.list[[1]])|netAnalysis_signalingRole_scatter(object.list[[2]])


groupSize <- as.numeric(table(object.list[[1]]@idents))

netVisual_circle(object.list[[1]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",idents.use = c('Macrophage-TREM2'))
netVisual_circle(object.list[[1]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",idents.use = c('Macrophage-TREM2'))

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram

pathways.show <- c("MIF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")


# Heatmap
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)


pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object


# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)


# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)


# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)


#Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 5, font.size = 4)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("NOTCH", "PTN"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("NOTCH", "PTN"))



