# singler customed ref create

library(Seurat)
library(SummarizedExperiment)
library(scater)
library(dplyr)



# primary -----------------------------------------------------------------

primary = merge(CID4290A,c(CID4471,CID4495))
primary = subset(primary,cells = primary.anno$cell.id)

primary.meta = primary@meta.data
primary.meta$index = rownames(primary.meta)
primary.pdata = primary.meta[,c('index','celltype.compartment')]
primary.pdata$index = NULL

primary.count = primary@assays$RNA@counts

ref.primary = SummarizedExperiment(assays=list(counts=primary.count),colData = primary.pdata) 
ref.primary = logNormCounts(ref.primary)
saveRDS(ref.primary,file = 'ref.primary.rds')

# brain -------------------------------------------------------------------
brain = merge(pid1,c(pid2,pid3))
brain = subset(brain,cells = cell.anno$cell.id)
brain.meta = brain@meta.data
brain.meta$index = rownames(brain.meta)
brain.pdata = brain.meta[,c('index','celltype.compartment')]
brain.pdata$index = NULL

brain.count = brain[['RNA']]@counts

ref.brain = SummarizedExperiment(assays=list(counts=brain.count),colData = brain.pdata) 
ref.brain = logNormCounts(ref.brain)
saveRDS(ref.brain,file = 'ref.brain.rds')

# bone --------------------------------------------------------------------

bone = merge(bom1,bom2,add.cell.ids = c('BoM1','BoM2'))
bone = subset(bone,cells = bone.anno$cell.id)

bone.meta = bone@meta.data
bone.meta$index = rownames(bone.meta)
bone.pdata = bone.meta[,c('index','celltype.compartment')]
bone.pdata$index = NULL

bone.count = bone[['RNA']]@counts

ref.bone = SummarizedExperiment(assays=list(counts=bone.count),colData = bone.pdata) 
ref.bone = logNormCounts(ref.bone)
saveRDS(ref.bone,file = 'ref.bone.rds')
