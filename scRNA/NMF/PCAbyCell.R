# NMF pca by cell

brain.subset = brain[,sample(dim(brain)[2],size = 1000)]
lymph.subset = lymph[,sample(dim(lymph)[2],size = 1000)]
primary.subset = primary[,sample(dim(primary)[2],size = 1000)]

meta.all = rbind(bone@meta.data,brain.subset@meta.data,lymph.subset@meta.data,primary.subset@meta.data)

program.count = t(meta.all[,c(27:36)])

program.meta = meta.all[,c(13,14)]

program.seu = CreateSeuratObject(program.count,meta.data = program.meta)

program.seu = program.seu %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()

DimPlot(program.seu,group.by = 'site')
