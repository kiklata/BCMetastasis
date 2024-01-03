library(Seurat)

setwd('~/PaperCD8/data/Tex')

CD8.Tex.harmony <- readRDS("~/CD8.Tex.harmony.rds")

before = subset(CD8.Tex.harmony,sample.timepoint == 'Before')
after = subset(CD8.Tex.harmony,sample.timepoint == 'After')

sample = names(table(before$sample.ID))

for(i in 1:length(sample)){
  
  seu = subset(before,sample.ID == sample[i])
  seu <- NormalizeData(seu, normalization.method = "RC",scale.factor = 1e6)
  
  data = as.data.frame(seu@assays$RNA@data)
  write.table(data,file = paste0('compass/before/',sample[i],'.tsv'),
              row.names = T,col.names = T,sep = '\t')
}


sample = names(table(after$sample.ID))

for(i in 1:length(sample)){
  
  seu = subset(after,sample.ID == sample[i])
  seu <- NormalizeData(seu, normalization.method = "RC",scale.factor = 1e6)
  
  data = as.data.frame(seu@assays$RNA@data)
  write.table(data,file = paste0('compass/after/',sample[i],'.tsv'),
              row.names = T,col.names = T,sep = '\t')
}

