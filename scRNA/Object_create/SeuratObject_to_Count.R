library(Seurat)
library(dplyr)

setwd('~')

# T cell ------------------------------------------------------------------


tcell = readRDS('~/Tcell.rds')

sample.names = names(table(tcell$sample))

for (i in 1:length(sample.names)) {
  
  #prefix = 'tumor'
  seu = subset(tcell,sample ==sample.names[i])
  seu = CreateSeuratObject(counts = seu@assays$RNA@counts,min.cells = 3)
  count = seu[['RNA']]@counts
  write.csv(x = as.data.frame(count),file = paste0("~/tcell/",sample.names[i],'.csv'),row.names = T)
  
}


# B cell ------------------------------------------------------------------


bcell = readRDS('~/Bcell.rds')

sample.names = names(table(bcell$sample))

for (i in 1:length(sample.names)) {
  
  #prefix = 'tumor'
  seu = subset(bcell,sample ==sample.names[i])
  seu = CreateSeuratObject(counts = seu@assays$RNA@counts,min.cells = 3)
  count = seu[['RNA']]@counts
  write.csv(x = as.data.frame(count),file = paste0("~/bcell/",sample.names[i],'.csv'),row.names = T)
  
}


# mono --------------------------------------------------------------------

mono = readRDS('~/Mono.Macro.rds')

sample.names = names(table(mono$sample))

for (i in 1:length(sample.names)) {
  
  #prefix = 'tumor'
  seu = subset(mono,sample ==sample.names[i])
  seu = CreateSeuratObject(counts = seu@assays$RNA@counts,min.cells = 3)
  count = seu[['RNA']]@counts
  write.csv(x = as.data.frame(count),file = paste0("~/mono/",sample.names[i],'.csv'),row.names = T)
  
}


# dc ----------------------------------------------------------------------

dc = readRDS('~/DC.rds')

sample.names = names(table(dc$sample))

for (i in 1:length(sample.names)) {
  
  #prefix = 'tumor'
  seu = subset(dc,sample ==sample.names[i])
  
  if(dim(seu)[2] > 3){
  seu = CreateSeuratObject(counts = seu@assays$RNA@counts,min.cells = 3)
  }else{
    seu = CreateSeuratObject(counts = seu@assays$RNA@counts)
  }
  
  count = as.matrix(seu[['RNA']]@counts)
  write.csv(x = as.data.frame(count),file = paste0("~/dc/",sample.names[i],'.csv'),row.names = T)
  
}


# ILC ---------------------------------------------------------------------

ilc = readRDS('~/ILC.rds')

sample.names = names(table(ilc$sample))

for (i in 1:length(sample.names)) {
  
  #prefix = 'tumor'
  seu = subset(ilc,sample ==sample.names[i])
  if(dim(seu)[2] > 3){
    seu = CreateSeuratObject(counts = seu@assays$RNA@counts,min.cells = 3)
  }else{
    seu = CreateSeuratObject(counts = seu@assays$RNA@counts)
  }
  
  count = as.matrix(seu[['RNA']]@counts)
  write.csv(x = as.data.frame(count),file = paste0("~/ilc/",sample.names[i],'.csv'),row.names = T)
  
}
