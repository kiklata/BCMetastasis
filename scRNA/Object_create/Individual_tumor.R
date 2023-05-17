library(Seurat)

setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/rds')

tumor = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds')

sample.names = names(table(tumor$sample))

for (i in 1:length(sample.names)) {
  
  seu = subset(tumor,sample ==sample.names[i])
  saveRDS(seu,file = paste0(sample.names[i],'.rds'))
}
