library(Seurat)
library(SeuratDisk)

setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/h5ad')

tumor = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds')

sample.names = names(table(tumor$sample))

for (i in 1:length(sample.names)) {
  
  prefix = 'tumor'
  seu = subset(tumor,sample ==sample.names[i])
  SaveH5Seurat(seu, filename = paste0(sample.names[i],".h5Seurat"),assay = 'RNA')
  Convert(paste0(sample.names[i],".h5Seurat"), dest = "h5ad")
  
}
