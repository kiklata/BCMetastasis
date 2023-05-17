library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SCENIC)

setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/h5ad/filter_NG')

tumor = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds')

sample.names = names(table(tumor$sample))

for (i in 1:length(sample.names)) {
  
  seu = subset(tumor,sample ==sample.names[i])
  # pySCENIC thresholds
  
  count = seu@assays$RNA@counts
  nCountsPerGene = apply(count,MARGIN = 1,sum)
  count[count>0] = 1
  nCellsPerGene = apply(count,MARGIN = 1,sum)
  
  nCells= ncol(seu)
  minCountsPerGene = round(.05*nCells,2) # 0.05 counts per cell
  minSamples=round(.05*nCells,0) # 5% of cells
  
  GeneSelectCounts = nCountsPerGene[nCountsPerGene>minCountsPerGene]
  GeneSelectSamples = nCellsPerGene[nCellsPerGene>minSamples]
  GeneSelect = intersect(names(GeneSelectCounts),names(GeneSelectSamples))
  
  exp = seu@assays$RNA@counts
  exp = exp[rownames(exp) %in% GeneSelect,]
  seufilter = CreateSeuratObject(exp)
  seufilter = AddMetaData(seufilter,metadata = seu@meta.data[,c(-1,-2,-3)])
  SaveH5Seurat(seufilter, filename = paste0(sample.names[i],".h5Seurat"),assay = 'RNA')
  Convert(paste0(sample.names[i],".h5Seurat"), dest = "h5ad")
  
}


