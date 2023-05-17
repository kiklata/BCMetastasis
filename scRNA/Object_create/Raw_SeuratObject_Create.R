library(Seurat)
library(dplyr)

setwd("")

# raw ---------------------------------------------------------------------
count = read.delim("",sep = ',',header = T,row.name = 1)

#count = Read10X(getwd())
seu = CreateSeuratObject(counts = count,min.cells = 3, min.features = 200)
seu$percent.mt = PercentageFeatureSet(seu,pattern = '^MT-')

# doubletfinder -----------------------------------------------------------

library(DoubletFinder)

seu = NormalizeData(seu) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 1.5)

## pK Identification (no ground-truth) 
sweep.res.list = paramSweep_v3(seu, PCs = 1:30, sct = F)
sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn = find.pK(sweep.stats)
pK_bcmvn = bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## Homotypic Doublet Proportion Estimate 
homotypic.prop = modelHomotypic(seu@meta.data$seurat_clusters)           
DoubletRate = ncol(seu)*8*1e-6
nExp_poi = round(DoubletRate*nrow(seu@meta.data))  ## Assuming 7.5% doublet formation rate 10000 cell - tailor for your dataset
nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies 
seu = doubletFinder_v3(seu, PCs = 1:30, pN = 0.25, 
                       pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
colnames(seu@meta.data)[8] = "DF.classifications"

seu = DietSeurat(seu)
saveRDS(seu,file = 'raw.rds')

# filter ------------------------------------------------------------------

seu = subset(seu, DF.classifications== 'Singlet'& percent.mt <20)

saveRDS(seu,file = 'filter.rds')
rm(list = ls())
