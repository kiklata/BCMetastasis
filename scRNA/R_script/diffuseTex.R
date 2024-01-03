library(Seurat)
library(dplyr)
library(destiny)

CD8.Tex.harmony <- readRDS("~/PaperCD8/data/Tex/CD8.downsample.rds")

GeneExp.mtx <- GetAssayData(CD8.Tex.harmony, assay = "RNA", slot = "data") %>% as.matrix() # normalized data matrix

dm <- DiffusionMap(t(GeneExp.mtx),n_pcs = 50)

saveRDS(dm,file = 'CD8.DM.rds')