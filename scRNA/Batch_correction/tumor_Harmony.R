options(warn = -1)

library(Seurat)
library(dplyr)
library(harmony)
library(glmGamPoi)
library(ggplot2)
library(patchwork)
source('~/scRNA_BC_metastases/Analysis/Function_definition/Batch_dim.R')

setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor')

# raw.batch ---------------------------------------------------------------

tumor = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds')

tumor = tumor %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(features = rownames(tumor)) %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

p1 = dim_batch_plot(tumor)
p3 = DimPlot(tumor,group.by = 'sample',label = T)+NoLegend()
ggsave('batch_dim.pdf',plot = p1,dpi = 72,width = 16,height = 6)
ggsave('batch_sample.pdf',plot = p3,dpi = 72, width = 16, height = 12)

saveRDS(tumor,file = 'tumor.batch.rds')

# harmony.batch -----------------------------------------------------------

tumor = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds')

tumor.list = SplitObject(tumor, split.by = 'sample')
tumor.list = lapply(X = tumor.list, FUN = SCTransform, method = "glmGamPoi",vars.to.regress = "percent.mt")
var_features_ind_sct = SelectIntegrationFeatures(tumor.list, nfeatures = 3000)

tumor = merge(x = tumor.list[[1]], y = tumor.list[2:length(tumor.list)])

VariableFeatures(tumor) = var_features_ind_sct

tumor = RunPCA(tumor)

tumor = RunHarmony(tumor,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                 kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

tumor = tumor %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

p2 = dim_batch_plot(tumor)
p4 = DimPlot(tumor,group.by = 'sample',label = T)+NoLegend()
ggsave('harmony_dim.pdf',plot = p2,dpi = 72,width = 16,height = 6)
ggsave('harmony_sample.pdf',plot = p4,dpi = 72, width = 16, height = 12)

saveRDS(tumor,file = 'tumor.harmony.rds')
