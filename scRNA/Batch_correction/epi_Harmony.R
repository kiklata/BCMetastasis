## need modified for sample with <200 epi cells

library(Seurat)
library(dplyr)
library(harmony)
library(glmGamPoi)
library(ggplot2)
library(patchwork)
source('~/scRNA_BC_metastases/Analysis/Function_definition/function.R')

setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/epi')

# raw.batch ---------------------------------------------------------------

epi = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/epi/epi.rds')

epi = epi %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(vars.to.regress = 'percent.mt') %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

p1 = dim_batch_plot(epi)
p3 = DimPlot(epi,group.by = 'sample',label = T)+NoLegend()
ggsave('batch_dim.pdf',plot = p1,dpi = 72,width = 16,height = 6)
ggsave('batch_sample.pdf',plot = p3,dpi = 72, width = 16, height = 12)

saveRDS(epi,file = 'epi.batch.rds')

# harmony.batch -----------------------------------------------------------

epi = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/epi/epi.rds')

epi.list = SplitObject(epi, split.by = 'sample')
epi.list = lapply(X = epi.list, FUN = SCTransform, method = "glmGamPoi",vars.to.regress = "percent.mt")
var_features_ind_sct = SelectIntegrationFeatures(epi.list, nfeatures = 3000)

epi = merge(x = epi.list[[1]], y = epi.list[2:length(epi.list)])

VariableFeatures(epi) = var_features_ind_sct

epi = RunPCA(epi)

epi = RunHarmony(epi,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                 kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

epi = epi %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

p2 = dim_batch_plot(epi)
p4 = DimPlot(epi,group.by = 'sample',label = T)+NoLegend()
ggsave('harmony_dim.pdf',plot = p2,dpi = 72,width = 16,height = 6)
ggsave('harmony_sample.pdf',plot = p4,dpi = 72, width = 16, height = 12)

saveRDS(epi,file = 'epi.harmony.rds')
