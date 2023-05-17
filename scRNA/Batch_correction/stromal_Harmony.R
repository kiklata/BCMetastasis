## need modified for sample with <200 stromal cells


library(Seurat)
library(dplyr)
library(harmony)
library(glmGamPoi)
library(ggplot2)
library(patchwork)
source('~/scRNA_BC_metastases/Analysis/Function_definition/function.R')

setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/stromal')

# raw.batch ---------------------------------------------------------------

stromal = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/stromal/stromal.rds')

stromal = stromal %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(vars.to.regress = 'percent.mt') %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

p1 = dim_batch_plot(stromal)
p3 = DimPlot(stromal,group.by = 'sample',label = T)+NoLegend()
ggsave('batch_dim.pdf',plot = p1,dpi = 72,width = 16,height = 6)
ggsave('batch_sample.pdf',plot = p3,dpi = 72, width = 16, height = 12)

saveRDS(stromal,file = 'stromal.batch.rds')

# harmony.batch -----------------------------------------------------------

stromal = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/stromal/stromal.rds')

stromal.list = SplitObject(stromal, split.by = 'sample')
stromal.list = lapply(X = stromal.list, FUN = SCTransform, method = "glmGamPoi",vars.to.regress = "percent.mt")
var_features_ind_sct = SelectIntegrationFeatures(stromal.list, nfeatures = 3000)

stromal = merge(x = stromal.list[[1]], y = stromal.list[2:length(stromal.list)])

VariableFeatures(stromal) = var_features_ind_sct

stromal = RunPCA(stromal)

stromal = RunHarmony(stromal,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                 kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

stromal = stromal %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

p2 = dim_batch_plot(stromal)
p4 = DimPlot(stromal,group.by = 'sample',label = T)+NoLegend()
ggsave('harmony_dim.pdf',plot = p2,dpi = 72,width = 16,height = 6)
ggsave('harmony_sample.pdf',plot = p4,dpi = 72, width = 16, height = 12)

saveRDS(stromal,file = 'stromal.harmony.rds')
