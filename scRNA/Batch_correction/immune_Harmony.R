library(Seurat)
library(dplyr)
library(harmony)
library(glmGamPoi)
library(ggplot2)
library(patchwork)
source('~/scRNA_BC_metastases/Analysis/scRNA/Function_definition/Batch_dim.R')

setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/immune')

# raw.batch ---------------------------------------------------------------

immune = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/immune.rds')

immune = immune %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(vars.to.regress = 'percent.mt') %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

p1 = dim_batch_plot(immune)
p3 = DimPlot(immune,group.by = 'sample',label = T)+NoLegend()
ggsave('batch_dim.pdf',plot = p1,dpi = 72,width = 16,height = 6)
ggsave('batch_sample.pdf',plot = p3,dpi = 72, width = 16, height = 12)

saveRDS(immune,file = 'immune.batch.rds')

# harmony.batch -----------------------------------------------------------

immune = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/immune.rds')

immune.list = SplitObject(immune, split.by = 'sample')
immune.list = lapply(X = immune.list, FUN = SCTransform, method = "glmGamPoi",vars.to.regress = "percent.mt")
var_features_ind_sct = SelectIntegrationFeatures(immune.list, nfeatures = 3000)

immune = merge(x = immune.list[[1]], y = immune.list[2:length(immune.list)])

VariableFeatures(immune) = var_features_ind_sct

immune = RunPCA(immune)

immune = RunHarmony(immune,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                 kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

immune = immune %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

p2 = dim_batch_plot(immune)
p4 = DimPlot(immune,group.by = 'sample',label = T)+NoLegend()
ggsave('harmony_dim.pdf',plot = p2,dpi = 72,width = 16,height = 6)
ggsave('harmony_sample.pdf',plot = p4,dpi = 72, width = 16, height = 12)

saveRDS(immune,file = 'immune.harmony.rds')


#tcell ---------------------

tcell.anno <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/tcell.anno.rds")
tcell <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/tcell.rds")
tcell = AddMetaData(tcell,tcell.anno)

tcell = tcell %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>% 
  RunHarmony(group.by.vars = 'sample',theta = 2, kmeans_init_nstart=20, kmeans_init_iter_max=5000)

tcell = tcell %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

saveRDS(tcell,"~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/tcell.harmony.rds")

p1 = DimPlot(tcell,group.by = 'celltypist.celltype.subset')
p2 = DimPlot(tcell,group.by = 'scibet.celltype.marker')
p = p1+p2
ggsave('harmony_sample.pdf',plot = p,dpi = 72, width = 16, height = 12)

#CD4--------------------
CD4 = subset(tcell,scibet.celltype.minor == 'CD4 T cell')

CD4 = CD4 %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>% 
  RunHarmony(group.by.vars = 'sample',theta = 2, kmeans_init_nstart=20, kmeans_init_iter_max=5000)

CD4 = CD4 %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

p2 = dim_batch_plot(CD4)
p4 = DimPlot(CD4,group.by = 'sample',label = T)+NoLegend()

ggsave('harmony_dim.pdf',plot = p2,dpi = 72,width = 16,height = 6)
ggsave('harmony_sample.pdf',plot = p4,dpi = 72, width = 16, height = 12)

saveRDS(CD4,file = 'cd4tcell.harmony.rds')

#CD8--------------------
CD8 = subset(tcell,scibet.celltype.minor == 'CD8 T cell')

CD8 = CD8 %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>% 
  RunHarmony(group.by.vars = 'sample',theta = 2, kmeans_init_nstart=20, kmeans_init_iter_max=5000)

CD8 = CD8 %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

p2 = dim_batch_plot(CD8)
p4 = DimPlot(CD8,group.by = 'sample',label = T)+NoLegend()

ggsave('harmony_dim.pdf',plot = p2,dpi = 72,width = 16,height = 6)
ggsave('harmony_sample.pdf',plot = p4,dpi = 72, width = 16, height = 12)

saveRDS(CD8,file = 'cd8tcell.harmony.rds')

