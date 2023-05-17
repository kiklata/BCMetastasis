library(Seurat)
library(dplyr)
library(harmony)
library(glmGamPoi)
library(ggplot2)
library(patchwork)
source('~/scRNA_BC_metastases/script/function.R')

# maybe two round use metacell first
# raw.batch ---------------------------------------------------------------

all = readRDS('~/scRNA_BC_metastases/merge/object/all.rds')

all = all %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(vars.to.regress = 'percent.mt') %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

p1 = dim_batch_plot(all)
p3 = DimPlot(all,group.by = 'sample',label = T)+NoLegend()
ggsave('batch_dim.pdf',plot = p1,dpi = 72,width = 16,height = 6)
ggsave('batch_sample.pdf',plot = p3,dpi = 72, width = 16, height = 12)


# harmony.batch -----------------------------------------------------------

all = readRDS('~/scRNA_BC_metastases/merge/object/all.rds')

all.list = SplitObject(all, split.by = 'sample')
all.list = lapply(X = all.list, FUN = SCTransform, method = "glmGamPoi",vars.to.regress = "percent.mt")
var_features_ind_sct = SelectIntegrationFeatures(all.list, nfeatures = 3000)

all = merge(x = all.list[[1]], y = all.list[2:length(all.list)])

VariableFeatures(all) = var_features_ind_sct

all = RunPCA(all)

all = RunHarmony(all,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                        kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

all = all %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

p2 = dim_batch_plot(all)
p4 = DimPlot(all,group.by = 'sample',label = T)+NoLegend()
ggsave('harmony_dim.pdf',plot = p2,dpi = 72,width = 16,height = 6)
ggsave('harmony_sample.pdf',plot = p4,dpi = 72, width = 16, height = 12)

saveRDS(all,file = '~/scRNA_BC_metastases/merge/object/all.harmony.rds')
