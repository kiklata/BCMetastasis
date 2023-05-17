library(Seurat)
library(NMF)
library(parallel)

ncores = 16
range = 2:25

# NMF ---------------------------------------------------------------------
atac = SCTransform(atac,vars.to.regress = 'percent.mt',method = "glmGamPoi")
data = as.matrix(GetAssayData(atac, assay = 'SCT', slot = 'scale.data'))
#data = data[VariableFeatures(atac),]
data[data < 0] = 0
data = data[apply(data, 1, var) > 0, ]
print(dim(data))

res.list = mclapply(range, function(r){
  nmf(data, rank = r, seed = 'ica', method = 'nsNMF')
}, mc.cores = ncores)

names(res.list) = range
saveRDS(res.list, file = paste0('.res.list.RData'))



# NMFmodule ---------------------------------------------------------------


NMFToModules = function(
    res,
    gmin = 5
){
  
  scores = basis(res)
  coefs = coefficients(res)
  
  # Remove if fewer than gmin genes
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  l = sapply(modules, length)
  keep = (l >= gmin)
  scores = scores[, keep]
  coefs = coefs[keep, ]
  
  # Find modules
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  
  names(modules) = sapply(modules, '[', 1)
  names(modules) = paste('m', names(modules), sep = '_')
  names(modules) = gsub('-','_',names(modules))
  
  return(modules)
}


# pbmcNMF -----------------------------------------------------------------


library(NMF) 
library(Seurat) 
pbmc <- readRDS("pbmc.rds") 
pbmc <- CreateSeuratObject(pbmc@assays$RNA@counts, meta.data = pbmc@meta.data)

pbmc <- NormalizeData(pbmc) %>% FindVariableFeatures() %>%
  ScaleData(do.center = F) 
vm <- pbmc@assays$RNA@scale.data 
saveRDS(vm, file = "pbmc_vm.rds") 
res <- nmf(vm, 12, method = "snmf/r") #很慢 
save(res, file = "pbmc_nmf_res.rda")

# 每个因子提取30个 
fs <- extractFeatures(res, 30L) 
fs <- lapply(fs, function(x) rownames(res)[x]) 
fs <- do.call("rbind", fs) 
rownames(fs) <- paste0("cluster", 1:12)
write.csv(t(fs), "pb,c_NMF_TopGenes.csv") 
DT::datatable(t(fs))

### 选择用于后续分析的因子 
s.f <- 1:12 # 因子 1 主要是线粒体和核糖体 
## 降维 
cell1 <- colnames(pbmc) 
cell2 <- colnames(coef(res)) 
cells <- intersect(cell1, cell2) 
pbmc <- pbmc[,cells] 
pbmc <- RunPCA(pbmc, verbose = F) 
pbmc@reductions$nmf <- pbmc@reductions$pca 
pbmc@reductions$nmf@cell.embeddings <- t(coef(res)[,cells]) 
pbmc@reductions$nmf@feature.loadings <- basis(res) 
pbmc <- RunUMAP(pbmc, reduction='nmf', dims=s.f) 

## 基于NMF降维矩阵的聚类 
pbmc <- FindNeighbors(pbmc, reduction='nmf', dims=s.f) %>% FindClusters() 

## 基于因子最大载荷分类 
pbmc$cluster <- apply(NMF::coefficients(res)[s.f,], 2, which.max)

p1 <- DimPlot(pbmc, label = T) + ggtitle("Clustered by Louvain") 
p2 <- DimPlot(pbmc, group.by = "cluster", label = T) + ggtitle("Clustered by max loading") 
pc <- p1|p2 
ggsave("pbmc_NMF_Cluster.pdf", pc, width = 10, height = 5)