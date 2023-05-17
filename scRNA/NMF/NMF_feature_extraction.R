library(NMF) 
library(Seurat) 


nmf.res = readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/100cell/cell.pid2.res.list.rds")

# NMFgene -------------------------------------------------------------------
NMFtoGenes = function(nmf.res,genenum = 30){
  
  fs.list = list()
  
  for (i in 1:length(nmf.res)) {
    res = nmf.res[[i]]
    fs = extractFeatures(res, genenum) 
    fs = lapply(fs, function(x) rownames(res)[x]) 
    fs = do.call("rbind", fs) 
    rownames(fs) = paste0("cluster", 1:(i+1))
    fs = as.data.frame(t(fs))
    fs.list[[i]] = fs
    names(fs.list)[i] = paste0('rank',i+1)
  }
  
  return(fs.list)
  
}

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


nmf.gene = NMFtoGenes(nmf.res = nmf.res)

# modual anno GO 

library(clusterProfiler)
library(org.Hs.eg.db)
#library(createKEGGdb)
#create_kegg_db('hsa')
#library(KEGG.db)

gmtfile = read.gmt("~/bioinfo/gmt/msigdb_v2022.1.Hs_GMTs/h.all.v2022.1.Hs.symbols.gmt")

# need convert gene symbol
go.res.list = data.frame()

for (i in 1:length(nmf.gene)) {
  
  for (k in 1:ncol(nmf.gene[[i]])) {

  genelist = nmf.gene[[i]][,k]
  
  genelist.ENTREZID = bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,2]
  
  ego = enrichGO(genelist.ENTREZID,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",
                 ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",
                 #universe = ,qvalueCutoff = 0.2,
                 minGSSize = 10,maxGSSize = 500,readable = FALSE,pool = FALSE)
  result = ego@result
  result = result[order(result$GeneRatio,decreasing = T),]
  go.res.list[i,k] = result$Description[1]
  rownames(go.res.list)[i] = paste0('rank',i+1)
  colnames(go.res.list)[k] = paste0('cluster',k)
  }
}
#ekg = enrichKEGG(genelist.ENTREZID,organism = "hsa",keyType = "kegg",
#                 pvalueCutoff = 0.05,pAdjustMethod = "BH",
                 #universe,qvalueCutoff = 0.2,
#                 minGSSize = 10,maxGSSize = 500,use_internal_data = F)





























### 选择用于后续分析的因子 
s.f <- 1:5 # 因子 1 主要是线粒体和核糖体 
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