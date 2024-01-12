library(SoupX)

paths = paste0('~/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM/',
               c('LiverM6'))

for (path in paths) {
  
setwd(path)

sc = load10X(path)

seu = CreateSeuratObject(sc$toc) %>% SCTransform(.) %>% RunPCA(., npcs = 40) %>%
  FindNeighbors(.,dims = 1:30) %>% FindClusters(., resolution = 0.8) %>% 
  RunUMAP(.,dims = 1:30)

matx = seu@meta.data

sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc,doPlot = F)
out = adjustCounts(sc,roundToInt = T)

new = CreateSeuratObject(out,min.cells = 3,min.features = 200)

saveRDS(new,file = 'soupx.rds')
}
