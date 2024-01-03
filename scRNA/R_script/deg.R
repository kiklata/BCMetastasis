seu <- readRDS("~/Project/MultiOmics/data/skin/res/cellbender_anno_count.rds")

seu = seu %>% NormalizeData() %>% ScaleData()

Idents(seu) = seu$cl_minor
markers = FindAllMarkers(seu, test.use = 'wilcox')

# subset
marker_list = list()
top_m_list = list()

for(i in c('T cell','Myeloid','Endothelial','Keratinocyte','Fibroblast')){
  
  seu_sub = subset(seu, cl_major == i)
  seu_sub = seu_sub %>% NormalizeData() %>% ScaleData()
  Idents(seu_sub) = seu_sub$cl_subset
  
  marker_list[[i]] = FindAllMarkers(seu_sub, test.use = 'wilcox')
  top_m_list[[i]] = marker_list[[i]] %>% group_by(cluster) %>% top_n(30,pct.1-pct.2)
}


saveRDS(marker_list, file = 'cellbender_cl_subset_marker.rds')
saveRDS(top_m_list, file = 'cellbender_cl_subset_topmarker.rds')


seu_sub = seu %>% subset(.,cl_major == 'Keratinocyte')
seu_sub = seu_sub %>% NormalizeData() %>% ScaleData()
Idents(seu_sub) = seu_sub$SampleType
markers = FindAllMarkers(seu_sub, test.use = 'wilcox')
