library(Seurat)
library(dplyr)

all_filter = all_filter %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() 
saveRDS(all_filter,file = 'merge/merge_filter_pca.rds')

print('pca saved')

# RPCA ---------------------------------------------------------------------

all_filter_list <- SplitObject(all_filter, split.by = 'sample')

all_filter_list <- lapply(X = all_filter_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = all_filter_list, nfeatures = 3000)
all_filter_list <- PrepSCTIntegration(object.list = all_filter_list, anchor.features = features)
all_filter_list <- lapply(X = all_filter_list, FUN = RunPCA, features = features)

print('pca')

all.anchors <- FindIntegrationAnchors(object.list = all_filter_list,normalization.method = "SCT",
                                      anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20)
print('anchor found')

all_filter_merge <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT", dims = 1:50)

DefaultAssay(all_filter_merge) <- "integrated"

all_filter_merge = all_filter_merge %>% 
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = 0.8)

print('rpca integrated')

saveRDS(all_filter_merge,file = 'merge/merge_rpca.rds')

print('rpca saved')


