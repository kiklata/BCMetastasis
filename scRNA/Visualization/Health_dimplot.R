health.brain <- readRDS("~/scRNA_BC_metastases/Data/raw_sc_Health/rds/health.brain.rds")

health.brain = health.brain %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(features = rownames(health.brain)) %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

#DimPlot(health.brain,group.by = 'subclass_label',cols = names(table(table(health.brain$subclass_color))))

brain.p = ClassDimPlot(
  srt = health.brain, group.by = 'subclass_label',
  reduction = "UMAP", theme_use = "theme_blank"
) + ggtitle('Brain')

TS_Bone_Marrow <- readRDS("~/scRNA_BC_metastases/Data/raw_sc_Health/rds/TS_Bone_Marrow.rds")

#DimPlot(TS_Bone_Marrow,reduction = 'scvi_umap',group.by = 'free_annotation')

bone.p = ClassDimPlot(
  srt = TS_Bone_Marrow, group.by = 'free_annotation',
  reduction = "scvi_umap", theme_use = "theme_blank"
) + ggtitle('Bone marrow')

TS_Mammary <- readRDS("~/scRNA_BC_metastases/Data/raw_sc_Health/rds/TS_Mammary.rds")

breast.p = ClassDimPlot(
  srt = TS_Mammary, group.by = 'free_annotation',
  reduction = "scvi_umap", theme_use = "theme_blank"
) + ggtitle('Mammary')

TS_Lymph_Node <- readRDS("~/scRNA_BC_metastases/Data/raw_sc_Health/rds/TS_Lymph_Node.rds")

lymph.p = ClassDimPlot(
  srt = TS_Lymph_Node, group.by = 'free_annotation',
  reduction = "scvi_umap", theme_use = "theme_blank"
) + ggtitle('Lymph node')

p = (breast.p|lymph.p)/(bone.p|brain.p) + plot_annotation(title = 'Human Cell Atlas',
                                                      subtitle = 'The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans. Science. 2022',
                                                      tag_levels = 'A')
ggsave('health.pdf',p,width = 16,height = 9)
