source("~/Project/scRNA_BC_metastases/Analysis/scRNA/func/obj_create.R")

#for(i in c('cellbender','cellranger','soupx')){
# obj_create(
# type = i,
# datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM',
# savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM',
# SampleID = 'BoM1',
# Site = 'Bone',
# Kit = 'sc',
# Age = '',
# ER = '',
# PR = '',
# HER = '',
# KI67 = ''
# )}

for(i in c('cellranger','decontx','soupx','cellbender')){
  
obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  SampleID = 'LiverM6',
  Site = 'Liver',
  Kit = 'sc5',
  Age = '55',
  meta_ER = '90',
  meta_PR = '90',
  meta_HER = '0',
  meta_KI67 = '8')

}

# test for decontamination effect -----------------

#seu = subset(seu, nCount_RNA <20000 & nCount_RNA >500 & nFeature_RNA<6000 & scDblFinder.class == 'singlet')
#seu = seu %>% NormalizeData(.) %>% FindVariableFeatures(.) %>% ScaleData(.) %>% RunPCA() %>%  RunUMAP(.,dims = 1:30)

#library(viridis)
#p4 = FeaturePlot(seu,c('KRT14','KRT10','KRT1'), ncol = 3,pt.size = 0.1)&scale_color_viridis_c()
#p0/p1/p2
