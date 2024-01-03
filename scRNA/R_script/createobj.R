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


for(i in c('cellranger','decontx')){
obj_create(
 type = i,
 datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM',
 savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM',
 SampleID = 'BoM1',
 Site = 'Bone',
 Kit = 'sc3',
 Age = '72',
 meta_ER = '80',
 meta_PR = '5',
 meta_HER = '2',
 meta_KI67 = '30')
  
obj_create(
 type = i,
 datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM',
 savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM',
 SampleID = 'BoM2',
 Site = 'Bone',
 Kit = 'sc3',
 Age = '48'
 )

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM',
  SampleID = 'BoM3',
  Site = 'Bone',
  Kit = 'sc5'
  )

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM',
  SampleID = 'BrM1',
  Site = 'Brain',
  Kit = 'sn3',
  Age = '64',
  meta_ER = '0',
  meta_PR = '0',
  meta_HER = '3',
  meta_KI67 = '25'
)

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM',
  SampleID = 'BrM2',
  Site = 'Brain',
  Kit = 'sn3',
  Age = '49',
  meta_ER = '10',
  meta_PR = '0',
  meta_HER = '3',
  meta_KI67 = '40')

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM',
  SampleID = 'BrM3',
  Site = 'Brain',
  Kit = 'sn3',
  Age = '44',
  meta_ER = '80',
  meta_PR = '80',
  meta_HER = '0',
  meta_KI67 = '70')

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM',
  SampleID = 'BrM4',
  Site = 'Brain',
  Kit = 'sn3',
  Age = '66',
  meta_ER = '0',
  meta_PR = '0',
  meta_HER = '1',
  meta_KI67 = '70')

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LungM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LungM',
  SampleID = 'LungM1',
  Site = 'Lung',
  Kit = 'sn3',
  Age = '50',
  meta_ER = '0',
  meta_PR = '0',
  meta_HER = '2',
  meta_KI67 = '90')

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  SampleID = 'LiverM1',
  Site = 'Liver',
  Kit = 'sc5',
  Age = '51',
  meta_ER = '40',
  meta_PR = '0',
  meta_HER = '0',
  meta_KI67 = '30')

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  SampleID = 'LiverM2',
  Site = 'Liver',
  Kit = 'sc5',
  Age = '46',
  meta_ER = '70',
  meta_PR = '0',
  meta_HER = '1',
  meta_KI67 = '30')

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  SampleID = 'LiverM3',
  Site = 'Liver',
  Kit = 'sc5',
  Age = '63',
  meta_ER = '0',
  meta_PR = '0',
  meta_HER = '0',
  meta_KI67 = '60')

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  SampleID = 'LiverM4',
  Site = 'Liver',
  Kit = 'sc5',
  Age = '44',
  meta_ER = '0',
  meta_PR = '0',
  meta_HER = '3',
  meta_KI67 = '40')

obj_create(
  type = i,
  datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM',
  SampleID = 'LiverM5',
  Site = 'Liver',
  Kit = 'sc5',
  Age = '56',
  meta_ER = '0',
  meta_PR = '0',
  meta_HER = '2',
  meta_KI67 = '30')

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
