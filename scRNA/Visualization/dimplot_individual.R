library(Seurat)
library(dplyr)
setwd('plot')
#bone-------
bom1 = readRDS("~/scRNA_BC_metastases/meta_Bone/BoM1/RNA/filter.rds")
bom2 = readRDS("~/scRNA_BC_metastases/meta_Bone/BoM2/RNA/filter.rds")

#brain------
#cell
cell.pid1 = readRDS('~/scRNA_BC_metastases/meta_Brain/cell/PID1/filter.rds')
cell.pid2 = readRDS('~/scRNA_BC_metastases/meta_Brain/cell/PID2/filter.rds')
cell.pid3 = readRDS('~/scRNA_BC_metastases/meta_Brain/cell/PID3/filter.rds')
#GSE143423
GSE143423 = readRDS('~/scRNA_BC_metastases/meta_Brain/GSE143423/filter.rds')
#GSE202501
GSE202501 = readRDS('~/scRNA_BC_metastases/meta_Brain/GSE202501/filter.rds')
#science
science.GSM4555888 = readRDS('~/scRNA_BC_metastases/meta_Brain/science/GSM4555888/filter.rds')
science.GSM4555889 = readRDS('~/scRNA_BC_metastases/meta_Brain/science/GSM4555889/filter.rds')
science.GSM4555891 = readRDS('~/scRNA_BC_metastases/meta_Brain/science/GSM4555891/filter.rds')

#lymph----------
#atac
atac = readRDS('~/scRNA_BC_metastases/meta_Lymph/atac/filter.rds')
#EMBOJ
EMBOJ.GSM4909308 = readRDS('~/scRNA_BC_metastases/meta_Lymph/EMBOJ/GSM4909308/filter.rds')
EMBOJ.GSM4909310 = readRDS('~/scRNA_BC_metastases/meta_Lymph/EMBOJ/GSM4909310/filter.rds')
EMBOJ.GSM4909312 = readRDS('~/scRNA_BC_metastases/meta_Lymph/EMBOJ/GSM4909312/filter.rds')
EMBOJ.GSM4909314 = readRDS('~/scRNA_BC_metastases/meta_Lymph/EMBOJ/GSM4909314/filter.rds')
EMBOJ.GSM4909316 = readRDS('~/scRNA_BC_metastases/meta_Lymph/EMBOJ/GSM4909316/filter.rds')
EMBOJ.GSM4909318 = readRDS('~/scRNA_BC_metastases/meta_Lymph/EMBOJ/GSM4909318/filter.rds')
EMBOJ.GSM4909321 = readRDS('~/scRNA_BC_metastases/meta_Lymph/EMBOJ/GSM4909321/filter.rds')
#oncogenesis
oncogenesis.B2 = readRDS('~/scRNA_BC_metastases/meta_Lymph/oncogenesis/B2/filter.rds')
oncogenesis.C2 = readRDS('~/scRNA_BC_metastases/meta_Lymph/oncogenesis/C2/filter.rds')
oncogenesis.D2 = readRDS('~/scRNA_BC_metastases/meta_Lymph/oncogenesis/D2/filter.rds')
oncogenesis.D3 = readRDS('~/scRNA_BC_metastases/meta_Lymph/oncogenesis/D3/filter.rds')
oncogenesis.E2 = readRDS('~/scRNA_BC_metastases/meta_Lymph/oncogenesis/E2/filter.rds')


# primary -----------------------------------------------------------------
CID3586 = readRDS("~/scRNA_BC_metastases/primary/CID3586/filter.rds")
CID3838 = readRDS("~/scRNA_BC_metastases/primary/CID3838/filter.rds")
CID3921 = readRDS("~/scRNA_BC_metastases/primary/CID3921/filter.rds")
CID3941 = readRDS("~/scRNA_BC_metastases/primary/CID3941/filter.rds")
CID3946 = readRDS("~/scRNA_BC_metastases/primary/CID3946/filter.rds")
CID3948 = readRDS("~/scRNA_BC_metastases/primary/CID3948/filter.rds")
CID3963 = readRDS("~/scRNA_BC_metastases/primary/CID3963/filter.rds")
CID4040 = readRDS("~/scRNA_BC_metastases/primary/CID4040/filter.rds")
CID4066 = readRDS("~/scRNA_BC_metastases/primary/CID4066/filter.rds")
CID4067 = readRDS("~/scRNA_BC_metastases/primary/CID4067/filter.rds")
CID4290A = readRDS("~/scRNA_BC_metastases/primary/CID4290A/filter.rds")
CID4398 = readRDS("~/scRNA_BC_metastases/primary/CID4398/filter.rds")
CID4461 = readRDS("~/scRNA_BC_metastases/primary/CID4461/filter.rds")
CID4463 = readRDS("~/scRNA_BC_metastases/primary/CID4463/filter.rds")
CID4465 = readRDS("~/scRNA_BC_metastases/primary/CID4465/filter.rds")
CID4471 = readRDS("~/scRNA_BC_metastases/primary/CID4471/filter.rds")
CID4495 = readRDS("~/scRNA_BC_metastases/primary/CID4495/filter.rds")
CID4513 = readRDS("~/scRNA_BC_metastases/primary/CID4513/filter.rds")
CID4515 = readRDS("~/scRNA_BC_metastases/primary/CID4515/filter.rds")
CID4523 = readRDS("~/scRNA_BC_metastases/primary/CID4523/filter.rds")
CID4530N = readRDS("~/scRNA_BC_metastases/primary/CID4530N/filter.rds")
CID4535 = readRDS("~/scRNA_BC_metastases/primary/CID4535/filter.rds")
CID44041 = readRDS("~/scRNA_BC_metastases/primary/CID44041/filter.rds")
CID44971 = readRDS("~/scRNA_BC_metastases/primary/CID44971/filter.rds")
CID44991 = readRDS("~/scRNA_BC_metastases/primary/CID44991/filter.rds")
CID45171 = readRDS("~/scRNA_BC_metastases/primary/CID45171/filter.rds")

# define plot ----------------------------------------------------------

run.plot = function(obj.list,sample.list,savepath){

  library(ggplot2)
  
  feature_gene = c(
'PTPRC','CD3D','CD4','CD8A','CD79A','MS4A1','CD68', # immune CD45:PTPRC
'EPCAM','KRT19','KRT18', # tumor/epithelial
'PECAM1','COL1A1'# stromal CD31:PECAM1, CD10:MME
  )

  for (i in 1:length(obj.list)) {
    
    seu = obj.list[[i]]
    seu = seu %>% FindVariableFeatures() %>% NormalizeData() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>% FindClusters()
    
    p1 = FeaturePlot(seu,features = feature_gene)
    p2 = DotPlot(seu,features = feature_gene)+coord_flip()
    p3 = DimPlot(seu,reduction = 'umap',label = T)
    ggsave(paste0(sample.list[i],'_feature.pdf'),p1,width = 12,height = 6,dpi = 72)
    ggsave(paste0(sample.list[i],'_dotplot.pdf'),p2,width = 8,height = 6,dpi = 72)
    ggsave(paste0(sample.list[i],'_dimplot.pdf'),p3,width = 8,height = 6,dpi = 72)
    
    markers = FindAllMarkers(seu,test.use = 'MAST')
    top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    meta = seu@meta.data
    meta$cell.id = colnames(seu)
    saveRDS(meta,file = paste0(sample.list[i],'_meta.rds'))
    write.csv(top10,file = paste0(sample.list[i],'_marker.csv'),row.names = F)
    
  }
}


# run ---------------------------------------------------------------------

obj.list = c(atac, 
             EMBOJ.GSM4909308,EMBOJ.GSM4909310,EMBOJ.GSM4909312,EMBOJ.GSM4909314,
             EMBOJ.GSM4909316,EMBOJ.GSM4909318,EMBOJ.GSM4909321,
             oncogenesis.B2,oncogenesis.C2,oncogenesis.D2,oncogenesis.D3,oncogenesis.E2, 
             cell.pid1,cell.pid2,cell.pid3, 
             GSE143423,
             GSE202501,
             science.GSM4555888,science.GSM4555889,science.GSM4555891,
             bom1,bom2,
             CID3586,CID3838,CID3921,CID3941,CID3946,CID3948,CID3963,
             CID4040,CID4066,CID4067,CID4290A,CID4398,CID44041,CID4461,CID4463,CID4465,CID4471,CID4495,CID44971,CID44991,
             CID4513,CID4515,CID45171,CID4523,CID4530N,CID4535)

sample.list = c('atac', 
                'EMBOJ.GSM4909308','EMBOJ.GSM4909310','EMBOJ.GSM4909312','EMBOJ.GSM4909314',
                'EMBOJ.GSM4909316','EMBOJ.GSM4909318','EMBOJ.GSM4909321',
                'oncogenesis.B2','oncogenesis.C2','oncogenesis.D2','oncogenesis.D3','oncogenesis.E2', 
                'cell.pid1','cell.pid2','cell.pid3', 
                'GSE143423',
                'GSE202501',
                'science.GSM4555888','science.GSM4555889','science.GSM4555891',
                'bom1','bom2',
                'CID3586','CID3838','CID3921','CID3941','CID3946','CID3948','CID3963',
                'CID4040','CID4066','CID4067','CID4290A','CID4398','CID44041','CID4461','CID4463','CID4465','CID4471','CID4495','CID44971','CID44991',
                'CID4513','CID4515','CID45171','CID4523','CID4530N','CID4535')

run.plot(obj.list = obj.list,sample.list = sample.list)
