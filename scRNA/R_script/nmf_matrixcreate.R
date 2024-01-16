path = '~/Project/scRNA_BC_metastases/Data/SingleCell/public/BrM/science/'

seu <- readRDS(paste0(path,'LungM_decontx_qc.rds'))
meta = read.delim(paste0(path,'LungM_decontx_celltype.csv'),sep = ',',row.names = 1)
seu = AddMetaData(seu, meta)

seu_sub = subset(seu, cl_major == 'Epi')
table(seu_sub$SampleID)

sample_list = seu_sub@meta.data[["SampleID"]] %>% unique()

for (sample in sample_list) {
  seu_sub_sample = subset(seu_sub, SampleID == sample)
  count = seu_sub_sample@assays$RNA$counts %>% as.data.frame()
  count = count[rowSums(count)>0,]
  count = count[!stringr::str_detect(rownames(count),pattern = '^MT-'),]
  count = t(count) %>% as.data.frame()
  dir.create(path = paste0(path,'/',sample),showWarnings = F,recursive = T)
  write.table(count, file = paste0(path,'/',sample,'/count.txt'),quote = F,sep = '\t',row.names = T, col.names = T)
}


seu <- readRDS("~/Project/scRNA_BC_metastases/Data/SingleCell/public/raw_sc_BC/BrM/cell/BrM3/seu.rds")

seu_sub = subset(seu, cellType == 'Epithelial')
count = seu_sub@assays$RNA$counts %>% as.data.frame()
count = count[rowSums(count)>0,]
count = count[!stringr::str_detect(rownames(count),pattern = '^MT-'),]
count = t(count) %>% as.data.frame()
write.table(count, file = '~/Project/scRNA_BC_metastases/Data/SingleCell/public/BoM/BoM2/count.txt',quote = F,sep = '\t',row.names = T, col.names = T)

seu <- readRDS("~/Project/scRNA_BC_metastases/Data/SingleCell/public/raw_sc_BC/BrM/GSE202501/cellranger_doublet.rds")
meta = read.delim("~/Project/scRNA_BC_metastases/Data/SingleCell/public/raw_sc_BC/BrM/GSE202501/cellranger_celltype.csv",sep = ',',row.names = 1)

seu = AddMetaData(seu, meta)
seu = seu[,!is.na(seu$cl_major)]
saveRDS(seu,"~/Project/scRNA_BC_metastases/Data/SingleCell/public/raw_sc_BC/BrM/GSE202501/seu.rds")

seu_sub = subset(seu, cl_major == 'Epi')
count = seu_sub@assays$RNA$counts %>% as.data.frame()
count = count[rowSums(count)>0,]
count = count[!stringr::str_detect(rownames(count),pattern = '^MT-'),]
count = t(count) %>% as.data.frame()
write.table(count, file = '~/Project/scRNA_BC_metastases/Data/SingleCell/public/raw_sc_BC/BrM/GSE202501/count.txt',quote = F,sep = '\t',row.names = T, col.names = T)
