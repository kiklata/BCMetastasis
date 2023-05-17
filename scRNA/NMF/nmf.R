#NMF_loop_100

library(Seurat)
library(NMF)
library(parallel)

range = 2:25
ncores = 8

tumor = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds')

sample.l = names(table(tumor$sample))
# filter sample with <100 tumor cell
sample.l.100 = sample.l[c(-1,-11,-16,-19,-44)]
# loop --------------------------------------------------------------------
setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/final')

for(i in 1:length(sample.l.100)){
  
seu = subset(tumor,sample == sample.l.100[i])
seu = NormalizeData(seu)
seu = ScaleData(seu, features = rownames(seu))
data = as.matrix(GetAssayData(seu, assay = 'RNA', slot = 'scale.data'))
data[data < 0] = 0
data = data[apply(data, 1, var) > 0, ]
print(dim(data))

   res.list2 = mclapply(range, function(r){
   nmf(data, rank = r, seed = 'ica', method = 'nsNMF',.options = 'v',maxIter = 5000)
   }, mc.cores = ncores)
   names(res.list2) = range
#  res.list = nmf(data, rank = range, seed = 'ica', method = 'nsNMF',.options = 'v')
  saveRDS(res.list2, file = paste0(sample.l.100[i],'.res.list.rds'))

}