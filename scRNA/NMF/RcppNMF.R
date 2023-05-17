library(Seurat)
library(RcppML)

tumor = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds')

sample.l = names(table(tumor$sample))
# filter sample with <100 tumor cell
sample.l.100 = sample.l[c(-1,-11,-16,-19,-44)]

rank.range = 2:25
setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/rcpp')

for (i in 1:length(sample.l.100)) {
 
seu = subset(tumor,sample==sample.l.100[i])
seu = SCTransform(seu)
data = seu@assays$SCT@scale.data
data[data < 0] = 0
data = data[apply(data, 1, var) > 0, ]

for (k in 1:length(rank.range)) {
  
  rank = rank.range[k]
  
  res = RcppML::nmf(data,rank)
  
  w <- res$w
  rownames(w)<- rownames(data)
  colnames(w)<- paste0("component", 1:rank)
  
  w = as.data.frame(w)
  new.w = data.frame(cluster1 = rownames(w)[order(w[,1],decreasing = T)][1:30])
  dir.create(paste0('rank',k+1))
  setwd(paste0('rank',k+1))
  for (j in 2:ncol(w)) {
    new.w[,j] = data.frame(rownames(w)[order(w[,j],decreasing = T)][1:30])
  }
  colnames(new.w) = colnames(w)
  saveRDS(new.w,file = paste0(sample.l.100[i],'.rds'))
  setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/rcpp')
}
}