# first time ks.cut = 0.1
# try ks.cut = 0.15 to increase tumor cell number

library(Seurat)
library(copykat)
library(dplyr)

source('~/copyKat/copykatnofilter.R')

all = readRDS('~/scRNA_BC_metastases/merge/object/all.rds')
sample.n = names(table(all$sample))
all.anno = all@meta.data
sample.n = sample.n[-c(22,25,28,32)]
  
for (i in 1:45) {
    
    seu = subset(all,sample == sample.n[i])
    
    seu1 = subset(seu,celltype.compartment == 'epi/cancer')
    
    if (ncol(seu1)<500) {
      
      anno = filter(all.anno,sample == sample.n[i])
      norm.cell = rownames(filter(anno,celltype.compartment !='epi/cancer'))
      
      if(length(norm.cell)>1000){
        norm.cell = norm.cell[sample(length(norm.cell),1000)]
      }
      setwd('~/copyKat')
      
      norm = subset(seu, cells =  norm.cell)
      count = cbind(seu1[['RNA']]@counts,norm[['RNA']]@counts)
    
      dir.create(sample.n[i])
      setwd(sample.n[i])
      
      copykat.test = copykat_no_filter(rawmat=count, ngene.chr=5, sam.name=sample.n[i],genome = 'hg20', id.type = "S", 
                                       norm.cell.names = norm.cell, n.cores = 16,KS.cut = 0.1)
    }else if(ncol(seu1)>500){
      anno = filter(all.anno,sample == sample.n[i])
      norm.cell = rownames(filter(anno,celltype.compartment !='epi/cancer'))
      
      if(length(norm.cell)>2000){
        norm.cell = norm.cell[sample(length(norm.cell),2000)]
      }
      setwd('~/copyKat')
      
      norm = subset(seu, cells =  norm.cell)
      count = cbind(seu1[['RNA']]@counts,norm[['RNA']]@counts)
      
      dir.create(sample.n[i])
      setwd(sample.n[i])
      
      copykat.test = copykat_no_filter(rawmat=count, ngene.chr=5, sam.name=sample.n[i],genome = 'hg20',  id.type = "S",
                                       norm.cell.names = norm.cell, n.cores = 16,KS.cut = 0.1)
    }}
    
    
    
