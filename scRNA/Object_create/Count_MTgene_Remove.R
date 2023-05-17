all.count = all@assays$RNA@counts
all.count[,'MT.CO3']

firstmt.name = grep ("^MT-", rownames(all.count),value = T)
secondmt.name = grep ("^MT.", rownames(all.count),value = T)[177:189]
new.name = gsub('MT.','MT-',secondmt.name)

for (i in  1:13) {
  
  value.1 = all.count[new.name[i],]
  value.2 = all.count[secondmt.name[i],]
  value = value.1 + value.2
  names(value) = NULL
  
  all.count = all.count[rownames(all.count) != secondmt.name[i],]
  all.count = all.count[rownames(all.count) != new.name[i],]
  
  rrr = matrix(value,nrow = 1,ncol = length(value),dimnames = list(new.name[i],names(value)))
  
  all.count = rbind(all.count,rrr)
  
  print(i)
}


all.new = CreateSeuratObject(all.count)

all.anno <- readRDS("~/scRNA_BC_metastases/merge/anno/all.anno.rds")
all.anno$percent.mt = NULL
saveRDS(all.anno,file = 'all.anno.rds')
