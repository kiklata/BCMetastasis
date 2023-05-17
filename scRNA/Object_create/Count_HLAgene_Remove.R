all.count = all@assays$RNA@counts
all.count['HLA.A',]

firstHLA.name = grep("^HLA-", rownames(all.count),value = T)
secondHLA.name = grep("^HLA.", rownames(all.count),value = T)[32:50]
new.name = gsub('HLA.','HLA-',secondHLA.name)

for (i in  1:19) {
  
  value.1 = all.count[new.name[i],]
  value.2 = all.count[secondHLA.name[i],]
  value = value.1 + value.2
  names(value) = NULL
  
  all.count = all.count[rownames(all.count) != secondHLA.name[i],]
  all.count = all.count[rownames(all.count) != new.name[i],]
  
  rrr = matrix(value,nrow = 1,ncol = length(value),dimnames = list(new.name[i],names(value)))
  
  all.count = rbind(all.count,rrr)
  
  print(i)
}


all.new = CreateSeuratObject(all.count)

saveRDS(all.new,file = '~/scRNA_BC_metastases/Data/integrate_sc_BC/all/all.rds')
