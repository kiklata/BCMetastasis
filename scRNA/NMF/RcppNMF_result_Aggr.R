
setwd("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/rcpp")
dir.list = dir(recursive = F)
for (i in 1:length(dir.list)) {
file.list = dir(dir.list[i])
setwd(dir.list[i])
rank.list = list()
for (k in 1:length(file.list)) {
nmf = readRDS(file.list[k])
rank.list[[k]] = nmf
names(rank.list)[k] = file.list[k]
}
setwd("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/rcpp")
saveRDS(rank.list,file = paste0(dir.list[i],'.rank.rds'))
}


name.list = dir(path = 'rank2')

file.list = dir(pattern = '.rds')

for (k in 1:length(name.list)) {

nmf.list = list()

for (i in 1:length(file.list)) {
  nmf = readRDS(file.list[i])  
  nmf.list[[i]] = nmf[[name.list[[k]]]]
  names(nmf.list)[i] = file.list[i]
}

saveRDS(nmf.list,file = paste0(name.list[k]))
}

# after remove previous 'rank' filefold

file.list = dir(pattern = '.rds')

for (i in 1:length(file.list)) {
  
nmf = readRDS(file.list[i])

new.name = character()

for (a in 1:length(names(nmf))) {
  
  new.name[a] = strsplit(names(nmf)[a],split = '\\.')[[1]][1]
  
}
names(nmf) = new.name
saveRDS(nmf,file = file.list[i])
}
