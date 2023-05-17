data.dir = '~/scRNA_BC_metastases/Data/integrate_sc_BC/epi/copyKat'
dir.list = list.dirs(data.dir,recursive = F)

sample.list = list.dirs(data.dir,recursive = F,full.names = F)

library(dplyr)

pred.list = list()

for (i in 1:length(dir.list)) {
  
  dir.list[i]
  prediction.file = list.files(dir.list[i],pattern = '*prediction.txt',full.names = T)
  
  pred.res = read.delim(prediction.file)
  pred.epi.res = filter(pred.res,copykat.pred != 'not.defined')
  pred.list[[i]] = data.frame(sample = sample.list[i],pred.epi.res)

}

pred = data.table::rbindlist(pred.list)

saveRDS(pred,file = '~/scRNA_BC_metastases/Data/integrate_sc_BC/epi/copyKat.res.rds')
