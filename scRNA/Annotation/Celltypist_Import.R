setwd("~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/ilc")

sample.n = dir(pattern = 'predicted')
data.list = list()
for (i in 1:length(sample.n)) {
  data = read.csv(sample.n[i])
  data.list[[i]] = data.frame(cell.id = data$X,celltypist.low = data$predicted_labels)
}

all.data.high = data.table::rbindlist(data.list)
saveRDS(all.data.high,file = '~/ilc.anno.rds')
