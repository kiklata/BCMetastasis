library(immunarch)

filepath = c("~/BCLM/mixcr/clonetype/")
#create metadata
mixcr = repLoad(filepath)
saveRDS(mixcr,file = 'mixcr_res.rds')

# add barcode