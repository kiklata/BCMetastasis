library(IOBR)
library(tidyr)

shiny::runApp('~/software/IOBRshiny')

propotion = read.csv("~/scRNA_BC_metastases/Data/bulk_BC/MET500/IOBR/cibersort_abs.csv")
colnames(mBC_metadata)[1] = 'ID'

propotion = left_join(propotion,mBC_metadata,by = 'ID')

names(table(mBC_metadata$biopsy_tissue))

boneM = filter(propotion,biopsy_tissue == 'bone_marrow')[,c(1:23)]
brain = filter(propotion,biopsy_tissue == 'brain')[,c(1:23)]
liver = filter(propotion,biopsy_tissue == 'liver')[,c(1:23)]
lung = filter(propotion,biopsy_tissue == 'lung')[,c(1:23)]
lymphNode = filter(propotion,biopsy_tissue == 'lymph_node')[,c(1:23)]
skin = filter(propotion,biopsy_tissue == 'skin')[,c(1:23)]

#save(boneM,brain,liver,lung,lymphNode,skin,file = 'cibersort_by_site.rdata')

cell_bar_plot(skin)
