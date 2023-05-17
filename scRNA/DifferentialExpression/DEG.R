library(Seurat)
library(COSG)

tumor <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds")

Idents(tumor) = tumor$site
marker.brain = FindMarkers(tumor,ident.1 = 'brain',ident.2 = 'primary')
marker.lymph = FindMarkers(tumor,ident.1 = 'lymph',ident.2 = 'primary')
marker.bone = FindMarkers(tumor,ident.1 = 'bone',ident.2 = 'primary')

save(marker.brain,marker.lymph,marker.bone,file = 'find.marker.rdata')

#cosg--------------------

brain = subset(tumor,site %in% c('primary','brain'))
cosg.brain = cosg(brain,groups = 'all',assay = 'RNA')

lymph = subset(tumor,site %in% c('lymph','primary'))
cosg.lymph = cosg(lymph,groups = 'all',assay = 'RNA')

bone = subset(tumor,site %in% c('bone','primary'))
cosg.bone = cosg(bone,groups = 'all',assay = 'RNA')

save(cosg.brain,cosg.lymph,cosg.brain,file = 'cosg.marker.rdata')