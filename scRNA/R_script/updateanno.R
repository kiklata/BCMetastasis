tcell <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/tcell_cluster.tsv", row.names=1) %>% .[c('Tcell_major','Tcell_minor')]
myeloid <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/myeloid_cluster.tsv", row.names=1) %>% .[c('myeloid_major','Myeloid_minor')]
fibro <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/fibro_cluster.tsv", row.names=1) %>% .[c('fibro_major','fibro_minor')]
stro <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/stro_cluster.tsv", row.names=1) %>% .[c('stro_major','stro_minor')]
endo <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/endo_cluster.tsv", row.names=1) %>% .[c('endo_minor')]


colnames(tcell) = c('celltype_major','celltype_minor')
colnames(myeloid) = c('celltype_major','celltype_minor')
colnames(fibro) = c('celltype_major','celltype_minor')
colnames(stro) = c('celltype_major','celltype_minor')
colnames(endo) = c('celltype_minor')

endo$celltype_major = 'Endothelial'

minor = rbind(tcell,myeloid,fibro,stro,endo)

write.table(minor,file = '~/Project/MultiOmics/data/snRNA/Object/summary/annotation/celltype_minor.tsv',sep = '\t',row.names = T,col.names = T,quote = F)
anno <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/anno.tsv", row.names=1)
anno$cor_bar = NULL

anno$id = rownames(anno)
celltype_minor$id = rownames(celltype_minor)

anno_update = full_join(anno,celltype_minor,by = 'id')

anno_update$celltype_major = if_else((anno_update$celltype_major %>% is.na()) == T, anno_update$manual_celltype_annotation,anno_update$celltype_major)
anno_update$celltype_minor = if_else((anno_update$celltype_minor %>% is.na()) == T, anno_update$manual_celltype_annotation,anno_update$celltype_minor)

write.table(anno_update,file = '~/Project/MultiOmics/data/snRNA/Object/summary/annotation/anno.tsv',sep = '\t',row.names = T,col.names = T,quote = F)
