library(compassR)
library(Seurat)

seu <- readRDS("~/PaperCD8/data/Tex/CD8.Tex.harmony.rds")
after = subset(seu,sample.timepoint == 'After')

a_sample = names(table(after$sample.ID))

count_list = list()
k = 1
for (i in a_sample) {
  seuobj = subset(after,sample.ID == i)
  seuobj = NormalizeData(seuobj,normalization.method = 'RC',scale.factor = 1e4)
  count = as.data.frame(seuobj@assays$RNA@data)
  count$symbol = rownames(count)
  rownames(count) = NULL
  count_list[[k]] = count
  k = k+1
}

finalcount = count_list[[1]]
library(dplyr)
for (i in 2:length(count_list)) {
  countnow = count_list[[i]]
  finalcount = left_join(finalcount,countnow,by = 'symbol')
} 

colnames(finalcount) = gsub('\\.','-',colnames(finalcount))

a_cell_metadata = after@meta.data
rownames(a_cell_metadata) = NULL
a_cell_metadata$cell.id = gsub('\\.','-', a_cell_metadata$cell.id)

a_cell_metadata = a_cell_metadata[order(a_cell_metadata$cell.id),]
newc = c('symbol',a_cell_metadata$cell.id)
finalcount = finalcount[,newc]

write.table(finalcount,file = 'compassTest/after/linear_gene_expression_matrix.tsv',
            quote = F,row.names = F,sep = '\t')
write.csv(a_cell_metadata,file = 'compassTest/after/cell_metadata.csv',row.names = F)


compass.after <- readRDS("~/compassTest/after/compass.after.rds")
rownames(compass.after) = compass.after$X
colnames(compass.after) = gsub('\\.','-',colnames(compass.after))
compass.after = compass.after[,colnames(compass.after) %in% colnames((finalcount))]
compass.after = compass.after[,a_cell_metadata$cell.id]
write.table(compass.after,file = 'compassTest/after/reactions.tsv',quote = F,row.names = T,sep = '\t')

colnames(reactions) = gsub('\\.','-',colnames(reactions))
reactions$reaction = rownames(reactions)
newc = colnames(reactions)
newc[length(newc)] = newc[1]
newc[1] = 'reaction'
reactions = reactions[,newc]

write.table(reactions,file = 'before/reactions.tsv',quote = F,row.names = F,sep = '\t')

