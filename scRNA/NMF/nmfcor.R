library(NMF)
setwd("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/SCT_100cell_RNMF")

all.res = dir(pattern = 'res.list.rds')

# func --------------------------------------------------------------------

NMFtoGenes.CELL = function(nmf.res,rank.select = 15,genenum = 50){
  
  res = nmf.res[[as.character(rank.select)]]
  fs = extractFeatures(res, genenum) 
  fs = lapply(fs, function(x) rownames(res)[x]) 
  fs = do.call("rbind", fs) 
  rownames(fs) = paste0("cluster", 1:as.numeric(rank.select))
  fs = as.data.frame(t(fs))
  
  return(fs)
  
}

NMFtoWeight.CELL = function(nmf.res,rank.select = 15){
  
  res = nmf.res[[as.character(rank.select)]]
  fs = res@fit@W
  fs = as.data.frame(fs)
  colnames(fs) = paste0('cluster',1:rank.select)
  
  return(fs) 
  
}


# import all modual -------------------------------------------------------

sign.cell.weight = list()

sign.cell.gene = list()

for (i in 1:length(all.res)) {
  
  nmf.res <- readRDS(all.res[i])
  
  nmf.res.weight = NMFtoWeight.CELL(nmf.res,rank.select = 10)
  sign.cell.weight[[i]] = nmf.res.weight
  names(sign.cell.weight)[i] = all.res[i]
  
  nmf.res.gene = NMFtoGenes.CELL(nmf.res,rank.select = 10,genenum = 50)
  colnames(nmf.res.gene) = paste0(all.res[i],colnames(nmf.res.gene))
  sign.cell.gene[[i]] = nmf.res.gene
  names(sign.cell.gene)[i] = all.res[i]
  
}

gene = rownames(sign.cell.weight[[1]])
gene = gene[!duplicated(gene)]

for (i in 2:length(sign.cell.gene)) {
  

gene1 = rownames(sign.cell.weight[[i]])

gene = intersect(gene,gene1)

#new.weight = sign.cell.weight[[i]][gene,]

}


new.weight = sign.cell.weight[[1]][gene,]
colnames(new.weight) = paste0(names(sign.cell.weight[1]),colnames(new.weight))
for (i in 2:length(sign.cell.weight)) {
  new.weight2 = sign.cell.weight[[i]][gene,]
  colnames(new.weight2) = paste0(names(sign.cell.weight[i]),colnames(new.weight2))
  new.weight = cbind(new.weight,new.weight2)
}

#library(vegan)

#distance = vegdist(t(new.weight),method = 'jaccard')
#distance.mat = as.matrix(distance)
#corrplot(distance.mat,method = 'color',order = 'hclust',hclust.method = 'ward.D2',
#         addrect = 8,tl.pos = 'n',tl.cex = 0.5,col = mycol,col.lim = c(0,1),cl.length = 6 )
#final = hclust(distance,method = 'ward.D2')
#table(cutree(final,k = 10))
#table(cutree(final,h = 2.2))


# metaProgram generate ----------------------------------------------------

cor.mat = cor(new.weight,method = 'pearson')
min(cor.mat)
max(cor.mat)

k.hclust = 8

library(corrplot)
mycol = rev(COL2('RdBu',200))
corrplot(cor.mat,method = 'color',order = 'hclust',hclust.method = 'ward.D2',
         addrect = k.hclust,tl.pos = 'n',tl.cex = 0.5,cl.pos = 'b',
         col = mycol,col.lim = c(-1,1),cl.length = 6)


hc = hclust(as.dist(1 - cor.mat), method = 'ward.D2')
tree = as.data.frame(cutree(hc,k = k.hclust))
table(tree$`cutree(hc, k = k.hclust)`)

clust.index = order.dendrogram(as.dendrogram(hc))
clust.modual = hc$labels[clust.index]
tree$modual = rownames(tree)
colnames(tree)[1] = 'cluster'

library(tidyverse)
tree = tree %>% mutate(modual = fct_relevel(modual,clust.modual)) %>% arrange(modual)
saveRDS(tree,file = '~/scRNA_BC_metastases/Result/NMF.cluster.res.rds')

# metaProgram annotate ----------------------------------------------------
genemodual = list()

for (i in 1:length(sign.cell.gene)) {
  
  genelist = data.frame()
  
  for (k in 1:ncol(sign.cell.gene[[i]])) {
    
    genelist[k,1:5] = sign.cell.gene[[i]][,k][1:5]
    rownames(genelist)[k] = colnames(sign.cell.gene[[i]])[k]
    
  }
  colnames(genelist)[1:5] = paste0('gene',1:5)
  genelist$modual = rownames(genelist)
  genemodual[[i]] = genelist
  names(genemodual)[i] = names(sign.cell.gene)[i]
}

gene.res = data.table::rbindlist(genemodual)
NMF.cluster.res = tree
#NMF.cluster.res = readRDS("~/scRNA_BC_metastases/Result/NMF.cluster.res.rds")
NMF.cluster.res = left_join(NMF.cluster.res,gene.res,by = 'modual')

go.res.list = data.frame()
gene.res.list = list()

# annotate with top3 gene for each modual ---------------------------------


library(clusterProfiler)
library(org.Hs.eg.db)

for(k in 1:k.hclust){
  
  genetoanno = filter(NMF.cluster.res,cluster == k)
  genetoanno = as.character(rbind(genetoanno$gene1,genetoanno$gene2,genetoanno$gene3))
  
  genetoanno = genetoanno[!grepl('^RPL',genetoanno)]
  genetoanno = genetoanno[!grepl('^RPS',genetoanno)]
  genetoanno = genetoanno[!grepl('^MT',genetoanno)]
  gene.res.list[[k]] = genetoanno
  names(gene.res.list)[k] = paste0('metaProgram',k)
  
  genelist.ENTREZID = bitr(genetoanno,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)[,2]
  
  ego = enrichGO(genelist.ENTREZID,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",
                 ont = "ALL",pvalueCutoff = 0.05,pAdjustMethod = "BH",
                 #universe = ,qvalueCutoff = 0.2,
                 minGSSize = 10,maxGSSize = 500,readable = FALSE,pool = FALSE)
  
  
  result = ego@result
  result = result[order(result$GeneRatio,decreasing = T),]

  go.res.list[k,1:3] = filter(result,ONTOLOGY =='BP')$Description[1:3]
  go.res.list[k,4:6] = filter(result,ONTOLOGY =='CC')$Description[1:3]
  go.res.list[k,7:9] = filter(result,ONTOLOGY =='MF')$Description[1:3]
  colnames(go.res.list)[1:3] = paste0('Biological Process',1:3)
  colnames(go.res.list)[4:6] = paste0('Cellular Component',1:3)
  colnames(go.res.list)[7:9] = paste0('Molecular Function',1:3)
  
  rownames(go.res.list)[k] = paste0('metaProgram',k)
}

# top5 gene for each modual ---------------------------------

top5gene = list()

for(k in 1:k.hclust){
  
  genetoanno = filter(NMF.cluster.res,cluster == k)
  genetoanno = as.character(cbind(genetoanno$gene1,genetoanno$gene2,
                                  genetoanno$gene3,genetoanno$gene4,
                                  genetoanno$gene5))
  
  genetoanno = genetoanno[!grepl('^RPL',genetoanno)]
  genetoanno = genetoanno[!grepl('^RPS',genetoanno)]
  genetoanno = genetoanno[!grepl('^MT',genetoanno)]
  genetoanno = genetoanno[!duplicated(genetoanno)]
  top5gene[[k]] = genetoanno
  names(top5gene)[k] = paste0('metaProgram',k)
  
}
saveRDS(top5gene,file = '~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/top5gene.rds')
fig.order = NMF.cluster.res$cluster[!duplicated(NMF.cluster.res$cluster)]

write.csv(go.res.list,file = '~/go.res.list.csv')


# heatmap of different site -----------------------------------------------

metaProgram.gene = top1gene$metaProgram1
names(metaProgram.gene)[1:length(top1gene[[1]])] = names(top1gene)[1]
  
for(i in 2:k.hclust){
  geneToadd = top1gene[[i]]
  names(geneToadd)[1:length(top1gene[[i]])] = names(top1gene)[i]
  metaProgram.gene = append(metaProgram.gene,geneToadd)
}

metaProgram.gene = as.data.frame(metaProgram.gene,row.names = names(metaProgram.gene))
metaProgram.gene$program = rownames(metaProgram.gene)
saveRDS(metaProgram.gene,file = 'metaProgram.gene.rds')

library(Seurat)
library(dplyr)

metaProgram.gene = readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/SCT_100cell_RNMF/metaProgram.gene.rds")
tumor = readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds")

metaProgram.gene$program = factor(metaProgram.gene$program,
                                     levels = c('metaProgram2','metaProgram7','metaProgram8',
                                                'metaProgram5','metaProgram3','metaProgram6',
                                                'metaProgram10','metaProgram1','metaProgram4',
                                                'metaProgram9'))

metaProgram.gene = metaProgram.gene[order(metaProgram.gene$program),]
rownames(metaProgram.gene) = NULL
metaProgram.gene$plot.pro = if_else(metaProgram.gene$program == 'metaProgram2','P1',
                                    if_else(metaProgram.gene$program == 'metaProgram7','P2',
                                    if_else(metaProgram.gene$program == 'metaProgram8','P3',
                                    if_else(metaProgram.gene$program == 'metaProgram5','P4',
                                    if_else(metaProgram.gene$program == 'metaProgram3','P5',
                                    if_else(metaProgram.gene$program == 'metaProgram6','P6',
                                    if_else(metaProgram.gene$program == 'metaProgram10','P7',
                                    if_else(metaProgram.gene$program == 'metaProgram1','P8',
                                    if_else(metaProgram.gene$program == 'metaProgram4','P9',
                                    if_else(metaProgram.gene$program == 'metaProgram9','P10','na'))))))))))



tumor = tumor %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(features = rownames(tumor))

library(ComplexHeatmap)
library(RColorBrewer)


# modual score ------------------------------------------------------------

# scoring by sample


top5gene = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/top5gene.rds')
tumor = readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds")

source('~/scRNA_BC_metastases/Analysis/scRNA/Function_definition/ProgramScoring.R')

metaprogram.score = ProgramScoring(tumor,top5gene)

saveRDS(metaProgramScore.list,file = 'metaProgram.score.rds')


#tumor = tumor %>% FindVariableFeatures() %>%
#  NormalizeData() %>% ScaleData(features = rownames(tumor)) %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
#  FindClusters() %>% RunUMAP(dims = 1:30)

# UMAP of metaProgram
tumor = readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds")

table(tumor$sample)

# pick a example
source("~/scRNA_BC_metastases/Analysis/scRNA/Function_definition/plotScore.R", echo=TRUE)
# brain cell.pid2

seu = subset(tumor,sample == 'cell.pid2')
seu = CreateSeuratObject(seu[['RNA']]@counts,meta.data = seu@meta.data,min.cells = 3,min.features = 200)

brain.score = ProgramScoring(seu,gene = top5gene)

seu.brain = seu %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(features = rownames(seu)) %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

score = brain.score[,5:14]

score$P1 = Stand(score$P1)
score$P2 = Stand(score$P2)
score$P3 = Stand(score$P3)
score$P4 = Stand(score$P4)
score$P5 = Stand(score$P5)
score$P6 = Stand(score$P6)
score$P7 = Stand(score$P7)
score$P8 = Stand(score$P8)
score$P9 = Stand(score$P9)
score$P10 = Stand(score$P10)

seu.brain = AddMetaData(seu.brain,score)

summary(seu.brain$P3)
  
plotScore(seu = seu.brain,savefile = "~/brain.pdf")

# lymph EMBOJ.GSM4909312
seu = subset(tumor,sample == 'EMBOJ.GSM4909312')
seu = CreateSeuratObject(seu[['RNA']]@counts,meta.data = seu@meta.data,min.cells = 3,min.features = 200)

lymph.score = ProgramScoring(seu,gene = top5gene)

seu.lymph = seu %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(features = rownames(seu)) %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

score = lymph.score[,5:14]

score$P1 = Stand(score$P1)
score$P2 = Stand(score$P2)
score$P3 = Stand(score$P3)
score$P4 = Stand(score$P4)
score$P5 = Stand(score$P5)
score$P6 = Stand(score$P6)
score$P7 = Stand(score$P7)
score$P8 = Stand(score$P8)
score$P9 = Stand(score$P9)
score$P10 = Stand(score$P10)

seu.lymph = AddMetaData(seu.lymph,score)

summary(seu.lymph$P3)

plotScore(seu = seu.lymph,savefile = "~/lymph.pdf")


# bone bom2
seu = subset(tumor,sample == 'BoM2')
seu = CreateSeuratObject(seu[['RNA']]@counts,meta.data = seu@meta.data,min.cells = 3,min.features = 200)

bone.score = ProgramScoring(seu,gene = top5gene)

seu.bone = seu %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(features = rownames(seu)) %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

score =bone.score[,5:14]

score$P1 = Stand(score$P1)
score$P2 = Stand(score$P2)
score$P3 = Stand(score$P3)
score$P4 = Stand(score$P4)
score$P5 = Stand(score$P5)
score$P6 = Stand(score$P6)
score$P7 = Stand(score$P7)
score$P8 = Stand(score$P8)
score$P9 = Stand(score$P9)
score$P10 = Stand(score$P10)

seu.bone = AddMetaData(seu.bone,score)

summary(seu.bone$P3)

plotScore(seu = seu.bone,savefile = "~/bone.pdf")


# primary primary.CID4067
seu = subset(tumor,sample == 'primary.CID4067')
seu = CreateSeuratObject(seu[['RNA']]@counts,meta.data = seu@meta.data,min.cells = 3,min.features = 200)

primary.score = ProgramScoring(seu,gene = top5gene)

seu.primary = seu %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(features = rownames(seu)) %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

score =primary.score[,5:14]

score$P1 = Stand(score$P1)
score$P2 = Stand(score$P2)
score$P3 = Stand(score$P3)
score$P4 = Stand(score$P4)
score$P5 = Stand(score$P5)
score$P6 = Stand(score$P6)
score$P7 = Stand(score$P7)
score$P8 = Stand(score$P8)
score$P9 = Stand(score$P9)
score$P10 = Stand(score$P10)

seu.primary = AddMetaData(seu.primary,score)

summary(seu.primary$P3)

plotScore(seu = seu.primary,savefile = "~/primary.pdf")


# boxplot for all metaProgram

# calc for each sample

tumor = readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds")
top5gene = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/top5gene.rds')

sample.n = names(table(tumor$sample))

Stand = function(data){
  new.data = (data-min(data))/(max(data)-min(data))
  new.data = if_else(new.data > 0.5,1,0)
}

data = data.frame()

for (i in 1:length(sample.n)) {
  
seu = subset(tumor,sample == sample.n[i])
seu = CreateSeuratObject(seu[['RNA']]@counts,meta.data = seu@meta.data,min.cells = 3,min.features = 200)

score = try(ProgramScoring(seu,gene = top5gene))

if('try-error' %in% class(score)){
  
  data[i,1:10] = sample.n[i]} 
  
else{
    
  score$P1 = Stand(score$P1)
  data[i,1] = sum(score$P1)/length(score$P1)
  data[i,2] = sum(score$P2)/length(score$P2)
  data[i,3] = sum(score$P3)/length(score$P3)
  data[i,4] = sum(score$P4)/length(score$P4)
  data[i,5] = sum(score$P5)/length(score$P5)
  data[i,6] = sum(score$P6)/length(score$P6)
  data[i,7] = sum(score$P7)/length(score$P7)
  data[i,8] = sum(score$P8)/length(score$P8)
  data[i,9] = sum(score$P9)/length(score$P9)
  data[i,10] = sum(score$P10)/length(score$P10)
  
  }
}
rownames(data) = sample.n
colnames(data)[1:10] = paste0('P',1:10)

saveRDS(data,file = '~/propor.rds')
