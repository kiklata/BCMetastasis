#SciBet

setwd("~/software/scibet/count")
library(scuttle)

file.path = dir()
file.path = file.path[c(-3,-4,-7,-8,-21,-22,-23)]

count.list = list()

for (i in 1:length(file.path)) {
  
  count = read.delim( file.path[i], row.names=1,header = T)
  #tpm = calculateTPM(count)
  
  count = as.data.frame(count)
  count$gene = rownames(count)
  
  count.list[[i]] = count
  names(count.list)[i] = file.path[i]
}


saveRDS(count.list, file = 'count.list.rds')

# subset cell for immune altas SciBet


#devtools::install_github("PaulingLiu/scibet")
library(scibet)

data = count.list[[1]]

for (i in 2:length(count.list)) {

  data1 = count.list[[i]]

  data = left_join(data,data1, by = 'gene')
  
  }

rownames(data) = data$gene
saveRDS(data,'data.rds')
seu = CreateSeuratObject(counts = data,meta.data = GSE156728_metadata.txt )
saveRDS(seu,file = 'seu.rds')

cluster.n = names(table(seu$meta.cluster))
# for cell>1000 subs for 1000cells
su.list = list()
for (i in 1:length(cluster.n)) {
  su = subset(seu,meta.cluster == cluster.n[i])
  if(dim(su)[2]>1000){
    subset.cell = sample(colnames(su),size = 1000)
    su = su[,subset.cell]
  }else{
      su = su
    }
  
  su.list[[i]] = su
  names(su.list)[i] = cluster.n[i]
    }
su.merge = merge(su.list[[1]],c(su.list[2:41]))

su.merge = as.SingleCellExperiment(su.merge)
tpm.ref = calculateTPM(su.merge)

tpm.ref.t = t(tpm.ref)

tpm.ref.t[,ncol(tpm.ref.t)+1] = rownames(tpm.ref.t)
# add label
saveRDS(SciBet.ref,file = '~/software/scibet/count/SciBet.ref.rds')
# save as ref 


# run SciBet
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))

SciBet.ref <- readRDS("~/software/scibet/count/SciBet.ref.rds")
tcell <- readRDS("~/tcell.rds")

expr = tcell[['RNA']]@counts
expr = scuttle::calculateTPM(expr)
expr = as.matrix(expr)
expr = t(expr)

#etest_gene <- SelectGene(SciBet.ref, k = 50)
#Marker_heatmap(SciBet.ref, etest_gene)

prd <- SciBet(SciBet.ref, expr)
prd = as.data.frame(prd)
prd$cellid = rownames(expr)

colnames(prd) = c('scibet.celltype.subset','cell.id')
rownames(prd) = prd$cell.id
prd$cell.id = NULL

tcell = AddMetaData(tcell,prd)

saveRDS(tcell.anno,file = '~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/immune.anno.rds')

# regroup

cluster.n = names(table(tcell.anno$scibet.celltype.subset))

CD4.n = cluster.n[c(1:24)]
CD8.n = cluster.n[c(25:41)]

cd4.n.n = CD4.n[1:4]
cd4.m.n = CD4.n[5:11]
cd4.em.n = CD4.n[12]
cd4.emra.n = CD4.n[13]
cd4.th17.n = CD4.n[c(14,15)]
cd4.tfh.n = CD4.n[16]
cd4.tfhth1.n = CD4.n[17]
cd4.treg.n = CD4.n[18:21]
cd4.isg.n = CD4.n[22]
cd4.mix.n = CD4.n[c(23,24)]

cd8.n.n = CD8.n[1]
cd8.m.n = CD8.n[c(2:4,17)]
cd8.em.n = CD8.n[c(5,6)]
cd8.emra.n = CD8.n[7]
cd8.k.n = CD8.n[c(8,9)]
cd8.rm.n = CD8.n[10]
cd8.ex.n = CD8.n[11:14]
cd8.isg.n = CD8.n[15]
cd8.mait.n = CD8.n[16]

tcell.anno$scibet.celltype.minor = if_else(tcell.anno$scibet.celltype.subset %in% CD4.n,'CD4 T cell','CD8 T cell')

tcell.anno$scibet.celltype.marker = tcell.anno$scibet.celltype.subset

tcell.anno$scibet.celltype.subset = if_else(tcell.anno$scibet.celltype.marker %in% cd4.n.n,'CD4 naive T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.m.n,'CD4 memory T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.em.n,'CD4 effector memory T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.emra.n,'CD4 terminally differentiated effector memory or effector cells T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.th17.n,'CD4 helper T cell 17',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.tfh.n,'CD4 follicular helper T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.tfhth1.n,'Tfh/Th1 dual-functional T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.treg.n,'CD4 T regulatory cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.isg.n,'CD4 interferon-stimulated genes related T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd4.mix.n,'CD4 mix T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.n.n ,'CD8 naive T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.m.n ,'CD8 memory T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.em.n,'CD8 effector memory T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.emra.n,'CD8 terminally differentiated effector memory or effector cells T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.k.n,'CD8 NK-like T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.rm.n,'CD8 tissue-resident memory T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.ex.n,'CD8 exhausted T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.isg.n,'CD8 interferon-stimulated genes related T cell',
                                            if_else(tcell.anno$scibet.celltype.marker %in% cd8.mait.n,'CD8 mucosal-associated invariant T cells','na')))))))))))))))))))

table(tcell$celltypist.celltype.minor,tcell$scibet.celltype.minor)
