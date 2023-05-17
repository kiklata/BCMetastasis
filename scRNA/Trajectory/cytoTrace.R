# add metadata based on NMF metaProgram
## import top5gene and func:ProgramScoring
source('~/scRNA_BC_metastases/Analysis/scRNA/Function_definition/ProgramScoring.R')

library(dplyr)
library(CytoTRACE)
# brain -------------------------------------------------------------------

table(brain$sample,brain$type)
brain.select = subset(brain,sample == 'GSE202501')   
brain.select = brain.select %>% SCTransform(vars.to.regress = c('percent.mt','S.Score','G2M.Score')) %>% 
  RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)

count = as.matrix(brain.select@assays$SCT@counts)
res <- CytoTRACE(mat = count,ncores = 8)


plotCytoTRACE(res, phenotype = brain.select$type, 
              #gene = "ESR1",
              emb = brain.select@reductions$umap@cell.embeddings,outputDir = "~/brain.select")



# primary -----------------------------------------------------------------


primary.score = ProgramScoring(primary)

score = apply(primary.score,2,Stand)

score = as.data.frame(score)

score$type = apply(score,1,
                        FUN = function(x){(colnames(score)[x == max(x)])[1]}) 
primary = AddMetaData(primary,score)

library(CytoTRACE)

count = as.matrix(primary@assays$RNA@counts)

res <- CytoTRACE(mat = count,ncores = 8)

plotCytoTRACE(res, phenotype = primary$type, 
              #gene = "CXCL14",
              #emb = primary@reductions$umap@cell.embeddings
              )
