library(Seurat)
library(dplyr)
top5gene = readRDS('~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/top5gene.rds')

ProgramScoring = function(sample = sample,gene = top5gene){
  
  seu = sample
  
  seu = AddModuleScore(seu,features = list(gene$metaProgram1),name = 'P8')
  seu = AddModuleScore(seu,features = list(gene$metaProgram2),name = 'P1')
  seu = AddModuleScore(seu,features = list(gene$metaProgram3),name = 'P5')
  seu = AddModuleScore(seu,features = list(gene$metaProgram4),name = 'P9')
  seu = AddModuleScore(seu,features = list(gene$metaProgram5),name = 'P4')
  seu = AddModuleScore(seu,features = list(gene$metaProgram6),name = 'P6')
  seu = AddModuleScore(seu,features = list(gene$metaProgram7),name = 'P2')
  seu = AddModuleScore(seu,features = list(gene$metaProgram8),name = 'P3')
  seu = AddModuleScore(seu,features = list(gene$metaProgram9),name = 'P10')
  seu = AddModuleScore(seu,features = list(gene$metaProgram10),name = 'P7')
  
  metaprogram.score = seu@meta.data[,c('site','study','sample','cell.names',
                                       'P11','P81','P51','P91','P41','P61','P21','P31','P101','P71')]
  colnames(metaprogram.score)[5:14] = c('P1','P8','P5','P9','P4','P6','P2','P3','P10','P7')
  metaprogram.score = metaprogram.score[,5:14]
  return(metaprogram.score)
  
}

Stand = function(data){
  new.data = (data-min(data))/(max(data)-min(data))
  #new.data = if_else(new.data > 0.5,1,0)
}


calcScoreType = function(obj){
  score = try(ProgramScoring(obj))
  if('try-error' %in% class(score)){
    print(obj$sample[1])
  } 
  else{ 
    
  
  score = apply(score,2,Stand)
  score = as.data.frame(score)
  
  score$type = apply(score,1,
                     FUN = function(x){(colnames(score)[x == max(x)])[1]}) 
  return(score)
}}
