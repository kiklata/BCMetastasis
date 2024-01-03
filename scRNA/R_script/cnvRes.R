library(infercnv)
options(scipen = 100)
library(infercna)# ref: An integrative model of cellular states, plasticity, and genetics for glioblastoma. 

workdir = '~/Project/MultiOmics/data/snRNA/Object/summary'

sample_list = c('P1013S2','P1015S2','P1018S1')
set.seed(42)

scores.list = list()
cna_mal.list = list()

for (i in sample_list) {

  run.final.infercnv_obj <- readRDS(paste0(workdir,"/cnv/",i,"/run.final.infercnv_obj"))
  cnv = run.final.infercnv_obj@expr.data
  
  scores.list[[i]]=apply(cnv,2,function(x){ sum(x < 0.95 | x > 1.05)/length(x) })
  
  df = data.frame(sig = cnaSignal(cnv,gene.quantile = 0.9),
                  cor = cnaCor(cnv,gene.quantile = 0.9))
  
  #cnaScatterPlot(cnv, gene.quantile = 0.9, main = 'threshold: 0.9')
  
  Modes = infercna::findMalignant(cnv - 1, gene.quantile = .9, verbose = T,plot = F)
  cna_mal.list[[i]] = Modes[['malignant']]

}

mal_index = cna_mal.list %>% unlist() %>% unname()

anno <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/cellranger_filter_anno.tsv", row.names=1)
anno$infercnv = if_else((rownames(anno) %in% mal_index) & (anno$manual_celltype_annotation == 'Epithelial'),
                               'malignant',anno$manual_celltype_annotation)

write.table(anno,file = 'infercnv.tsv',sep = '\t',quote = F,col.names = T,row.names = T)
#saveRDS(anno.matrix, file = file.path(workdir,'cnv', "mal.anno.matrix.rds"))
