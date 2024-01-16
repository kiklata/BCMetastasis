cnv_obj <- readRDS("~/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/LiverM/cnv/LiverM2/run.final.infercnv_obj")
cnv = cnv_obj@expr.data

scores=apply(cnv,2,function(x){ sum(x < 0.95 | x > 1.05)/length(x) })
hist(scores, breaks = 500)
table(scores>0.05)






# infercna based ------------------------------------------------------

# ref: An integrative model of cellular states, plasticity, and genetics for glioblastoma. 
library(infercna) # kiklata/infercna
new = cnv -1

df = data.frame(sig = cnaSignal(new,gene.quantile = 0.9),
                cor = cnaCor(new,gene.quantile = 0.9))

cnaScatterPlot(new, gene.quantile = 0.9, main = 'threshold: 0.9')


Modes = infercna::findMalignant(new, gene.quantile = .9, verbose = T)
mal_index = Modes$malignant

anno.matrix$infercnv = if_else((rownames(anno.matrix) %in% mal_index) & (anno.matrix$manual_celltype_annotation == 'Epithelial'),
                               'malignant',anno.matrix$manual_celltype_annotation)

saveRDS(anno.matrix, file = file.path(workdir,'cnv', "mal.anno.matrix.rds"))
