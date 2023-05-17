setwd("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/SCENIC")

tumor <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/tumor.rds")
brain <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/brain.rds")
lymph <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/lymph.rds")
primary <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/primary.rds")
bone <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/bone.rds")
meta.all = tumor@meta.data
meta.all = meta.all[,c(12:14,16,17,27)]

meta.all$primary.subtype = if_else(meta.all$primary.subtype %in% c('ER+','HER2+','Luminal B','TNBC','unspecific'),
                                                                   meta.all$primary.subtype,'unspecific')
library(SCENIC)
aucfile.list = dir(pattern = 'aucell.csv')
sample.list = gsub(".aucell.csv",'',aucfile.list)

aucres.list = list()
for (i in 1:length(aucfile.list)) {
  
aucres = as.data.frame(read_csv(aucfile.list[i]))

aucres.list[[i]] = aucres
names(aucres.list)[i] = sample.list[i]
}

aucres.merge = data.table::rbindlist(aucres.list,fill = T)

type.all.meta = rbind(brain@meta.data,bone@meta.data,lymph@meta.data,primary@meta.data)
type.all.meta = type.all.meta[,c(27:37)]

# TF dimplot -----------------------------------------------------------------

rownames(aucres.merge) = aucres.merge$Cell
tumor = AddMetaData(tumor,aucres.merge)

colnames(tumor@meta.data)[20:245] = gsub(pattern = "\\...",replacement = "(+)",colnames(tumor@meta.data)[20:245])

tumor = tumor %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData(features = rownames(tumor)) %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(dims = 1:30)
library(SCP)
#FeaturePlot(tumor,features = c("AR(+)", "BRCA1(+)", "HDAC2(+)", "SMAD1"))

p1 = ClassDimPlot(
  srt = tumor, group.by = 'site',
  reduction = "UMAP", theme_use = "theme_blank",palette = 'Set1'
)
p2 = ClassDimPlot(
  srt = tumor, group.by = 'study',
  reduction = "UMAP", theme_use = "theme_blank",palette = 'Paired'
)
p3 = ExpDimPlot(
  srt = tumor, features = c("AR(+)", "BRCA1(+)", "HDAC2(+)", "SMAD1(+)"),
  reduction = "UMAP", theme_use = "theme_blank"
)

((p1/p2)|p3) + plot_layout(widths = c(1,2))


# all--------------------
tumor.scenic = tumor[,colnames(tumor) %in% rownames(type.all.meta)]
tumor.scenic = AddMetaData(tumor.scenic,type.all.meta)


p4 = ClassDimPlot(
  srt = tumor.scenic, group.by = 'type',
  reduction = "UMAP", theme_use = "theme_blank",palette = 'Set1'
)

p4

rownames(aucres.merge) = aucres.merge$Cell
aucres.merge$Cell = NULL
aucres.merge = as.matrix(t(aucres.merge))

rss = as.data.frame(calcRSS(AUC = aucres.merge, cellAnnotation = tumor.scenic$type))
rss=na.omit(rss) 

rss = rss[,c('P1','P2','P3','P4','P5',
             'P6','P7','P8','P9','P10')]

fc = apply(rss, 1, sum)-apply(rss, 1, median)  

topgene = names(fc[order(fc,decreasing = T)][1:20])

select.gene = c('AR(+)','CEBPB(+)','ATF3(+)','FOXO1(+)','HES1(+)','HMGA1(+)','RB1(+)','NR2F2(+)',
                'SOX4(+)','SOX9(+)','TCF3(+)','STAT1(+)','BRCA1(+)','E2F1(+)','JUN(+)','NFKB1(+)')

Stand = function(data){
  new.data = (data-min(data))/(max(data)-min(data))
  #new.data = if_else(new.data > 0.5,1,0)
}
plot.rss = rss[select.gene,] %>% Stand() %>% as.matrix()

mycol = rev(corrplot::COL2('RdBu',200))

ComplexHeatmap::pheatmap(plot.rss,
         cluster_rows = F,cluster_cols = F,
         legend = T,
         color = mycol,
         show_colnames = T,angle_col = '0',
         main = 'Tumor')

#bone---------
aucres.merge = data.table::rbindlist(aucres.list,fill = T)

bone.aucres.merge = aucres.merge[aucres.merge$Cell %in% colnames(bone),]
rownames(bone.aucres.merge) = bone.aucres.merge$Cell
bone.aucres.merge$Cell = NULL
bone.aucres.merge = as.matrix(t(bone.aucres.merge))
rss = as.data.frame(calcRSS(AUC = bone.aucres.merge, cellAnnotation = bone$type))
rss=na.omit(rss) 

rss = rss[,c('P1','P2','P3','P4','P5',
             'P6','P7','P8','P9','P10')]

fc = apply(rss, 1, sum)-apply(rss, 1, median)  

topgene = names(fc[order(fc,decreasing = T)][1:20])

select.gene = c('AR(+)','CEBPB(+)','ATF3(+)','FOXO1(+)','HES1(+)','HMGA1(+)','RB1(+)','NR2F2(+)',
                'SOX4(+)','SOX9(+)','TCF3(+)','STAT1(+)','BRCA1(+)','E2F1(+)','JUN(+)','NFKB1(+)')

plot.rss = rss[select.gene,] %>% Stand() %>% as.matrix()

mycol = rev(corrplot::COL2('RdBu',n = 200))
ComplexHeatmap::pheatmap(plot.rss,
         cluster_rows = F,cluster_cols = F,
         legend = T, color = mycol,
         show_colnames = T,angle_col = '0',
         main = 'Bone')

#brain
aucres.merge = data.table::rbindlist(aucres.list,fill = T)

brain.aucres.merge = aucres.merge[aucres.merge$Cell %in% colnames(brain),]
rownames(brain.aucres.merge) = brain.aucres.merge$Cell
brain.aucres.merge$Cell = NULL
brain.aucres.merge = as.matrix(t(brain.aucres.merge))
rss = as.data.frame(calcRSS(AUC = brain.aucres.merge, cellAnnotation = brain$type))
rss=na.omit(rss) 

rss = rss[,c('P1','P2','P3','P4','P5',
             'P6','P7','P8','P9','P10')]

fc = apply(rss, 1, sum)-apply(rss, 1, median)  

topgene = names(fc[order(fc,decreasing = T)][1:20])

select.gene = c('AR(+)','CEBPB(+)','ATF3(+)','FOXO1(+)','HES1(+)','HMGA1(+)','RB1(+)','NR2F2(+)',
                'SOX4(+)','SOX9(+)','TCF3(+)','STAT1(+)','BRCA1(+)','E2F1(+)','JUN(+)','NFKB1(+)')

plot.rss = rss[select.gene,] %>% Stand() %>% as.matrix()

mycol = rev(corrplot::COL2('RdBu',n = 200))
ComplexHeatmap::pheatmap(plot.rss,
                         cluster_rows = F,cluster_cols = F,
                         legend = T, color = mycol,
                         show_colnames = T,angle_col = '0',
                         main = 'Brain')

#lymph----------------
aucres.merge = data.table::rbindlist(aucres.list,fill = T)

lymph.aucres.merge = aucres.merge[aucres.merge$Cell %in% colnames(lymph),]
rownames(lymph.aucres.merge) = lymph.aucres.merge$Cell
lymph.aucres.merge$Cell = NULL
lymph.aucres.merge = as.matrix(t(lymph.aucres.merge))
rss = as.data.frame(calcRSS(AUC = lymph.aucres.merge, cellAnnotation = lymph$type))
rss=na.omit(rss) 

rss = rss[,c('P1','P2','P3','P4','P5',
             'P6','P7','P8','P9','P10')]

fc = apply(rss, 1, sum)-apply(rss, 1, median)  

topgene = names(fc[order(fc,decreasing = T)][1:20])

select.gene = c('AR(+)','CEBPB(+)','ATF3(+)','FOXO1(+)','HES1(+)','HMGA1(+)','RB1(+)','NR2F2(+)',
                'SOX4(+)','SOX9(+)','TCF3(+)','STAT1(+)','BRCA1(+)','E2F1(+)','JUN(+)','NFKB1(+)')

plot.rss = rss[select.gene,] %>% Stand() %>% as.matrix()

mycol = rev(corrplot::COL2('RdBu',n = 200))
ComplexHeatmap::pheatmap(plot.rss,
                         cluster_rows = F,cluster_cols = F,
                         legend = T, color = mycol,
                         show_colnames = T,angle_col = '0',
                         main = 'lymph')
#primary------------------

aucres.merge = data.table::rbindlist(aucres.list,fill = T)

primary.aucres.merge = aucres.merge[aucres.merge$Cell %in% colnames(primary),]
rownames(primary.aucres.merge) = primary.aucres.merge$Cell
primary.aucres.merge$Cell = NULL

primary.aucres.merge = as.matrix(t(primary.aucres.merge))

rss = as.data.frame(calcRSS(AUC = primary.aucres.merge, cellAnnotation = primary$type))
rss=na.omit(rss) 

rss = rss[,c('P1','P2','P3','P4','P5',
             'P6','P7','P8','P9','P10')]

fc = apply(rss, 1, sum)-apply(rss, 1, median)  

topgene = names(fc[order(fc,decreasing = T)][1:20])

select.gene = c('AR(+)','CEBPB(+)','ATF3(+)','FOXO1(+)','HES1(+)','HMGA1(+)','RB1(+)','NR2F2(+)',
                'SOX4(+)','SOX9(+)','TCF3(+)','STAT1(+)','BRCA1(+)','E2F1(+)','JUN(+)','NFKB1(+)')

plot.rss = rss[select.gene,] %>% Stand() %>% as.matrix()

mycol = rev(corrplot::COL2('RdBu',n = 200))
ComplexHeatmap::pheatmap(plot.rss,
                         cluster_rows = F,cluster_cols = F,
                         legend = T, color = mycol,
                         show_colnames = T,angle_col = '0',
                         main = 'Primary site')
