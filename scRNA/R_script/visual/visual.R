
# 0.1 pkg load ----------------------------------------------------------

library(Seurat)
library(dplyr)

library(ggplot2)
library(ggplotify)
library(ggthemes)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(ggsci)
library(scCustomize) # dimplot feature plot
library(patchwork)
#library(showtext)
#font_add('arial',regular = '~/software/fonts/arial.ttf')
#showtext_auto()

library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 2400 * 1024^2)
options(future.rng.onMisuse="ignore")

# 0.2 color set -----------------------------------------------------------

# Naive effector.memory memory MHCII resi.memory exhausted interferon cycling nk.like mait gdt
CD8.dim.major.col = c('#8dd3c7','#ffffb3','#80b1d3','#8da0cb','#b3de69','#fb8072','#bebada','#d9d9d9','#fc8d62','#ccebc5','#fccde5')
CD8.dim.minor.col = c('#8dd3c7','#ffffb3','#fdb462','#ffed6f','#80b1d3','#8da0cb','#b3de69','#fb8072','#bebada','#d9d9d9','#66c2a5','#fc8d62','#ccebc5','#fccde5','#bc80bd','#e78ac3')

prop.col = '#00868b'

FC.high.col = '#8b1a1a'
FC.low.col = '#4682b4'
FC.line.col = '#bebebe'

CD8.stack.col = CD8.dim.major.col

heatmap.col = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7",
                                 "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)

Tex.dim.col = c('#4dbbd5','#3c5488','#f39b7f','#00a087')

CD8.downsample.major.col = c('#8dd3c7','#ffffb3','#80b1d3','#8da0cb','#b3de69','#fb8072')
CD8.downsample.minor.col = c('#8dd3c7','#ffffb3','#ffed6f','#80b1d3','#8da0cb','#b3de69','#4dbbd5','#3c5488','#f39b7f','#00a087')

# bassez zhang liu caushi rcc bcc scc luoma
study.col = c('#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')
tissue.col = c('#1B9E77','#D95F02')
timepoint.col = c('#7FC97F','#BEAED4')

Tex.bar.col = Tex.dim.col

RvsNR.box.col = c('#bf2c2d','#1b91bb')

# E.pre NE.pre E.post NE.post
Texc02.box.col = c('#de7597','#bf2c2d','#4e8abc','#1b91bb')

# R.post R.pre NR.post NR.pre 
Tex.volcano.point.col = c('#ff3030','#009acd','#cd853f','#2e8b57')
Tex.volcano.margin.col = c('#4d4d4d')

Tex.score.box.col = Tex.dim.col

T.trajectory.col = CD8.dim.minor.col[c(1,2,4,7,8)]
Tex.trajectory.col = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99)

# 0.3 obj load ------------------------------------------------------------
selected.all.T <- readRDS("~/Project/PaperCD8/data/selected.all.T.rds")

CD8 <- readRDS("~/PaperCD8/data/reCD8.finished.rds")

Tprop <- read.csv("~/PaperCD8/data/Tprop.csv")

CD8.Tex.harmony <- readRDS("~/PaperCD8/data/Tex/CD8.Tex.harmony.rds")
Tex.scenic.rss <- readRDS("~/PaperCD8/data/Tex/Tex.scenic.rss.rds")
volcanodf <- readRDS("~/PaperCD8/data/Tex/volcanodf.rds")

CD8.downsample <- readRDS("~/PaperCD8/data/Tex/CD8.downsample.rds")
cytotrace <- readRDS("~/cytotrace.rds")
CD8.genes <- readRDS("~/CD8.genes.rds")
CD8.cds <- readRDS("~/CD8.cds.rds")
CD8.DM <- readRDS("~/CD8.DM.rds")


# 1.1 CD8 dimplot ---------------------------------------------------------


# 1.2 CD8 featureplot -----------------------------------------------------
gene = c('CCR7','IL7R','ZNF683','CXCL13','PDCD1','ISG15','SLC4A10','GZMB','GZMK')

# 1.3 CD8 prop plot -------------------------------------------------------
source("~/PaperCD8/code/visual/plotprop.R")


# 1.4 CD8 FC plot ---------------------------------------------------------
source("~/PaperCD8/code/visual/plotfc.R")


# 1.5 CD8 Bassez dimplot --------------------------------------------------


# 1.6 CD8 stackplot -------------------------------------------------------
source("~/PaperCD8/code/visual/majorstack.R")


# 1.7 CD8 avgheatmap ------------------------------------------------------
gene = c('CCR7','SELL','LEF1','TCF7','IL7R',
         'GZMK','ITM2C','CST7','CD28','KLRG1',
         'NR4A1','EGR1','BAG3','UBE2S',
         'HLA-DRA','HLA-DRB1','HLA-DQA1','HLA-DRB5','HLA-DPA1',
         'XCL1','ZNF683','ANXA1','VIM','LMVA',
         'CXCL13','KRT86','LAYN','GZMB','CTLA4','PDCD1','HAVCR2','TIGIT',
         'IFI6','MX1','ISG15','IFIT3','IFIT1',
         'STMN1','TUBB','LINC00861','S100A4','HMGN2',
         'FGFBP2','KKLRD1','KLRC3','KLRF1','TYROBP',
         'KLRB1','CCL20','CEBPD','CCR6','SLC4A10',
         'TRDC','TRDV2','TRGV9')


# 1.8 CD8 batch effect ----------------------------------------------------
p1 = DimPlot(CD8,group.by = 'Study',cols = study.col,raster = F)+
  labs(title = '',x = 'UMAP1',y = 'UMAP2')+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10))
p2 = DimPlot(CD8,group.by = 'sample.Tissue',cols = tissue.col, raster = F)+
  labs(title = '',x = 'UMAP1',y = '')+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10))
p3 = DimPlot(CD8,group.by = 'sample.timepoint',cols = timepoint.col, raster = F)+
  labs(title = '',x = 'UMAP1',y = '')+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10))

p = p1|p2|p3
ggsave('1.8.CD8.batch.png',p,width = 16,height = 4,dpi = 300)

# 2.1 Tex dimplot ---------------------------------------------------------
p = DimPlot(CD8.Tex.harmony,group.by = 'manual.celltype.Tex',cols = Tex.dim.col,raster = F)+ 
  labs(title = '',x = 'UMAP1',y = 'UMAP2')+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10) )
ggsave('2.1.Tex.dimplot.pdf',p,width = 6,height = 4,dpi = 300, device = cairo_pdf) 

# 2.2 Tex dotplot ---------------------------------------------------------
gene = c('NR4A2','ZNF331','NR4A3','CXCR4','MYADM','CREM', # c04
         'TNFRSF4','KLRB1','SELL','CCR7','LMNA','LTB','IL7R', # c03
         'KLRG1','FCGR3A','CX3CR1','GZMH','NKG7','PLEK','FGFBP2', #c02
         'CXCL13','TNFRSF9','LAG3','IFNG','GZMB','CCL4','CCL4L2','CCL3' #c01 
         )

# 2.3 Tex subcluster heatmap ----------------------------------------------
gene = c('NR4A2','ZNF331','NR4A3','CXCR4','MYADM','CREM', # c04
         'TNFRSF4','KLRB1','SELL','CCR7','LMNA','LTB','IL7R', # c03
         'HAVCR2','CD63','CXCR6','CX3CR1','ZNF683','GNLY','GZMH','CCL5', #c02
         'CXCL13','TNFRSF9','LAG3','IFNG','GZMB','CCL4','CCL4L2','CCL3' #c01
         )


# 2.4 Tex subcluster RvsNR boxplot ----------------------------------------
source("~/PaperCD8/code/visual/clusterbox.R")


# 2.5 Tex.c02.GZMH RvsNR boxplot ------------------------------------------
source("~/PaperCD8/code/visual/texclusterbox.R")


# 2.6 Tex RvsNR forestplot ------------------------------------------------


# 2.7 Tex velcano ---------------------------------------------------------


# 2.8 Tex ssgsea score boxplot & radar plot -------------------------------------------
library(ggradar2)

# 2.9 Tex avgheatmap ------------------------------------------------------
stem.feature = c('CCR7','TCF7','SELL','LEF1','CCR5')
resident.feature = c('NR4A1','NR4A3','CD69','CXCR6','ITGAE')
cyto.feature = c('IFNG','GNLY','GZMB','GZMK','GZMH','GZMA','NKG7','FGFBP2')
exhaust.feature =c('TOX2','SOX4','TIGIT','PDCD1','CTLA4','HAVCR2','LAG3','CXCL13')
costi.feature = c('ICOS','TNFSF14','TNFRSF25','TNFRSF9','CD28','TNFSF4')

gene = c(stem.feature, resident.feature, cyto.feature, exhaust.feature, costi.feature)

# 2.10 trajectory ----------------------------------------------------------

# diffusemap
library(destiny)
meta.all = CD8.downsample@meta.data

tmp <- data.frame(DC1 = eigenvectors(CD8.DM)[, 1],
                  DC2 = eigenvectors(CD8.DM)[, 2],
                  celltype = meta.all$manual.celltype.minor)
tmp$celltype = factor(tmp$celltype,levels = c('Tn.c01.CCR7','Tem.c02.GZMK',
                                              'Tem.c04.GZMH','Trm.c07.ZNF683','Tex.c08.CXCL13'))

ggplot(tmp, aes(x = DC1, y = DC2, colour = celltype)) +
  geom_point() + scale_color_manual(values = T.trajectory.col) +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic()

# monocle2
library(monocle,lib.loc = '/home/shpc_100839/R/x86_64-pc-linux-gnu-library/4.2')
plot_cell_trajectory(CD8.cds, color_by = "manual.celltype.minor",size=1,show_backbone=T)

plotdf = pData(CD8Tex.cds)

library(ggridges)

ggplot(plotdf, aes(x=Pseudotime,y=manual.celltype.major,fill=manual.celltype.major))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
    )

ggplot(plotdf, aes(x=Pseudotime,y=manual.celltype.major,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = Tex.trajectory.col)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )


plotdf2=as.data.frame(t(CD8.downsample.cds@reducedDimS))
colnames(plotdf2)=c("component1","component2")
plotdf2$Pseudotime=CD8.downsample.cds$Pseudotime

ggplot(data = plotdf2, aes(component1,component2,color=Pseudotime))+
  geom_point()+scale_color_gradientn(colours = Tex.trajectory.col)+
  theme_classic()


# cytotrace
cytotrace.score = data.frame(names = names(cytotrace$CytoTRACE),
                             score = cytotrace$CytoTRACE)
plotdf2$names = rownames(plotdf2)
plotdf3 = left_join(plotdf2,cytotrace.score,by= 'names')

ggplot(data = plotdf3, aes(component1,component2,color=score))+
  geom_point()+scale_color_gradientn(colours = Tex.trajectory.col)+
  theme_classic()

plotCytoTRACE(res,emb = CD8.Tex.harmony@reductions$umap_harmony@cell.embeddings)


# 2.11 Tex scenic ---------------------------------------------------------


# 2.12 Tex compass --------------------------------------------------------

beforeres = list.files('~/PaperCD8/data/Tex/compass/before')
beforeres = gsub('.tsv','',beforeres)
beforeres = gsub('.pre','',beforeres)

beforeres = paste('~/PaperCD8/data/Tex/compass/res',beforeres,'reactions.tsv',sep = '/')

res1 = read.delim(beforeres[1])

for (i in 2:length(beforeres)) {
  res = read.delim(beforeres[i])
  res1 = cbind(res1,res)
}

afterres = list.files('~/PaperCD8/data/Tex/compass/after')
afterres = gsub('.tsv','',afterres)
afterres = gsub('.post','',afterres)

afterres = paste('~/PaperCD8/data/Tex/compass/res',afterres,'reactions.tsv',sep = '/')

res2 = read.delim(afterres[1])

for (i in 2:length(afterres)) {
  res = read.delim(afterres[i])
  res2 = cbind(res2,res)
}

# add pseudotime Compass score like SCPA https://jackbibby1.github.io/SCPA/articles/pseudotime.html
