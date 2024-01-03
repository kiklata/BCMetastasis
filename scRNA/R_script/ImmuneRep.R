library(scRepertoire)
setwd("~/Project/MultiOmics/data/scImmune")

kit = c('TCR')
sample = c('P1013S2','P1015S2','P1018S1')
path = c('contig/filtered_contig_annotations.csv')

anno.matrix <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/anno.tsv",sep = '\t',row.names = 1)

# obj create --------------------------------------------------------------
tcr.list = list()

for (i in sample) {
  tcr.list[[i]] = read.csv(file.path(kit,i,path))
}

contiglist = loadContigs(tcr.list, format = "10X")

Tcell_index = dplyr::filter(anno.matrix, manual_celltype_annotation == 'T cell') %>% rownames()

for (i in 1:length(sample)) {
  contiglist[[i]]$barcode = gsub('-1',replacement = paste0('-',c(i-1)),contiglist[[i]]$barcode)
  contiglist[[i]]$contig_id = gsub('-1',replacement = paste0('-',c(i-1)),contiglist[[i]]$contig_id)
  
  contiglist[[i]] = contiglist[[i]] %>% dplyr::filter(., barcode %in% Tcell_index)
}


combined.TCR <- combineTCR(contiglist, 
                           samples = c("P1013S2","P1015S2","P1018S1"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)
combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "timepoint", 
                            variables = c('S2','S2','S1'))

#combined.BCR <- combineBCR(bcr[bcr$barcode %in% Bcell_index,], 
#                           samples = "P1013S2", 
#                           threshold = 0.85)

# basic visual ----------------------------------------------------------
clonalQuant(combined.TCR, cloneCall="strict",chain = "both",scale = TRUE)+ 
  scale_fill_manual(values = c('#1d76b4','#ff983c','#269e67'))
clonalQuant(combined.TCR, cloneCall = "gene", group.by = "timepoint", scale = TRUE)+
  scale_fill_manual(values = c('#dbe4c7','#be7b92'))

clonalAbundance(combined.TCR, cloneCall = "gene", scale = FALSE)+
  scale_color_manual(values = c('#1d76b4','#ff983c','#269e67'))
clonalAbundance(combined.TCR, group.by = "sample",  scale = TRUE)

clonalCompare(combined.TCR, top.clones = 10, 
              samples = c("P1013S2","P1018S1"), 
              relabel.clones = TRUE,
              cloneCall="aa", graph = "alluvial")

clonalScatter(combined.TCR, 
              cloneCall ="gene", 
              x.axis = "P1013S2", 
              y.axis = "P1018S1",
              dot.size = "total",
              graph = "proportion")

clonalHomeostasis(combined.TCR, cloneCall = "gene",cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded =1))


df.genes <- percentGenes(combined.TCR, 
                         chain = "TRB", 
                         gene = "Vgene", 
                         exportTable = TRUE)

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic() 


p = clonalDiversity(combined.TCR, metrics = c("shannon", "gini.simpson"),
                cloneCall = "gene", 
                n.boots = 100)
p$layers[[1]] = NULL

p + scale_fill_manual(values = c('#1d76b4','#ff983c','#269e67'))+
  labs(caption = i, y = "Frequence", x = "")+
  theme_bw()+theme(plot.caption = element_text(hjust=0.5, size=8,),
                   panel.grid = element_blank(),
                   axis.text.x = element_blank(),legend.title = element_blank(),axis.ticks.x = element_blank(),
                   panel.border = element_blank(), axis.line = element_line(colour = 'black'),
                   axis.text = element_text(colour = 'black'),
                   legend.text = element_text(size = 6))+
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    #strip.text.x = element_blank()
  )
  

clonalRarefaction(combined.TCR, 
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)

clonalSizeDistribution(combined.TCR, cloneCall = "gene", method= "ward.D2")

tcell_cluster <- readRDS("~/Project/MultiOmics/data/snRNA/Object/summary/Immune/tcell_cluster.rds")

colnames(tcell_cluster) = paste0(tcell_cluster$SampleID,'_',colnames(tcell_cluster))

sce <- combineExpression(combined.TCR, 
                         tcell_cluster, 
                         cloneCall="gene", 
                         group.by = "sample", 
                         proportion = TRUE)
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

DimPlot(sce, group.by  = "cloneSize")+
  scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))

seu <- combineExpression(combined.TCR, 
                                   tcell_cluster, 
                                   cloneCall="gene", 
                                   group.by = "sample", 
                                   proportion = FALSE, 
                                   cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

DimPlot(seu, group.by = "cloneSize") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))+theme(legend.text = element_text(size = 8))

alluvialClones(seu, 
               cloneCall = "gene", 
               y.axes = c('Tcell_minor','SampleTimepoint'),color = 'Tcell_minor')
