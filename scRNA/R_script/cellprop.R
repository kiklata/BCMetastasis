anno = read.delim('~/Project/MultiOmics/data/skin/res/cellbender_celltype.csv',row.names = 1,sep = ',')

anno = dplyr::filter(anno, cl_major == 'Keratinocyte')

ptexpan = as.data.frame(table(anno$SampleID,anno$cl_subset))
ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
ptexpan$total = apply(ptexpan[,2:ncol(ptexpan)],1,sum)

for (h in 2:(ncol(ptexpan)-1)) {
  ptexpan[,h] = ptexpan[,h]/ptexpan$total
}

pttimepoint = as.data.frame(table(anno$SampleID,anno$SampleType))
pttimepoint = tidyr::spread(pttimepoint,key = 'Var2',value = 'Freq')
pttimepoint$type = if_else(pttimepoint$H!=0,'ARD','Normal')

ptexpan$type = pttimepoint$type
ptexpan$pair = c('a','b','c','a','b','c')


ptexpan$type = factor(ptexpan$type,levels = c('Normal','ARD'))

celltype_list = c("KC-Basal-COL17A1","KC-Basal-ITGA6", "KC-Suprabasal-DSC3" ,    
                  "KC-Suprabasal-KRT6A","KC-Suprabasal-LYPD3","KC-Suprabasal-PLD1",
                  "KC-Spinous-AZGP1","KC-Spinous-SPINK5",           
                  "KC-Proliferating-DIAPH3","KC-Cycling" )

plist = list()

for (i in celltype_list) {

  plist[[i]] = 
    ggplot(data = ptexpan, aes(x = .data[['type']],y = .data[[i]]))+ 
      geom_point()+
      geom_line(aes(group = pair), color = 'lightgrey')+
      labs(title = i, y = "Frequence", x = "")+
      theme_bw()+
      theme(plot.title = element_text(size = 8,hjust = 0.5),
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 4,angle = 0,vjust = 0.5,hjust = 0.7),
            axis.text.y = element_text(size = 4),
            axis.title = element_text(size = 6),
            panel.border = element_blank(), axis.line = element_line(colour = 'black'),
            axis.text = element_text(colour = 'black'))
}
for (i in celltype_list) {
  ggsave(filename = paste0(i,'_line.png'),plist[[i]],width = 1.82,height = 1.36,dpi = 300)
}

library(patchwork)
#p = 
  ((plist[['T cell']]+plist[['Myeloid']]+plist[['Mast cell']]+plist[['Keratinocyte']]+plist[['Eccrine gland']])+plot_layout(ncol = 5))/
  ((plist[['Melanocyte']]+plist[['Endothelial']]+plist[['PVL']]+plist[['Fibroblast']]+plist[['Schwann']])+plot_layout(ncol = 5))

ggsave('porpchange.png',p,width = 8,height = 2.8)


# tcell -------------------------------------------------------------------

tcell_cluster <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/tcell_cluster.tsv", row.names=1)

ptexpan = as.data.frame(table(tcell_cluster$SampleID,tcell_cluster$Tcell_minor))
ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
ptexpan$total = apply(ptexpan[,2:ncol(ptexpan)],1,sum)

for (h in 2:(ncol(ptexpan)-1)) {
  ptexpan[,h] = ptexpan[,h]/ptexpan$total
}

pttimepoint = as.data.frame(table(tcell_cluster$SampleID,tcell_cluster$SampleTimepoint))
pttimepoint = tidyr::spread(pttimepoint,key = 'Var2',value = 'Freq')
pttimepoint$type = if_else(pttimepoint$S1!=0,'pre','post')

ptexpan$type = pttimepoint$type
ptexpan$pair = c('a','b','b')

ptexpan$type = factor(ptexpan$type,levels = c('pre','post'))

celltype_list = c( "CD4-Tn-IL7R","CD4-Treg-FOXP3","CD8-Tem","CD8-Tex-DUSP2","CD8-Tex-ITM2C","CD8-Trm-ZNF683","NKT")
library(ggpubr)
library(patchwork)

plist = list()
for (i in celltype_list) {
  
  plist[[i]] = 
    ggplot(data = ptexpan, aes(x = .data[['type']],y = .data[[i]]))+ 
    #geom_point()+
    geom_boxplot(aes(color = .data[['type']],fill = .data[['type']]))+
    scale_color_manual(values = c('#aab19a','#8b596a'))+
    scale_fill_manual(values = c('#dbe4c7','#be7b92'))+
    stat_compare_means(comparisons = list(c('Normal','ARD')),
                       method = "wilcox.test",
                       label = "p.signif",tip.length = 0,size = 2,lwd = 0.5)+
                       
    labs(caption = i, y = "Frequence", x = "")+
    theme_bw()+theme(plot.caption = element_text(hjust=0.5, size=8,),
                     panel.grid = element_blank(),
                     axis.text.y = element_text(size = 4),
                     axis.title.y = element_text(size = 6),
                     axis.text.x = element_blank(),legend.title = element_blank(),axis.ticks.x = element_blank(),
                     panel.border = element_blank(), axis.line = element_line(colour = 'black'),
                     axis.text = element_text(colour = 'black'),
                     legend.text = element_text(size = 6))
}
p = (plist[['CD4-Tn-IL7R']]+plist[['CD4-Treg-FOXP3']]+plist[['CD8-Tem']]+plist[['CD8-Tex-DUSP2']]+plist[['CD8-Tex-ITM2C']]+plist[['CD8-Trm-ZNF683']])+
  plot_layout(widths = c(1,1,1,1,1,1), guides = 'collect')
for (i in celltype_list) {
  ggsave(filename = paste0(i,'_box.png'),plist[[i]],width = 1.8,height = 1.70,dpi = 300)
}


# myeloid -------------------------------------------------------------------

anno <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/cellranger_filter_anno.tsv", row.names=1)
mye_cluster = dplyr::filter(anno, manual_celltype_annotation == 'Myeloid')

ptexpan = as.data.frame(table(mye_cluster$SampleID,mye_cluster$celltype_minor))
ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
ptexpan$total = apply(ptexpan[,2:ncol(ptexpan)],1,sum)

for (h in 2:(ncol(ptexpan)-1)) {
  ptexpan[,h] = ptexpan[,h]/ptexpan$total
}

pttimepoint = as.data.frame(table(mye_cluster$SampleID,mye_cluster$SampleTimepoint))
pttimepoint = tidyr::spread(pttimepoint,key = 'Var2',value = 'Freq')
pttimepoint$type = if_else(pttimepoint$S1!=0,'pre','post')

ptexpan$type = pttimepoint$type
ptexpan$pair = c('a','b','b')

ptexpan$type = factor(ptexpan$type,levels = c('pre','post'))

celltype_list = c( "Macrophage-IGFBP7","Macrophage-MAML2","Macrophage-TREM2","Monocyte-FCN1","Monocyte-CXCL8","cDC1-CLEC10A","cDC2-CD1C")
library(ggpubr)
library(patchwork)

plist = list()
for (i in celltype_list) {
  
  plist[[i]] = 
    ggplot(data = ptexpan, aes(x = .data[['type']],y = .data[[i]]))+ 
    #geom_point()+
    geom_boxplot(aes(color = .data[['type']],fill = .data[['type']]))+
    scale_color_manual(values = c('#aab19a','#8b596a'))+
    scale_fill_manual(values = c('#dbe4c7','#be7b92'))+
    stat_compare_means(comparisons = list(c('pre','post')),
                       method = "wilcox.test",
                       label = "p.format",tip.length = 0,size = 2,lwd = 0.5)+
    
    labs(caption = i, y = "Frequence", x = "")+
    theme_bw()+theme(plot.caption = element_text(hjust=0.5, size=5,),
                     panel.grid = element_blank(),
                     axis.text.x = element_blank(),legend.title = element_blank(),axis.ticks.x = element_blank(),
                     panel.border = element_blank(), axis.line = element_line(colour = 'black'),
                     axis.text = element_text(colour = 'black'),
                     legend.text = element_text(size = 6))
}
p = (plist[[1]]+plist[[2]]+plist[[3]]+plist[[4]]+plist[[5]]+plist[[6]])+
  plot_layout(widths = c(1,1,1,1,1,1), guides = 'collect')


