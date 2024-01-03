
stackplot = function(obj, ylabuse = FALSE, mycol= c("#8DD3C7", "#FFFFB3", "#80B1D3", "#8DA0CB", "#B3DE69",
                                                    "#FB8072", "#BEBADA", "#D9D9D9", "#FC8D62", "#CCEBC5",
                                                    "#FCCDE5")){
  p <- ggplot(data = obj) +
    geom_col(aes(x = xname,y = value/CD8,fill = variable),
             position = position_stack(),#stat = 'identity',
             width = 0.9) +
    scale_fill_manual(values = mycol) +
    #scale_y_continuous(labels = scales::percent_format()) +
    theme_gray(base_size = 18) + guides(fill = guide_legend(title = NULL,label.position = 'right',nrow = 2))+
    xlab('')
  if (ylabuse == T) {
    p = p + ylab('Proportion') 
  }
  if (ylabuse == F) {
    p = p + ylab('') + theme(axis.ticks.y = element_blank())
  }
  p = p + 
    theme(axis.text.x = element_text(family = 'arial',size = 10,angle = 45,vjust = 1,hjust = 0.5,colour = 'black'),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          panel.background = element_blank(),
          legend.position = 'bottom',legend.direction = 'horizontal'
          #plot.margin = margin(t = 3,b = 3,r = 1,unit = 'cm')
    )
  return(p)
}

anno = read.delim('~/Project/MultiOmics/data/skin/res/cellbender_celltype.csv',row.names = 1,sep = ',')

#anno = dplyr::filter(anno, cl_major == 'Keratinocyte')

ptexpan = as.data.frame(table(anno$SampleID,anno$cl_minor))
ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
ptexpan$total = apply(ptexpan[,2:ncol(ptexpan)],1,sum)

celltype_list = c("KC-Basal-COL17A1","KC-Basal-ITGA6", "KC-Suprabasal-DSC3" ,    
                  "KC-Suprabasal-KRT6A","KC-Suprabasal-LYPD3","KC-Suprabasal-PLD1",
                  "KC-Spinous-AZGP1","KC-Spinous-SPINK5",           
                  "KC-Proliferating-DIAPH3","KC-Cycling" )

stackbar_cl_order = c('T-CD4','T-CD8','NK','Macrophage','DC','LC','Mast cell',
                     'KC-Basal','KC-Suprabasal','KC-Spinous','KC-Proliferating','KC-Cycling',
                     'Eccrine gland','Melanocyte',
                     'Endo','Endo-Lymph','Endo-Vas','PVL','Fibroblast',
                     'Schwann')


df = ptexpan %>% reshape2::melt(.,c('Var1','total'))
df$Var1 = factor(df$Var1,levels = stackbar_cl_order)

p = ggplot(data = df) +
  geom_col(aes(x = Var1,y = value/total,fill = variable),
           position = position_stack(),
           width = 0.9,color = 'grey25') + 
  coord_flip()+ # cl on, sample off
  #scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 6,angle = 0,hjust = 1,colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 6, colour = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2,'cm'))+ # cl 0.2, sample 0.5
  labs(x = '', y = 'Proportion')+guides(fill = guide_legend(ncol= 2))

ggsave('sample_stackbar.png',p,width = 5.42,height = 2.84,dpi = 300)
ggsave('cl_minor_stackbar.png',p,width = 5.66,height = 1.71,dpi = 300)
