
source("~/Project/MultiOmics/code/func/rungsea.R")

# enrichment res
plotdf = df %>% dplyr::filter(., NES>1.72 | NES< (-0.5) )

p = 
  ggplot(plotdf)+
  geom_col(aes(reorder(Description,NES),y = NES,fill=group))+
  scale_fill_manual(values = c("#8D4873","#1084A4"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5))+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=nrow(plotdf)+1))+
  geom_text(data = plotdf[plotdf$NES> 0.5,],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=2)+
  geom_text(data = plotdf[plotdf$NES< -0.5,],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=2)+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_continuous(breaks = c(-1, -0.5, 0, .5, 1))+
  labs(title = '',x='', y='Normalized Enrichment Score')

ggsave('h_gsea_bar.png',width = 6.04,height = 3.74,dpi = 300)
