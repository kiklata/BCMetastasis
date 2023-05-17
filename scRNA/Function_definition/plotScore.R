
plotScore = function(seu = seu.brain,savefile = 'file.pdf'){
  
  myggplot = theme_light(base_size = 15)+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.line.x = element_blank(), axis.line.y = element_blank(),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_blank())
  
  p1 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P1)) +
    geom_point(size=.5) +labs(title = "P1")+
    scale_color_gradient(low = 'grey',high = 'darkred')+myggplot
  
  p2 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P2)) +
    geom_point(size=.5) +labs(title = "P2")+
    scale_color_gradient(low = 'grey',high = 'darkorchid')+myggplot
  
  p3 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P3)) +
    geom_point(size=.5) +labs(title = "P3")+
    scale_color_gradient(low = 'grey',high = 'darkgreen')+myggplot
  
  p4 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P4)) +
    geom_point(size=.5) +labs(title = "P4")+
    scale_color_gradient(low = 'grey',high = 'darkgoldenrod')+myggplot
  
  p5 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P5)) +
    geom_point(size=.5) +labs(title = "P5")+
    scale_color_gradient(low = 'grey',high = 'darkcyan')+myggplot
  
  p6 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P6)) +
    geom_point(size=.5) +labs(title = "P6")+
    scale_color_gradient(low = 'grey',high = 'deepskyblue')+myggplot
  
  p7 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P7)) +
    geom_point(size=.5) +labs(title = "P7")+
    scale_color_gradient(low = 'grey',high = 'deeppink')+myggplot
  
  p8 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P8)) +
    geom_point(size=.5) +labs(title = "P8")+
    scale_color_gradient(low = 'grey',high = 'firebrick')+myggplot
  
  p9 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
              aes(UMAP_1, UMAP_2, color= P9)) +
    geom_point(size=.5) +labs(title = "P9")+
    scale_color_gradient(low = 'grey',high = 'darkorange')+myggplot
  
  p10 = ggplot(data.frame(seu@meta.data, seu@reductions$umap@cell.embeddings), 
               aes(UMAP_1, UMAP_2, color= P10)) +
    geom_point(size=.5) +labs(title = "P10")+
    scale_color_gradient(low = 'grey',high = 'darkslategray')+myggplot
  
  p = (p1|p2|p3|p4|p5)/(p6|p7|p8|p9|p10)
  
  ggsave(savefile,p,width = 18,height = 5,dpi = 72)
}
