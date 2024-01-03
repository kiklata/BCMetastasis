# boxsubcluster------------------

clusterbox = function(obj,tissues,sampletime,tcr,denominator,removezero,titles,
                      mycolor = c( "#BF2C2D", "#1B91BB")){
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  
  df = filter(obj, tissue == tissues)
  tmp = df[,c("Sample", "timepoint", "study", "cancertype", "response",'CD8','Exhausted','Tex.c01.CCL4','Tex.c02.GZMH','Tex.c03.IL7R','Tex.c04.CREM')]
  
  long.pt.pre = 
    tidyr::gather(tmp,key=textype,value = count,-c("Sample", "timepoint", "study", "cancertype", "response",'CD8','Exhausted'))
  long.pt.pre$prop = long.pt.pre$count/long.pt.pre[,denominator]
  
  if(removezero == T){
    long.pt.pre = filter(long.pt.pre, prop>0)}
  
  if(tcr == T){
    plot.df = filter(long.pt.pre,study == 'bassez')
  }else {
    plot.df = filter(long.pt.pre,study != 'bassez')
  }
  
  plot.df = filter(plot.df,timepoint == sampletime)
  
  if(tcr == T){
    plot.df$response = factor(plot.df$response,levels = c('E','NE'))
  }else {
    plot.df$response = factor(plot.df$response,levels = c('R','NR'))
  }
  
  ytitle =expression(atop('Frequency of exhausted cell subclusters',
                          'in all CD8 T cells'))
  
  
  p = ggboxplot(plot.df, x = 'textype', y = 'prop',
            color = 'response',outlier.shape = 19)+scale_color_manual(values = mycolor)+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
    stat_compare_means(aes(group = response),
                       method = "wilcox.test",
                       label = "p.format",
                       #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                      #                  symbols = c("***", "**", "*", ""))
                      )+
    labs(title = titles,x = '', y = ytitle)+
    theme(plot.title = element_text(hjust = 0.5),
          
          legend.position = 'bottom',)
  
  return(p)
}

