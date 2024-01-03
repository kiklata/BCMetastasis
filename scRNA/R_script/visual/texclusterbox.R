
clusteronetexbox = function(obj,tissues,studys,denominator,removezero,ylabuse = FALSE,
                            titles,cluster,mycolor = c("#DE7597", "#BF2C2D", "#4E8ABC","#1B91BB")){
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(ggpubr)
  library(tidyr)
  
  df = filter(obj, tissue == tissues)
  tmp = df[,c("Sample", "timepoint", "study", "cancertype", "response",'ALLT','CD8','Exhausted','Tex.c01.CCL4','Tex.c02.GZMH','Tex.c03.IL7R','Tex.c04.CREM')]
  
  long.pt.pre = 
    gather(tmp,key=textype,value = count,-c("Sample", "timepoint", "study", "cancertype", "response",'ALLT','CD8','Exhausted'))
  
  long.pt.pre$prop = long.pt.pre$count/long.pt.pre[,denominator]
  
  if(removezero == T){
    long.pt.pre = filter(long.pt.pre, prop>0)}
  
  plot.df = filter(long.pt.pre,study == studys)

  plot.df$xname = paste0(plot.df$timepoint,plot.df$response)
  
  plot.df$xname = if_else(plot.df$xname == 'afterE','E (Post)',
                  if_else(plot.df$xname == 'afterNE','NE (Post)',
                  if_else(plot.df$xname == 'beforeE','E (Pre)',
                  if_else(plot.df$xname == 'beforeNE','NE (Pre)',
                  if_else(plot.df$xname == 'beforeNR','NR (Pre)',
                  if_else(plot.df$xname == 'beforeR','R (Pre)',
                  if_else(plot.df$xname == 'afterNR','NR (Post)',
                  if_else(plot.df$xname == 'afterR','R (Post)','a'))))))))
  if(studys == 'bassez'){
    plot.df$response = factor(plot.df$response,levels = c('E','NE'))
    plot.df$xname = factor(plot.df$xname,levels = c('E (Pre)','E (Post)','NE (Pre)','NE (Post)'))
    mycomprassion = list(c('E (Post)','E (Pre)'),c('NE (Pre)','E (Pre)'),c('NE (Post)','E (Post)'))
    
  }else {
    plot.df$response = factor(plot.df$response,levels = c('R','NR'))
    plot.df$xname = factor(plot.df$xname,levels = c('R (Pre)','R (Post)','NR (Pre)','NR (Post)'))
    if(studys == 'liu'){
      mycomprassion = list(c('R (Post)','R (Pre)'),c('NR (Post)','R (Post)'))
    }else if(studys %in% c('caushi','kevin')){
      mycomprassion = list(c('NR (Post)','R (Post)'))}
    else{
    mycomprassion = list(c('R (Post)','R (Pre)'),c('NR (Pre)','R (Pre)'),c('NR (Post)','R (Post)'))}
  }
  
  plot.df = filter(plot.df,textype == cluster)
  
  ytitle =expression(atop('Frequency of Tex.c02.GZMH',
                          'in all CD8 T cells'))

  p = ggboxplot(plot.df, x = 'xname', y = 'prop',
            color = 'xname',add = 'jitter',outlier.shape = 19)+scale_color_manual(values = mycolor)+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.8))+
    stat_compare_means(comparisons = mycomprassion,
                       method = "wilcox.test",
                       label = "p.format",tip.length = 0,size = 2,lwd = 0.5
                       #symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1,100),
                      #                  symbols = c("***", "**", "*", "",""))
                       
                      )+
    guides(color = 'none')+
  labs(title = titles,x = '', y = '')+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black'),
        legend.position = 'bottom')

  if (ylabuse == T) {
    p = p + ylab(label = ytitle) 
  }
  if (ylabuse == F) {
    p = p + ylab('')
  }
  
  return(p)
}
