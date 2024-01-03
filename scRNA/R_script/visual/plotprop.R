
plot.prop = function(obj,cluster,timepoint,bar.col = 'turquoise4'){
  
  library(ggplot2)

  # celltype.minor by sample before 
  minorsample = as.data.frame(table(obj$sample.ID,obj@meta.data[,cluster]))
  minorsample = tidyr::spread(minorsample,key = 'Var2',value = 'Freq')
  minorsample$totalcd8 = apply(minorsample[,2:ncol(minorsample)],1,sum)
  
  for (i in 2:(ncol(minorsample)-1)) {
    minorsample[,i] = minorsample[,i]/minorsample$totalcd8
  }
  
  minormean = apply(minorsample[,2:(ncol(minorsample)-1)], 2, mean)
  minorsd = apply(minorsample[,2:(ncol(minorsample)-1)], 2, sd)
  
  minordf = data.frame(celltype = names(minormean),mean = minormean,sd = minorsd)
  
  if(cluster == 'manual.celltype.minor') {
    minordf$celltype = factor(minordf$celltype,levels = c('Tn.c01.CCR7','Tem.c02.GZMK','Tem.c03.HMGB2','Tem.c04.GZMH',
                                                        'Tm.c05.NR4A1','T.c06.MHCII','Trm.c07.ZNF683','Tex.c08.CXCL13',
                                                        'T.c09.IFN','T.c10.STMN1','T.c11.ASPN','T.c12.NK-like',
                                                        'MAIT.c13','γδT.c14.GNLY','γδT.c15.TRDV1','γδT.c16.TRDV2'))
  }else if(cluster == 'manual.celltype.major'){
    minordf$celltype = factor(minordf$celltype,levels = c('Naive','Effector memory','Memory','MHC II','Resident memory','Exhausted',
                                                                                'Interferon','Cycling','NK-like','MAIT','γδT'))
    
  }
  
  
  p1 = 
    ggplot(minordf, aes(x = celltype,y = mean))+ geom_bar(stat = "identity",fill = bar.col,color = bar.col)+
    geom_errorbar(aes(ymin = mean,ymax = mean + sd), width = .2, position = position_dodge(0.6),color = bar.col)+
    labs(title = timepoint, y = "Frequence", x = "")+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 60,vjust = 0.5,hjust = 0.7),
                     panel.border = element_blank(), axis.line = element_line(colour = 'black'),
                     axis.text = element_text(colour = 'black'))
  return(p1)
}
