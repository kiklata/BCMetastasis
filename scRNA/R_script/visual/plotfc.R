
plot.fc = function(obj,cluster,norm.col = 'steelblue',sign.col = 'firebrick4',mytitle = 'Title'){
  
  library(ggplot2)

  ptexpan = as.data.frame(table(obj$sample.ID,obj@meta.data[,cluster]))
  ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
  rownames(ptexpan) = ptexpan$Var1
  ptexpan$Var1 = NULL
  ptexpan$total = apply(ptexpan,1,sum)
  
  for (h in 1:(ncol(ptexpan)-1)) {
    ptexpan[,h] = ptexpan[,h]/ptexpan$total
  }
  
  ptefficacy = as.data.frame(table(obj$sample.ID,obj$treatment.efficacy))
  ptefficacy = tidyr::spread(ptefficacy,key = 'Var2',value = 'Freq')
  rownames(ptefficacy) = ptefficacy$Var1
  ptefficacy$Var1 = NULL
  ptefficacy$type = if_else(ptefficacy$NR==0,'R','NR')
  
  ptexpan$type = ptefficacy$type
  
  Response = dplyr::filter(ptexpan,type == 'R')
  NonR = dplyr::filter(ptexpan,type == 'NR')
  
  celltype.r = apply(Response[,1:(ncol(Response)-2)],2,mean)
  celltype.nr = apply(NonR[,1:(ncol(Response)-2)],2,mean)
  
  celltype.fc = celltype.r/celltype.nr
  
  for (i in 1:length(celltype.fc)) {
    celltype.fc[i] = dplyr::if_else(celltype.fc[i]<1,-(1/celltype.fc[i]),celltype.fc[i])
  }
  
  celltype.fc.d = as.data.frame(celltype.fc)
  celltype.fc.d$celltype = rownames(celltype.fc.d)
  pos.threhold = head(celltype.fc.d$celltype.fc[order(celltype.fc.d$celltype.fc,decreasing = T)],2)[2]
  neg.threhold = head(rev(celltype.fc.d$celltype.fc[order(celltype.fc.d$celltype.fc,decreasing = T)]),2)[2]
  celltype.fc.d$type = if_else((celltype.fc.d$celltype.fc < pos.threhold)&
                                 (celltype.fc.d$celltype.fc > neg.threhold),'norm','sign')
 
  
  if(cluster == 'manual.celltype.minor') {
    celltype.fc.d$celltype = factor(celltype.fc.d$celltype,levels = celltype.fc.d$celltype[order(celltype.fc.d$celltype.fc,decreasing = T)])
    #celltype.fc.d$celltype = factor(celltype.fc.d$celltype,levels = c('Tn.c01.CCR7','Tem.c02.GZMK','Tem.c03.HMGB2','Tem.c04.GZMH',
    #                                                'Tm.c05.NR4A1','T.c06.MHCII','Trm.c07.ZNF683','Tex.c08.CXCL13',
    #                                                'T.c09.IFN','T.c10.STMN1','T.c11.ASPN','T.c12.NK-like',
    #                                                'MAIT.c13','γδT.c14.GNLY','γδT.c15.TRDV1','γδT.c16.TRDV2'))
  }else if(cluster == 'manual.celltype.major'){
    celltype.fc.d$celltype = factor(celltype.fc.d$celltype,levels = celltype.fc.d$celltype[order(celltype.fc.d$celltype.fc,decreasing = T)])
    #celltype.fc.d$celltype = factor(celltype.fc.d$celltype,levels = c('Naive','Effector memory','Memory','MHC II','Resident memory','Exhausted',
    #                                                'Interferon','Cycling','NK-like','MAIT','γδT'))
  }
  
  #celltype.fc.d = celltype.fc.d[order(celltype.fc.d$celltype.fc,decreasing = T),]
  #celltype.fc.d$celltype.fc = factor(celltype.fc.d$celltype.fc,levels = celltype.fc.d$celltype.fc)
  #norm = filter(celltype.fc.d,type=='norm')
  #sign = filter(celltype.fc.d,type=='sign')

  # or use ggpubr::ggbarplot 
  p = ggplot(celltype.fc.d,aes(x = celltype,y = celltype.fc))+
    geom_bar(aes(fill = type),stat = "identity",show.legend = F)+
    scale_fill_manual(values = c(norm.col,sign.col))+
    geom_hline(yintercept = c(1,-1),linetype = 4,col = 'grey')+
    labs(title = mytitle, y = "Fold Change (R/NR)", x = "")+
    theme_bw()+theme(plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 60,vjust = 0.5,hjust = 0.7),
                     panel.border = element_blank(), axis.line = element_line(colour = 'black'),
                     axis.text = element_text(colour = 'black')
                     )
  return(p)
}
