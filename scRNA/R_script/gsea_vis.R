df = res$res@result
df$group = if_else(df$NES>0, 'ARD','Normal')

library(GseaVis)

sigs = c(
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_IL2_STAT5_SIGNALING')

for(i in sigs){
  p = gseaNb(object = res$res,geneSetID = i,addPval = T,subPlot = 2)
  ggsave(paste0(i,'_gsea.png'),p, width = 4.98, height = 4.26)
}

