run.gsea = function(marker,clusters,gmt,method){
  
  meta = marker %>% dplyr::filter(.,cluster == clusters)
  meta = meta %>% arrange(desc(avg_log2FC))
  meta = meta[,c('gene','avg_log2FC')]
  rank = tibble::deframe(meta)    
  
  for (i in 1:length(rank)) {
    if(rank[i] == Inf){
      rank[i] = 999
    }else if(rank[i] == -Inf){
      rank[i] = -999
    }
    
  }
  
  if(method == 'gsea'){
    library(clusterProfiler)
    
    gmts = read.gmt(gmt)
    res = GSEA(rank,TERM2GENE = gmts,pvalueCutoff = 100)
    
  }else if(method == 'fgsea'){
    library(fgsea)
    library(cogena)
    
    gmts = cogena::gmt2list(gmt)
    res = fgsea::fgsea(gmts,rank)
  }
  return(list(rank = rank, res = res))
}
