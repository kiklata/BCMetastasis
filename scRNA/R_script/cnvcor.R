
cnvIdent = function(cnv_obj, cor.method = 'spearman') {
  
  cnv = as.matrix(cnv_obj@expr.data)
  
  refCells = cnv_obj@observation_grouped_cell_indices %>% unlist() %>% unname()
  
  t_cnv = cnv[,refCells]
  
  #scores=apply(cnv,2,function(x){ sum(x < 0.95 | x > 1.05)/length(x) })
  scores = apply(cnv, 2, function(x){sum(abs(x-1))})
  
  scores %>% hist(., breaks = 100)
  
  top5per_index = scores[scores>=quantile(scores, 0.95)] %>% names
  top_tcnv = cnv[,colnames(cnv) %in% top5per_index] %>% rowMeans(.)
 
  cellcors = suppressWarnings(stats::cor(top_tcnv, cnv, method = cor.method)) %>% unlist(as.data.frame(.)) %>% t()
  #cellcors %>% hist(., breaks = 100)
  cbind(cellcors,scores) %>% plot()
  
  res = cluster::pam(cbind(cellcors,scores)[refCells,],k = 2)
  table(res$clustering)

  res$clustering[res$clustering == 2] %>% names
}

DimPlot(seu, cells.highlight = res$clustering[res$clustering == 2] %>% names)
