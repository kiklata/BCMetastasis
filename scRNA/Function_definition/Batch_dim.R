# batch dim ---------------------------------------------------------------
library(Seurat)
library(SCP)
library(patchwork)

dim_batch_plot = function(obj){
  
  
  p1 = DimPlot(obj, group.by = 'study',reduction = "umap")
  p2 = DimPlot(obj, split.by = 'study',reduction = "umap")+NoLegend()
  p3 = DimPlot(obj,reduction = 'umap',group.by = 'site')
  p4 = DimPlot(obj,reduction = 'umap',split.by = 'site')+NoLegend()
  
  p5 = p1 + p2 + plot_layout(widths = c(1, 10))
  p6 = p3+p4 + plot_layout(widths = c(1, 4))
  p7 = p5/p6+plot_layout(heights = c(1,2))
  
  return(p7)
}

