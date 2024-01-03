convertSeu5Format = function(seu,savepaths){
  seu[["RNA3"]] <- as(object = seu[["RNA"]], Class = "Assay")
  DefaultAssay(seu) <- "RNA3"
  seu[["RNA"]] <- NULL
  seu <- RenameAssays(object = seu, RNA3 = 'RNA')
  sceasy::convertFormat(seu, from="seurat", to="anndata",outFile = savepaths, drop_single_values = F, main_layer = 'counts')
}

