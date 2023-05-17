# demo use SCP package

library(SCP)
library(Seurat)

tcell <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/tcell.harmony.rds")

tcell$scibet.celltype.subset = factor(tcell$scibet.celltype.subset,levels = c('CD4 naive T cell','CD4 memory T cell','CD4 effector memory T cell',
                                                                              'CD4 terminally differentiated effector memory or effector cells T cell',
                                                                              'CD4 helper T cell 17','CD4 follicular helper T cell',
                                                                              "Tfh/Th1 dual-functional T cell",
                                                                              'CD4 interferon-stimulated genes related T cell','CD4 T regulatory cell',
                                                                              'CD4 mix T cell',
                                                                              'CD8 naive T cell',"CD8 memory T cell","CD8 effector memory T cell",
                                                                              "CD8 terminally differentiated effector memory or effector cells T cell",
                                                                              "CD8 interferon-stimulated genes related T cell","CD8 exhausted T cell",
                                                                              "CD8 mucosal-associated invariant T cells","CD8 NK-like T cell",
                                                                              "CD8 tissue-resident memory T cell"))

ClassDimPlot(
  srt = tcell, group.by = c("celltypist.celltype.minor", "scibet.celltype.subset"),
  reduction = "UMAP", theme_use = "theme_blank")
ExpDimPlot(
  srt = tcell, features = c("CCR7", "CXCL13", "PDCD1", "CD274"),
  reduction = "UMAP", theme_use = "theme_blank"
)

source("~/scRNA_BC_metastases/Analysis/scRNA/Function_definition/Geneset_by_Zhang.R", echo=TRUE)

p = ExpDotPlot(
  srt = tcell,
  features = c(
    chemokines[1:3],
    costimulatory.molecules[1:3], 
    effector.memory.molecules[1:3], 
    HLA.genes[1:3], 
    inhibitory.receptors[1:3], 
    transcription.factors[1:3]
  ),
  cell_split_by =  "scibet.celltype.marker"
)

p+ggplot2::coord_flip()

tcell <- RunDEtest(srt = tcell, group_by = "scibet.celltype.marker", only.pos = FALSE, fc.threshold = 1)
VolcanoPlot(srt = tcell, group_by = "scibet.celltype.marker")
#DEGs <- tcell@tools$DEtest_CellType$AllMarkers_wilcox
#DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]

manual_gene = data.frame(gene = c(chemokines,
                chemokine.receptors,
  costimulatory.molecules, 
  effector.memory.molecules, 
  HLA.genes, 
  inhibitory.receptors, 
  transcription.factors,
  integrins),genegroup = c('chemokines','chemokine.receptors','costimulatory.molecules','effector.memory.molecules','HLA.genes',
                                'inhibitory.receptors','transcription.factors','integrins'))

ht <- ExpHeatmap(
  srt = tcell, features = manual_gene$gene, feature_split = manual_gene$genegroup, cell_split_by = "scibet.celltype.marker",
  species = "Homo_sapiens", anno_terms = F, anno_keys = F, anno_features = T,
  row_title_size = 0, height = 5, width = 7
)
print(ht$plot)
