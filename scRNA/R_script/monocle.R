# monocle
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)

source('~/scRNA_BC_metastases/Analysis/scRNA/Function_definition/ProgramScoring.R')


# select -------------------------------------------------------------------

# brain

brain.select = subset(brain,sample == 'GSE202501')   
brain.select@active.assay = 'RNA'

brain.select.cds = as.cell_data_set(brain.select)

brain.select.cds = preprocess_cds(brain.select.cds, num_dim = 50)
brain.select.cds = reduce_dimension(brain.select.cds)
brain.select.cds = cluster_cells(brain.select.cds,reduction_method ="UMAP")
brain.select.cds = align_cds(brain.select.cds,alignment_group = 'sample')

brain.select.cds <- learn_graph(brain.select.cds)
brain.select.cds <- order_cells(brain.select.cds, reduction_method = "UMAP",root_cells = colnames(subset(brain.select,type == 'P4')))


brain.select.cds@colData@listData[["type"]] = factor(brain.select.cds@colData@listData[["type"]],levels = c('P1','P2','P3','P4','P5',
                                                                                                            'P6','P7','P8','P9','P10'))
plot_cells(brain.select.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           reduction_method = 'UMAP',
           show_trajectory_graph = T,
           cell_size = 0.5)

Track_genes = graph_test(brain.select.cds,neighbor_graph="principal_graph", cores=6)

pr_deg_ids <- row.names(subset(Track_genes, q_value < 0.05))

rowData(brain.select.cds)$gene_name <- rownames(brain.select.cds)
rowData(brain.select.cds)$gene_short_name <- rowData(brain.select.cds)$gene_name

p = plot_cells(brain.select.cds,genes= Track_genes_sig, show_trajectory_graph=FALSE,
               label_cell_groups=FALSE,label_leaves=FALSE)

p$facet$params$ncol <- 5
ggsave("Genes_Featureplot.pdf",plot = p, width = 20, height = 8)


# tcell -------------------------------------

tcell.anno <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/tcell.anno.rds")
tcell <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/immune/tcell.rds")

tcell = AddMetaData(tcell,tcell.anno)

# after harmony
CD4@active.assay = 'RNA'

CD4.cds = as.cell_data_set(CD4)

#CD4.cds = preprocess_cds(CD4.cds, num_dim = 50)
#CD4.cds = reduce_dimension(CD4.cds)
CD4.cds = cluster_cells(CD4.cds,reduction_method ="UMAP")

CD4.cds <- learn_graph(CD4.cds)
CD4.cds <- order_cells(CD4.cds, reduction_method = "UMAP",root_cells = colnames(subset(CD4,scibet.celltype.subset == 'CD4 naive T cell')))

CD8.cds = as.cell_data_set(CD8)
CD8.cds = cluster_cells(CD8.cds,reduction_method ="UMAP")
CD8.cds <- learn_graph(CD8.cds)
CD8.cds <- order_cells(CD8.cds, reduction_method = "UMAP",root_cells = colnames(subset(CD8,scibet.celltype.subset == 'CD8 naive T cell')))


p = plot_cells(CD4.cds,
           color_cells_by = "scibet.celltype.marker",
           label_cell_groups = FALSE,
           label_groups_by_cluster = FALSE,
           label_roots = FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.5,
           trajectory_graph_color = 'grey28',
           rasterize = TRUE)
ggsave('cd4.pdf',p,height = 10,width = 18,dpi = 72)


p1 = plot_cells(CD8.cds,
               color_cells_by = "scibet.celltype.marker",
               label_cell_groups = FALSE,
               label_groups_by_cluster = FALSE,
               label_roots = FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE,
               cell_size = 0.5,
               trajectory_graph_color = 'grey28',
               rasterize = TRUE)
ggsave('cd8.pdf',p1,height = 10,width = 18,dpi = 72)
