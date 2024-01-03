library(slingshot)

embed_seu = readRDS("~/Project/MultiOmics/data/skin/res/cellbender_keratinocyte_subset.rds")
seu = readRDS("~/Project/MultiOmics/data/skin/res/cellbender_keratinocyte_count.rds")

seu@reductions = embed_seu@reductions

sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$cl_subset)

index_highlight = subset(seu, cl_subset %in% subset_high) %>% colnames()

unique(ggplot_build(p)$data[[1]]["colour"])





p = DimPlot(seu, group.by = 'cl_subset')+scale_color_manual(values = c('KC-Basal-COL17A1' = '#1f77b4','KC-Basal-ITGA6'= 'grey',
                                                                       'KC-Cycling' = 'grey',
                                                                       'KC-Proliferating-DIAPH3' = '#d62728',
                                                                       'KC-Suprabasal-KRT6A'='grey','KC-Suprabasal-LYPD3'= 'grey',
                                                                       'KC-Suprabasal-DSC3' = 'grey','KC-Suprabasal-PLD1' = '#aec7e8',
                                                                       'KC-Spinous-AZGP1'= '#aa40fc','KC-Spinous-SPINK5'= '#8c564b'))

p

curves = slingCurves(sds, as.df = T)
p2 = p + geom_path(data = curves %>% dplyr::filter(.,Lineage %in% c(4)),aes(x = UMAP_1,y = UMAP_2,group = Lineage),color = 'grey25')+NoLegend()
ggsave('/home/zhepan/Project/MultiOmics/code/figures/kera_slingshot.png',width = 3.82,height = 3.46,dpi = 300)

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

lin1 <- getLineages(rd, cl, start.clus = '1')
crv1 <- getCurves(lin1)
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')


library(Polychrome)
library(ggbeeswarm)
library(ggthemes)

# this define the cluster color. You can change it with different color scheme.
my_color <- createPalette(length(levels(sce$ident)), c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(sce$ident))

slingshot_df <- data.frame(colData(sce))

# re-order y-axis for better figure: This should be tailored with your own cluster names
# slingshot_df$ident = factor(slingshot_df$ident, levels=c(4,2,1,0,3,5,6))

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = ident, 
                         colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)


# SCP ---------------------------------------------------------------------
library(SCP)

embed_seu = readRDS("~/Project/MultiOmics/data/skin/res/cellbender_keratinocyte_subset.rds")
seu = readRDS("~/Project/MultiOmics/data/skin/res/cellbender_keratinocyte_count.rds")

seu@reductions = embed_seu@reductions

seu_ss <- RunSlingshot(srt = seu, group.by = "cl_subset", reduction = "UMAP")

FeatureDimPlot(seu_ss, features = paste0("Lineage", 4), reduction = "UMAP", theme_use = "theme_blank")

CellDimPlot(seu_ss, group.by = "cl_subset", reduction = "UMAP",
            lineages = c('Lineage4'),lineages_palcolor = 'black')+ 
  labs(title = '',caption = '')+theme(legend.position = 'bottom')



pt_seu = seu_ss %>% subset(.,cl_subset %in% c('KC-Basal-COL17A1','KC-Proliferating-DIAPH3','KC-Suprabasal-PLD1','KC-Spinous-AZGP1','KC-Spinous-SPINK5'))
pt_seu$cl_subset = droplevels(pt_seu$cl_subset)

pt_seu <- RunDynamicFeatures(srt = pt_seu, lineages = c("Lineage4"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = pt_seu,
  lineages = c("Lineage4"),feature_split_by = 'Lineage4',
  use_fitted = TRUE, n_split = 1, cell_annotation = 'cl_subset', 
  heatmap_palette = "viridis",
  pseudotime_label = 25, 
  height = 5, width = 2
)
ggsave('/home/zhepan/Project/MultiOmics/code/figures/kera_ht.png',ht$plot,width = 5.1,height = 6.47,dpi = 300)

pd = DynamicPlot(
  srt = pt_seu,
  lineages = c("Lineage4"), group.by = "cl_subset",
  features = c("KRT1",'KRT14','IL18'),
  compare_lineages = TRUE, compare_features = FALSE
)
ggsave('/home/zhepan/Project/MultiOmics/code/figures/kera_dp.png',pd,width = 10.28,height = 3.02,dpi = 300)
