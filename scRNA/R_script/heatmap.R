
## scp
ht1 <- GroupHeatmap(pancreas_sub,
                    features = c(
                      "Sox9", "Anxa2", "Bicc1", # Ductal
                      "Neurog3", "Hes6", # EPs
                      "Fev", "Neurod1", # Pre-endocrine
                      "Rbp4", "Pyy", # Endocrine
                      "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
                    ),
                    group.by = c("CellType", "SubCellType")
)

count =myeloid_count@assays$RNA@counts

mt.name = grep ("^MT-", rownames(count),value = T)
rp.name = grep("^RP[L|S]",rownames(count),value = T)
hb.name = grep("^HB[^(P)]",rownames(count),value = T)

new.count = count[!rownames(count) %in% hb.name,]

new_myeloid = CreateSeuratObject(new.count, meta.data = myeloid_count@meta.data)

new_myeloid = NormalizeData(new_myeloid,assay = 'RNA')
new_myeloid = ScaleData(new_myeloid)

Idents(new_myeloid) = new_myeloid$Myeloid_minor
all_marker = FindAllMarkers(new_myeloid)

tes = all_marker %>% dplyr::group_by(cluster) %>% dplyr::filter(p_val_adj <0.05) %>% 
  top_n(n = 5, wt = avg_log2FC)

removeg = HGNChelper::checkGeneSymbols(x = tes$gene) %>% dplyr::filter(Approved == T) 

tes_l = tes[tes$gene %in% removeg$x,]

vag_exp = AverageExpression(obj, 
                            assays = "RNA", 
                            features = tes_l$gene,
                            group.by = "Myeloid_minor",
                            layer = "data")

celltype = unique(FetchData(obj, vars = c("myeloid_major","Myeloid_minor")))

celltype$myeloid_major = factor(celltype$myeloid_major, 
                                levels = c("Macrophage",'Monocyte','DC'))

celltype = celltype[order(celltype$myeloid_major), ]

celltype$Myeloid_minor = factor(celltype$Myeloid_minor, 
                                levels = unique(celltype$Myeloid_minor))

tes_l$cluster = factor(tes_l$cluster, levels = levels(celltype$Myeloid_minor))
tes_l = tes_l[order(tes_l$cluster),]
tes_l$gene = factor(tes_l$gene, levels = tes_l$gene)

dat = base::apply(vag_exp$RNA, 1, function(x) (x - mean(x)) / sd(x)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene') %>% 
  reshape2::melt()

dat$Gene = factor(dat$Gene, levels = rev(levels(tes_l$gene)))

dat$variable = factor(dat$variable, levels = levels(celltype$Myeloid_minor))

dat = dat[order(dat$Gene), ]

heatmap_color = RColorBrewer::brewer.pal(name = "RdBu",n = 11)

pal = rev(colorRampPalette(heatmap_color)(500))
label1 = levels(celltype$Myeloid_minor)[seq(1, length(levels(celltype$Myeloid_minor)), by = 2)]
label2 = levels(celltype$Myeloid_minor)[seq(2, length(levels(celltype$Myeloid_minor)), by = 2)]


p1 = ggplot(dat, aes(as.numeric(variable), Gene, fill=value))+
  geom_tile() +
  scale_y_discrete(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(1, length(levels(celltype$Myeloid_minor)), by = 2),
                     labels = label1,
                     sec.axis = dup_axis(
                       breaks = seq(2, length(levels(celltype$Myeloid_minor)), by = 2),
                       labels = label2)
  ) +
  scale_fill_gradientn(colors = pal, limits = c(-2.5, 3), name = "Z Score") +
  geom_hline(yintercept = as.numeric(cumsum(rev(table(tes_l$cluster)[-1])) + .5), linetype = 2)+
  geom_vline(xintercept = as.numeric(cumsum(table(celltype$myeloid_major)) + .5), linetype = 2)+
  theme(text = element_text(face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = .5),
        axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = .5)
  )

tes_l$x = rep(c(1,2,3),31)[1:92]
tes_l$y = tes_l$x
initindex = 0
for ( i in names(table(tes_l$cluster))) {
  for (k in 1:nrow(tes_l[tes_l$cluster == i,])) {
    tes_l[tes_l$cluster == i,][k,'y'] = 31-initindex
    if( tes_l[tes_l$cluster == i,][k,'x'] == 3){
      initindex = initindex+1
    }else{NULL}
  }
}


p2 = ggplot(tes_l, aes(x,y,fill = cluster))+
  geom_tile()+
  geom_text(aes(label = gene), family = "serif", fontface = "italic") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_brewer(palette = "Pastel2") +
  theme(text = element_text(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.ticks = element_blank(),
        legend.position="none")+
  scale_x_continuous(expand = c(0,0)) + 
  geom_hline(yintercept = as.numeric(cumsum(rev(table(tes_l$cluster)[-1]/3))) + .5, color = "white")

## markergene number need manula select

library(patchwork)

pic.heatmap = p2 + p1 + plot_layout(ncol = 2, widths  = c(1, 3))

pdf("test.pdf",w=10.5,h=10)
pic.heatmap & theme(plot.margin = margin(0,0,0,0))
dev.off()