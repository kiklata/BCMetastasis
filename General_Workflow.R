# import pkg ---------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

#try to use 'singlet' for NMF
#and 'SCP' for visual demo

setwd('scRNA_BC_metastases')

source('script/function.R')

# filtering ---------------------------------------------------------------

# individual data filter scDbifinder doublet

all = readRDS('merge/merge_unfilter.rds')

print('loaded')

all$percent_mt = PercentageFeatureSet(all,pattern = '^MT-')

VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
        ncol = 3,group.by = 'study',raster = F)

FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(all, feature1 = "percent_mt", feature2 = "nFeature_RNA")


all_filter = subset(all, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 
                    & nCount_RNA < 80000 & percent_mt <20)

VlnPlot(all_filter, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
        ncol = 3,group.by = 'study',raster = F)

FeatureScatter(all_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(all_filter, feature1 = "percent_mt", feature2 = "nFeature_RNA")

# PCA ---------------------------------------------------------------------

all_filter = all_filter %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() 
saveRDS(all_filter,file = 'merge/merge_filter_pca.rds')

print('pca saved')

# RPCA ---------------------------------------------------------------------

all_filter_list <- SplitObject(all_filter, split.by = 'sample')

all_filter_list <- lapply(X = all_filter_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = all_filter_list, nfeatures = 3000)
all_filter_list <- PrepSCTIntegration(object.list = all_filter_list, anchor.features = features)
all_filter_list <- lapply(X = all_filter_list, FUN = RunPCA, features = features)

print('pca')

all.anchors <- FindIntegrationAnchors(object.list = all_filter_list,normalization.method = "SCT",
                                      anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20)
print('anchor found')

all_filter_merge <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT", dims = 1:50)

DefaultAssay(all_filter_merge) <- "integrated"

all_filter_merge = all_filter_merge %>% 
RunPCA(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = 0.8)

print('rpca integrated')

saveRDS(all_filter_merge,file = 'merge/merge_rpca.rds')

print('rpca saved')



# harmony -----------------------------------------------------------------
library(harmony)
library(glmGamPoi)

all_filter_list <- SplitObject(all_filter, split.by = 'sample')
all_filter_list <- lapply(X = all_filter_list, FUN = SCTransform, method = "glmGamPoi")
var_features_ind_sct <- SelectIntegrationFeatures(all_filter_list, nfeatures = 3000)

all_filter = merge(x = all_filter_list[[1]], y = all_filter_list[2:length(all_filter_list)])
VariableFeatures(all_filter) <- var_features_ind_sct

saveRDS(all_filter,'merge/merge_filter_sct.rds')

#all_filter = readRDS('merge/merge_filter_sct.rds')
all_filter <- RunPCA(all_filter)

all_filter = RunHarmony(all_filter,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                        kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

all_filter = all_filter %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

saveRDS(all_filter,file = 'merge/merge_harmony.rds')
print('harmony saved')

rm(list = ls())


# LISI --------------------------------------------------------------------

# evaluate data integration

library(lisi)
library(ggthemes)
library(ggpubr)

ifnb.data.raw = readRDS('merge/merge_filter_pca.rds')
ifnb.data.harmony = readRDS("merge/merge_harmony.rds")
ifnb.data.RPCA = readRDS("merge/merge_rpca.rds")

# raw
ifnb.data.raw = ifnb.data.raw %>% FindNeighbors(dims = 1:30)  %>% FindClusters()


raw.score =  lisi::compute_lisi(X = Embeddings(ifnb.data.raw,reduction = "pca"),
                                meta_data = ifnb.data.raw@meta.data,
                                c('study','sample','RNA_snn_res.0.8')) %>% 
  dplyr::mutate(type = "Raw")

saveRDS(raw.score,file = 'merge/LISI/raw.score.rds')

#SCT


pca.score =  lisi::compute_lisi(X = Embeddings(ifnb.data.harmony,reduction = "pca"),
                                meta_data = ifnb.data.harmony@meta.data,
                                c('study','sample','SCT_snn_res.0.8')) %>% 
  dplyr::mutate(type = "SCT")

#RPCA

RPCA.score = lisi::compute_lisi(X = Embeddings(ifnb.data.RPCA,reduction = "pca"),
                                meta_data = ifnb.data.RPCA@meta.data,
                                c('study','sample','integrated_snn_res.0.8')) %>% 
  dplyr::mutate(type = "RPCA")

#harmony


harmony.score = lisi::compute_lisi(X = Embeddings(ifnb.data.harmony,reduction = "harmony"),
                                   meta_data = ifnb.data.harmony@meta.data,
                                   c('study','sample','SCT_snn_res.0.8')) %>% 
  dplyr::mutate(type = "Harmony")

colnames(raw.score)[3] = 'cluster'
colnames(sct.score)[3] = 'cluster'
colnames(rpca.score)[3] = 'cluster'
colnames(harmony.score)[3] = 'cluster'

lisi.res = rbind(raw.score,sct.score,rpca.score,harmony.score)
lisi_sample = data.frame(val = lisi.res$sample,key = 'sample',method = lisi.res$type)
lisi_study = data.frame(val = lisi.res$study,key = 'study',method = lisi.res$type)
lisi_sample$method = factor(lisi_sample$method,levels = c('Raw','SCT','RPCA','Harmony'))
lisi_study$method = factor(lisi_study$method,levels = c('Raw','SCT','RPCA','Harmony'))

p1 = ggboxplot(lisi_study, x = "method", y = "val",color = 'method',
               palette = "jco",short.panel.labs = FALSE,width = 0.2)+
  labs(y=("LISI score for Study"),x=NULL)+
  # coord_cartesian(ylim = c(0,8))+
  stat_compare_means(comparisons =list(c('Raw','SCT'),c('SCT','RPCA'),c('RPCA','Harmony')),
                     label = "p.format")


p2 = ggboxplot(lisi_sample, x = "method", y = "val",color = 'method',
               palette = "jco",short.panel.labs = FALSE,width = 0.2)+
  labs(y=("LISI score for Sample"),x=NULL)+
  # coord_cartesian(ylim = c(0,8))+
  stat_compare_means(comparisons =list(c('Raw','SCT'),c('SCT','RPCA'),c('RPCA','Harmony')),
                     label = "p.format")
p3 = ggarrange(p1,p2,common.legend = T)
ggsave('LISI.pdf',plot = p3,width = 12,height = 6,dpi = 72)


# clusting and annotation ---------------------------------------------------

library(Nebulosa)
library(clustree)

merge_harmony = readRDS('merge/merge_harmony.rds')

merge_harmony = FindClusters(merge_harmony,resolution = c(1,1.5,2,2.5,3,3.5))
p0 = clustree(merge_harmony,prefix = 'SCT_snn_res.')
ggsave('cluster_res.pdf',plot = p0,width=12,height = 10,dpi = 72)

p1 = FeaturePlot(merge_harmony,features = feature_gene,cols = c("lightgrey" ,"#DE1F1F"))
ggsave('dim.pdf',plot = p1,width=12,height = 10,dpi = 72)

p2 = plot_density(merge_harmony, c("CD8A", "CCR7"),joint = TRUE)[[3]]
p3 = plot_density(merge_harmony, c("CD4", "CCR7"),joint = TRUE)[[3]]
ggsave('CD8A_CD4_CCR1_density.pdf',plot = p2|p3,width=12,height = 4,dpi = 72)

p = DotPlot(merge_harmony,features = manual_label_gene)+coord_flip()
ggsave('manual_label.pdf',plot = p,width = 45,height = 8,dpi = 72,limitsize = F)

# set res as 2.0
Idents(merge_harmony) = factor(merge_harmony$SCT_snn_res.2,levels = c(0:63))

p1 = DotPlot(merge_harmony,features = manual_label_gene)+coord_flip()
ggsave('manual_dot.pdf',plot = p1,width = 49,height = 6,dpi = 72)

# First round manual label------------------

# if existed, use pre-defined annotation in original study


saveRDS(merge_harmony,file = 'merge/merge_label.rds')

merge.metadata = merge_label@meta.data

merge.metadata = merge.metadata[,c(4:8,11,12,14)]

if_else(merge.metadata$cellType)

table(merge.metadata$celltype_major)
table(merge.metadata$Cell_Type)
table(merge.metadata$cellType)

merge.metadata$paper_label = if_else(merge.metadata$celltype_major %in% c('B-cells','Myeloid','Plasmablasts','T-cells'),'immune',
                                     if_else(merge.metadata$celltype_major %in% c('CAFs','Endothelial','PVL'),'stromal',
                                             if_else(merge.metadata$celltype_major %in% c('Cancer Epithelial','Normal Epithelial'),'epi/cancer',
                                                     if_else(merge.metadata$Cell_Type %in% c('Astrocytes','EC-1','EC-2','EC-3','MSC-like-c1','MSC-like-c2','PC-1','PC-2','PC-3','vSMCs'),'stromal',
                                                             if_else(merge.metadata$Cell_Type %in% c('B-c1','B-c2','cDC2:CD1C+/CLEC10A+','MAMs:APOE+','MAMs:S100A8+','T:CD4+:CM1','T:CD4+:CM2','T:CD8+:EM','T:CM','Treg'),'immune',
                                                                     if_else(merge.metadata$Cell_Type %in% c('MTC'),'epi/cancer',
                                                                             if_else(merge.metadata$cellType %in% c('Endothelial','Fibroblast','MSC','Osteoclast'),'stromal',
                                                                                     if_else(merge.metadata$cellType %in% c('Epithelial'),'epi/cancer',
                                                                                             if_else(merge.metadata$cellType %in% c('Macrophage','T'),'immune','not.defined')))))))))
table(merge_label$paper_label)

# atac,cell(part),emboj,both gse,oncogenesis have not defined cell
library(Seurat)
library(dplyr)

merge = readRDS("~/scRNA_BC_metastases/merge/merge_label.rds")
merge.not.define = subset(merge,paper_label =='not.defined')
remove(merge)

sample.name = c('atac',#atac
                'GSM4909308', 'GSM4909310', 'GSM4909312', 'GSM4909314', 'GSM4909316', 'GSM4909318', 'GSM4909321', #emboj
                'B2','C2','D2','D3','E2',#oncogenesis
                'PID1','PID2','PID3', #cell
                'GSE143423','GSE202501' )

library(SingleR)

# ref load
ref.bone <- readRDS("~/scRNA_BC_metastases/singleref/ref.bone.rds")
ref.brain <- readRDS("~/scRNA_BC_metastases/singleref/ref.brain.rds")
ref.primary <- readRDS("~/scRNA_BC_metastases/singleref/ref.primary.rds")

singler.list = list()

library(BiocParallel)

# loop
for (i in 1:18) {
  sample.now = sample.name[i]
  
  # seu object
  seu = subset(merge.not.define, sample == sample.now)
  DefaultAssay(seu) = 'RNA'
  seu = DietSeurat(seu,assays = 'RNA')
  
  print('object loaded')
  
  # set ref by site
  
  if(seu$site[[1]] == 'bone'){
    
    ref = list(primary=ref.primary, bone=ref.bone)
    labels = list(ref.primary$primary.celltype.major,ref.bone$bone.celltype.major)

  }else if(seu$site[[1]] == 'brain'){
    
    ref = list(primary=ref.primary, brain=ref.brain)
    labels = list(ref.primary$primary.celltype.major,ref.brain$brain.celltype.major)

  }else if(seu$site[[1]] == 'lymph'){
    
    ref = list(primary=ref.primary)
    labels = list(ref.primary$primary.celltype.major)

  }else if(seu$site[[1]] == 'primary'){
    
    ref = list(primary=ref.primary)
    labels = list(ref.primary$primary.celltype.major)

  }
  print('reference selected')
  
  # run
  
  seu = as.SingleCellExperiment(seu)
  seu.singler = SingleR(test = seu, ref = ref, labels = labels,
                        de.method="wilcox",de.n = 50, BPPARAM=MulticoreParam(16))
  print('SingleR predicted')
  
  singler.list[[i]] = seu.singler
  names(singler.list[[i]]) = sample.now
  
}

for (i in 1:18) {
  
  study = singler.list[[i]]
  
  study.label = data.frame(x = study@listData[3],y = study@listData[4],study = sample.name[i],cell.name = study@rownames) 
  colnames(study.label)[c(1,2)] = c('label','pruned.label')
  study.label$pruned.label[is.na(study.label$pruned.label)] = 'not.defined'
  
  label.list[[i]] = study.label
  
}

library(data.table)

label.all = rbindlist(label.list)
rownames(label.all) = label.all$cell.name
saveRDS(meta.all,file = 'merge/all.anno.rds')

meta.all$compartment = if_else(meta.all$paper_label %in% c('epi/cancer'),'epi/cancer',
                               if_else(meta.all$paper_label %in% c('immune'),'immune',
                                       if_else(meta.all$paper_label %in% c('stromal'),'stromal',
                                               if_else(meta.all$singler.transfer.label %in% c('Astrocytes','B cells','B-cells','Myeloid','Plasmablasts','T cells','T-cells'),'immune',
                                                       if_else(meta.all$singler.transfer.label %in% c('CAFs','Endothelial','Mesenchymal Stromal Cell-like cells','pericytes','PVL','Vascular smooth muscle cells'),'stromal',
                                                               if_else(meta.all$singler.transfer.label %in% c('Cancer Epithelial','MTC','Normal Epithelial'),'epi/cancer','not.defined'))))))



# second round  -----------------------------------------------------------

# use celltypist!

AddMetadata(merge_label,all.anno)

# immune--------
immune = subset(merge,manual_label=='immune')
DefaultAssay(immune) = 'RNA'
immune = immune %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 2)

saveRDS(immune,file = 'merge/immune.rds')

p1 = dim_batch_plot(immune)
ggsave('batch_immune.pdf',plot = p1,dpi = 72,width = 16,height = 6)

library(harmony)
library(glmGamPoi)

immune_list <- SplitObject(immune, split.by = 'sample')
immune_list <- lapply(X = immune_list, FUN = SCTransform, method = "glmGamPoi")
var_features_ind_sct <- SelectIntegrationFeatures(immune_list, nfeatures = 3000)

immune = merge(x = immune_list[[1]], y = immune_list[2:length(immune_list)])
VariableFeatures(immune) <- var_features_ind_sct

immune <- RunPCA(immune)

immune = RunHarmony(immune,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                        kmeans_init_nstart=20, kmeans_init_iter_max=5000)

immune = immune %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

saveRDS(immune,file = 'merge/immune_harmony.rds')

p2 = dim_batch_plot(immune)
ggsave('batch_immune_harmony.pdf',plot = p2,dpi = 72,width = 16,height = 6)


## use pre-defined anno in paper

meta = immune@meta.data

meta$paper_label_immune = if_else(meta$celltype_major %in% c('B-cells','Plasmablasts'),'B cells',
                                     if_else(meta$celltype_major %in% c('T-cells'),'T cells',
                                             if_else(meta$celltype_major %in% c('Myeloid'),'Myeloid',
                                                     if_else(meta$Cell_Type %in% c('B-c1','B-c2'),'B cells',
                                                             if_else(meta$Cell_Type %in% c('cDC2:CD1C+/CLEC10A+','MAMs:APOE+','MAMs:S100A8+'),'Myeloid',
                                                                     if_else(meta$Cell_Type %in% c('T:CD4+:CM1','T:CD4+:CM2','T:CD8+:EM','T:CM','Treg'),'T cells',
                                                                             if_else(meta$cellType %in% c('T'),'T cells',
                                                                                     if_else(meta$cellType %in% c('Macrophage'),'Myeloid','not.defined'))))))))


immune@meta.data = meta

not.defined.immune = subset(immune,paper_label_immune=='not.defined')
table(not.defined.immune$seurat_clusters)


# cell cycle
immune = CellCycleScoring(immune,
                         s.features = cc.genes$s.genes,
                         g2m.features = cc.genes$g2m.genes)


## t cell subcluster---------

### t cell

table(immune$immune.sub.label)

t.cell = subset(immune,immune.sub.label == 'T cells')
t.cell = t.cell %>% FindNeighbors() %>% FindClusters()

table(t.cell$seurat_clusters)


### tcell manual

meta = t.cell@meta.data
table(meta$celltype_minor)
table(meta$Cell_Type)

meta$paper.label = if_else(meta$celltype_minor %in% c('T cells CD4+'),'T cell CD4+',
                           if_else(meta$celltype_minor %in% c('T cells CD8+'),'T cell CD8+',
                                   if_else(meta$celltype_minor %in% c('NKT cells'),'NK cell',
                                           if_else(meta$celltype_minor %in% c('Cycling T-cells'),'T cell MKI67',
                                                   if_else(meta$celltype_minor %in% c('NKT cells'),'NKT',
                                                           if_else(meta$Cell_Type %in% c('T:CD4+:CM1', 'T:CD4+:CM2'),'T cell CD4+',
                                                                   if_else(meta$Cell_Type %in% c('T:CD8+:EM'),'T cell CD8+',
                                                                           if_else(meta$Cell_Type %in% c('T:CM'),'T cell center memory',
                                                                                   if_else(meta$Cell_Type %in% c('Treg'),'Treg','not.defined')))))))))
t.cell@meta.data = meta

t.not.define = subset(t.cell,paper.label =='not.defined')

meta.not.define = t.not.define@meta.data

## SingleR
library(celldex)
library(SingleR)
ref = celldex::MonacoImmuneData()

t.not.define = t.not.define %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 2)

t.not.define.singler = SingleR(test = t.not.define@assays$RNA@data, ref = ref,labels = ref$label.fine, 
                               clusters = t.not.define$seurat_clusters, fine.tune = TRUE)

new.cluster.ids = t.not.define.singler@listData[["pruned.labels"]]
names(new.cluster.ids) = levels(t.not.define)

t.not.define = RenameIdents(t.not.define,new.cluster.ids)
t.not.define$singler.label = Idents(t.not.define)
Idents(t.not.define) = t.not.define$seurat_clusters


table(t.not.define$singler.label)

t.singler.not.define = subset(t.not.define,singler.label == 'not.defined')
table(t.singler.not.define$site)

t.singler.not.define = t.singler.not.define %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 2)

t.singler.not.define.singler = SingleR(test = t.singler.not.define@assays$RNA@data, ref = ref,labels = ref$label.fine, 
                               clusters = t.singler.not.define$seurat_clusters, fine.tune = TRUE)

plotDeltaDistribution(t.singler.not.define.singler, ncol = 3)


new.cluster.ids = t.singler.not.define.singler@listData[["pruned.labels"]]
names(new.cluster.ids) = levels(t.singler.not.define)

t.singler.not.define = RenameIdents(t.singler.not.define,new.cluster.ids)
t.singler.not.define$singler.label = Idents(t.singler.not.define)
Idents(t.singler.not.define) = t.singler.not.define$seurat_clusters

table(meta.t.nd$singler.label)
meta.t.nd.d = filter(meta.t.nd,singler.label != 'not.defined')

table(meta.t.singler.nd$singler.label)
meta.t.nd = rbind(meta.t.nd.d,meta.t.singler.nd)

t.not.define@meta.data = meta.t.nd

t.define = subset(t.cell,paper.label !='not.defined')

meta.define = t.define@meta.data

meta$t.cell.label = if_else(meta$singler.label %in% 'aa',meta$paper.label,meta$singler.label)


# cd4 ---------------------------------------------------------------------


# cd8 ---------------------------------------------------------------------


# DP ----------------------------------------------------------------------


# dn ----------------------------------------------------------------------



# gd ----------------------------------------------------------------------




# B cell -------------------------------------------------------------------


b.cell = subset(immune,immune.sub.label=='B cells')

meta = b.cell@meta.data
table(meta$celltype_minor)
table(meta$Cell_Type)

meta$paper.label = if_else(meta$celltype_minor %in% c('B cells Memory'),'B cells Memory',
                           if_else(meta$celltype_minor %in% c('B cells Naive'),'B cells Naive',
                                   if_else(meta$celltype_minor %in% c('Plasmablasts'),'Plasmablasts','not.defined')))

b.cell@meta.data = meta

b.not.define = subset(b.cell,paper.label =='not.defined')

#singler
library(celldex)
library(SingleR)
ref = celldex::MonacoImmuneData()
DefaultAssay(b.not.define) = 'RNA'
b.not.define = b.not.define %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 2)

b.not.define.singler = SingleR(test = b.not.define@assays$RNA@data, ref = ref,labels = ref$label.fine, 
                               clusters = b.not.define$seurat_clusters, fine.tune = TRUE)

new.cluster.ids = b.not.define.singler@listData[["pruned.labels"]]
names(new.cluster.ids) = levels(b.not.define)

b.not.define = RenameIdents(b.not.define,new.cluster.ids)
b.not.define$singler.label = Idents(b.not.define)
Idents(b.not.define) = b.not.define$seurat_clusters

meta.singler = b.not.define@meta.data
meta.d = filter(meta,paper.label !='not.defined')
meta= rbind(meta.d,meta.singler)

meta$b.cell.label = if_else(meta$singler.label %in% 'aa',meta$paper.label,meta$singler.label)

table(b.cell$b.cell.label)


# myeloid -----------------------------------------------------------------

# complicate!!

myeloid = subset(immune,immune.sub.label=='Myeloid')

meta = myeloid@meta.data
table(meta$celltype_minor)
table(meta$Cell_Type)

meta$paper.label = if_else(meta$celltype_minor %in% c('Cycling_Myeloid'),'Cycling_Myeloid',
                           if_else(meta$celltype_minor %in% c('DCs'),'DC',
                                   if_else(meta$celltype_minor %in% c('Macrophage'),'Macrophage',
                                           if_else(meta$celltype_minor %in% c('Monocyte'),'Monocyte',
                                                   if_else(meta$Cell_Type %in% c('cDC2:CD1C+/CLEC10A+'),'DC',
                                                           if_else(meta$Cell_Type %in% c('MAMs:APOE+','MAMs:S100A8+'),'Macrophage',
                                                                   if_else(meta$cellType %in% c('Macrophage'),'Macrophage','not.defined')))))))

myeloid@meta.data = meta

myeloid.not.define = subset(myeloid,paper.label =='not.defined')

#singler
library(celldex)
library(SingleR)
ref = celldex::MonacoImmuneData()
DefaultAssay(myeloid.not.define) = 'RNA'
myeloid.not.define = myeloid.not.define %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 2)

myeloid.not.define.singler = SingleR(test = myeloid.not.define@assays$RNA@data, ref = ref,labels = ref$label.fine, 
                               clusters = myeloid.not.define$seurat_clusters, fine.tune = TRUE)

new.cluster.ids = myeloid.not.define.singler@listData[["pruned.labels"]]
names(new.cluster.ids) = levels(myeloid.not.define)

myeloid.not.define = RenameIdents(myeloid.not.define,new.cluster.ids)
myeloid.not.define$singler.label = Idents(myeloid.not.define)
Idents(myeloid.not.define) = myeloid.not.define$seurat_clusters

meta.singler = myeloid.not.define@meta.data
meta.d = filter(meta,paper.label !='not.defined')
meta= rbind(meta.d,meta.singler)

meta$myeloid.label = if_else(meta$singler.label %in% 'aa',meta$paper.label,meta$singler.label)

myeloid@meta.data = meta

saveRDS(myeloid,file = 'merge/myeloid.rds')


# scCATCH 
library(scCATCH)

primary.tissue = c('Breast','Mammary epithelium')
blood.tissue = c('Blood', 'Blood vessel','Plasma')
lymph.tissue = c('Lymph', 'Lymph node', 'Lymphoid tissue')
brain.tissue = c('Brain')
bone.tissue = c('Bone', 'Bone marrow')

myeloid = myeloid %>% RunUMAP(reduction ='harmony',dims = 1:30) %>% 
  FindNeighbors(reduction ='harmony',dims = 1:30) %>% FindClusters(resolution = 2)

immune.catch = createscCATCH(myeloid[['SCT']]@data,cluster = as.character(myeloid$seurat_clusters))
immune.catch <- findmarkergene(object = immune.catch, marker = cellmatch,species = 'Human',
                               tissue = c(primary.tissue,lymph.tissue,brain.tissue,blood.tissue,bone.tissue))
immune.catch <- findcelltype(object = immune.catch)

singler.not.defined.meta = singler.not.defined@meta.data


# merge--------------
meta.mye = myeloid@meta.data
meta.t = t.cell@meta.data
meta.b = b.cell@meta.data

colnames(meta.mye)
colnames(meta.t)
colnames(meta.b)

meta.mye = meta.mye[,c(1:18,23)]
meta.t = meta.t[,c(1:18,24)]
meta.b = meta.b[,c(1:18,23)]

colnames(meta.mye)[19] = 'immune.label'
colnames(meta.t)[19] ='immune.label'
colnames(meta.b)[19] = 'immune.label'

immune.meta = rbind(meta.mye,meta.t,meta.b)



# epi/cancer------------------
epi = subset(merge,manual_label=='epi/cancer')
DefaultAssay(epi) = 'RNA'
epi = epi %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 2)
saveRDS(epi,file = 'merge/epi.rds')

p1 = dim_batch_plot(epi)
ggsave('batch_epi.pdf',plot = p1,dpi = 72,width = 16,height = 6)

library(harmony)
library(glmGamPoi)

epi_list <- SplitObject(epi, split.by = 'sample')
epi_list <- lapply(X = epi_list, FUN = SCTransform, method = "glmGamPoi")
var_features_ind_sct <- SelectIntegrationFeatures(epi_list, nfeatures = 3000)

epi = merge(x = epi_list[[1]], y = epi_list[2:length(epi_list)])
VariableFeatures(epi) <- var_features_ind_sct

epi <- RunPCA(epi)

epi = RunHarmony(epi,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                    kmeans_init_nstart=20, kmeans_init_iter_max=5000)

epi = epi %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 2) %>% RunUMAP(reduction = "harmony", dims = 1:30)

saveRDS(epi,file = 'merge/epi_harmony.rds')
p2 = dim_batch_plot(epi)
ggsave('batch_epi_harmony.pdf',plot = p2,dpi = 72,width = 16,height = 6)

## copykat
library(copykat)
library(tidyverse)

epi.count <- as.matrix(epi@assays$RNA@counts)
epi.cnv <- copykat(rawmat=epi.count, ngene.chr=5, sam.name="epi")

epi <- AddMetaData(epi, metadata = mallignant)
p1 <- DimPlot(sce, group.by = "celltype", label = T)
p2 <- DimPlot(sce, group.by = "copykat.pred") + scale_color_manual(values = c("red", "gray"))
pc <- p1 + p2
ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)

# stromal------------------
stromal = subset(merge,manual_label=='stromal')

DefaultAssay(stromal) = 'RNA'
stromal = stromal %>% FindVariableFeatures() %>%
  NormalizeData() %>% ScaleData() %>% RunPCA() %>%  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 2)
saveRDS(stromal,file = 'merge/stromal.rds')

p1 = dim_batch_plot(stromal)
ggsave('batch_stromal.pdf',plot = p1,dpi = 72,width = 16,height = 6)

library(harmony)
library(glmGamPoi)

stromal_list <- SplitObject(stromal, split.by = 'sample')
stromal_list <- lapply(X = stromal_list, FUN = SCTransform, method = "glmGamPoi")
var_features_ind_sct <- SelectIntegrationFeatures(stromal_list, nfeatures = 3000)

stromal = merge(x = stromal_list[[1]], y = stromal_list[2:length(stromal_list)])
VariableFeatures(stromal) <- var_features_ind_sct

stromal <- RunPCA(stromal)

stromal = RunHarmony(stromal,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                 kmeans_init_nstart=20, kmeans_init_iter_max=5000)

stromal = stromal %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 2) %>% RunUMAP(reduction = "harmony", dims = 1:30)

saveRDS(stromal,file = 'merge/stromal_harmony.rds')
p2 = dim_batch_plot(stromal)
ggsave('batch_stromal_harmony.pdf',plot = p2,dpi = 72,width = 16,height = 6)


# DEG ---------------------------------------------------------------------

library(scDEA)
library(presto)

# cell composition --------------------------------------------------------
library(vegan)

# GSEA --------------------------------------------------------------------
library(GseaVis)
# signature calculate
library(hacksig)
library(CBNplot)

# gene regulation ---------------------------------------------------
# using SCENIC(py) or DoRothEA(R or py) or PISCES

# velocity ----------------------------------------------------------------
# using cellrank or dynamo, or scShaper
library(scshaper)

# cell communcation ---------------------------------------------------
library(CellChat)
library(nichenetr)
library(scseqcomm)

# bulk --------------------------------------------------------------------
library(Scissor)

# CNV ---------------------------------------------------
# using inferCNV or copyKat
#devtools::install_github("broadinstitute/infercnv")

library(inferCNV)
library(copykat)


# umap_visual -----------------------------------------------------------
library(scRNAtoolVis)

# try scHCL ref: Single-cell dissection of intratumoral heterogeneity and lineage diversity in metastatic gastric adenocarcinoma

# also n-gene signature like ref: Single-cell transcriptomic analysis suggests two molecularly distinct subtypes of intrahepatic cholangiocarcinoma

# intratumor heterogeneity? using ref: An algorithm to quantify intratumor heterogeneity based on alterations of gene expression profiles
