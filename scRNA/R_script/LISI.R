
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
