#devtools::install_github("xl0418/ggradar2",dependencies=TRUE)
library(ggradar2)

curated_meta = read.delim('~/Project/MultiOmics/data/skin/res/cellbender_celltype.csv',row.names = 1,sep = ',')
scores = read.delim('~/Project/MultiOmics/data/skin/res/ssgsea_kera.txt')


celltype = 'Keratinocyte'
celltype_index = curated_meta[curated_meta$cl_major %in% c(celltype),] %>% rownames()

celltype_score = scores[celltype_index,] %>% as.data.frame()
celltype_score$CellID = rownames(celltype_score)

sub_meta = curated_meta[c('SampleID','SampleType','cl_minor','cl_subset')]
sub_meta$CellID = rownames(sub_meta)

select_celltype = 'KC-Spinous-SPINK5'

df = left_join(celltype_score,sub_meta, by = 'CellID') #%>% dplyr::filter(., cl_subset == select_celltype)
df$SampleType = if_else(df$SampleType == 'H','ARD','Normal') %>% factor(., levels = c('ARD','Normal'))

inf_sig = c(
  'HALLMARK_UV_RESPONSE_UP',
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_IL2_STAT5_SIGNALING',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  'HALLMARK_TGF_BETA_SIGNALING',
  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
  'HALLMARK_INTERFERON_ALPHA_RESPONSE')


plotdf = df[c(inf_sig,'SampleID','SampleType','cl_minor','cl_subset')] 

group1_plotdf = dplyr::filter(plotdf, SampleType == 'ARD')
group1_plotdf2 = aggregate(group1_plotdf[inf_sig],FUN = mean,by = list(group1_plotdf[,'cl_subset']))
group1_plotdf2$type = 'ARD'

group2_plotdf = dplyr::filter(plotdf, SampleType == 'Normal')
group2_plotdf2 = aggregate(group2_plotdf[inf_sig],FUN = mean,by = list(group2_plotdf[,'cl_subset']))
group2_plotdf2$type = 'Normal'

total_plotdf2 = rbind(group1_plotdf2,group2_plotdf2)

total_plotdf2[inf_sig] = 
  total_plotdf2[inf_sig] %>% apply(.,2,FUN = 
                                     function(data) {
                                       new.data = (data - min(data)) / (max(data) - min(data))
                                       })


colnames(total_plotdf2)[10] = 'group'
rownames(total_plotdf2) = paste0(total_plotdf2$group)

total_plotdf2$Group.1 = NULL


mycolor = c('#4dbbd5','#3c5488','#f39b7f','#00a087')

plot_sig = c(
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_IL2_STAT5_SIGNALING',
  'HALLMARK_TGF_BETA_SIGNALING',
  'HALLMARK_INTERFERON_ALPHA_RESPONSE')

plist = list()
for( i in c('KC-Basal-ITGA6','KC-Suprabasal-PLD1','KC-Spinous-AZGP1','KC-Spinous-SPINK5','KC-Suprabasal-DSC3')){
final = total_plotdf2 %>% dplyr::filter(.,Group.1 == i)
final = final[c(plot_sig,'group')] 

colnames(final)[1:5] = c('Inflammatory','IL6-JAK-STAT3','IL2-STAT5','TGF-β','INF-α')

plist[[i]] = ggradar2(final, polygonfill = FALSE, fullscore = rep(x = 0.75, ncol(final)-1), 
         grid.label.size = 0, gridline.label = c(0,50,100), plot.title = i,
         legend.text.size = 10,axis.label.size = 4
         )+ theme(plot.title = element_text(size = 10))

ggsave(paste0('radar_',i,'.png'),plist[[i]],width = 9.04,height = 5.59, scale = 1, dpi = 300)
}
#library(patchwork)
