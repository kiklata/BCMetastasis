## ssGSEA
library(ggpubr)

# RUN scMetabolism
curated_meta = read.delim('~/Project/MultiOmics/data/snRNA/Object/summary/annotation/myeloid_celltype.csv',row.names = 1,sep = ',')
scores = read.delim('~/Project/MultiOmics/data/skin/res/ssgsea_kera.txt')


celltype = c('Macrophage','Monocyte')
celltype_index = curated_meta[curated_meta$myeloid_major %in% c(celltype),] %>% rownames()

celltype_score = scores[celltype_index,] %>% as.data.frame()
celltype_score$CellID = rownames(celltype_score)

sub_meta = curated_meta[c('SampleTimepoint','Myeloid_minor')]
sub_meta$CellID = rownames(sub_meta)

select_celltype = 'KC-Spinous-SPINK5'

df = left_join(celltype_score,sub_meta, by = 'CellID') #%>% dplyr::filter(., cl_subset == select_celltype)
#df$SampleType = if_else(df$SampleType == 'H','ARD','Normal') %>% factor(., levels = c('ARD','Normal'))

# scores by celltype---------------------
meta_interest = 'Fatty acid elongation'
ggplot(data = df,aes(x = Myeloid_minor, y = .data[[meta_interest]]))+
  geom_violin(aes(fill = Myeloid_minor),width = 0.7)+
  geom_boxplot(aes(fill = Myeloid_minor),width = 0.2, position = position_dodge(0.7),outlier.shape = NA)+
  #scale_fill_manual(values = timepoint_color_p)+
  labs(x = '',y = 'Score', title = meta_interest)+
  theme_classic()+
  NoLegend()+
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 8))

library(ggridges)
ggplot(data = df,aes(x = .data[[meta_interest]], y = Myeloid_minor, fill = Myeloid_minor))+
  geom_density_ridges()+
  scale_fill_manual(values = minor_color_p)+
  labs(x = 'Score',y = '', title = meta_interest)+
  theme_classic()+
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 8))



# scores by sampletype ----------------------------------------------------
limmstest = function(df,SampleType){
  require(limma)
  group <- df[['SampleType']] %>% as.factor()
  desigN <- model.matrix(~ 0 + group) 
  colnames(desigN) <- levels(group)
  
  contrast.matrix <- makeContrasts(ARD-Normal,levels = desigN)
  
  fit = lmFit(df[,c(1:134)] %>% t(), desigN)  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  diff=topTable(fit2,coef=1,number=Inf,adjust = 'fdr',p.value = 0.001)
}

diff = limmstest(df,'SampleType')
diff[inf_sig,]

df %>% group_by(SampleType) %>% summarise(., mean = mean(HALLMARK_TGF_BETA_SIGNALING))
 
inf_sig = c(
  'HALLMARK_UV_RESPONSE_UP',
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_IL2_STAT5_SIGNALING',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  'HALLMARK_TGF_BETA_SIGNALING',
  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
  'HALLMARK_INTERFERON_ALPHA_RESPONSE')

plist = list()
for(meta_interest in inf_sig){

  plist[[meta_interest]] = 
  ggplot(data = df,aes(x = SampleType, y = .data[[meta_interest]]))+
  geom_violin(aes(fill = SampleType),width = 0.7)+
  geom_boxplot(aes(fill = SampleType),width = 0.2, position = position_dodge(0.7),outlier.shape = NA)+
  #scale_fill_manual(values = timepoint_color_p)+
  labs(x = '',y = 'Score', title = meta_interest %>% gsub('HALLMARK_','',.))+
  stat_compare_means(comparisons = list(c('ARD','Normal')),
                     method = "wilcox.test",
                     label = "p.value",tip.length = 0,size = 2,lwd = 0.5)+
  theme_classic()+
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 8))
  ggsave(filename = paste0('~/Project/MultiOmics/code/figures/all_kera_',meta_interest,'_violin.png'),plist[[meta_interest]],width = 2.27, height = 2.38, scale = 1)
}
