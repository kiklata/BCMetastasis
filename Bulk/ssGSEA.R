library(GSVA)
library(GSEABase)
library(ggpubr)

# MET500 ------------------------------------------------------------------

FPKM_mBC <- readRDS("~/scRNA_BC_metastases/Data/bulk_BC/MET500/FPKM_mBC.rds")
mBC_metadata <- readRDS("~/scRNA_BC_metastases/Data/bulk_BC/MET500/mBC_metadata.rds")
#top5gene <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/top5gene.rds")
metaProgram.gene <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/metaProgram.gene.rds")
metaProgram.gene.list = list()

for(i in 1:10){
  gene = filter(metaProgram.gene,plot.pro == paste0('P',i))$metaProgram.gene
  metaProgram.gene.list[[i]] = gene
  names(metaProgram.gene.list)[i] = paste0('P',i)
}

FPKM_mBC = FPKM_mBC[!duplicated(FPKM_mBC$SYMBOL),]
FPKM_mBC = FPKM_mBC[!is.na(FPKM_mBC$SYMBOL),]

rownames(FPKM_mBC) = FPKM_mBC$SYMBOL

FPKM_mBC$ENSEMBL = NULL
FPKM_mBC$SYMBOL = NULL

es<-gsva(as.matrix(FPKM_mBC),metaProgram.gene.list,method='ssgsea')

es = as.data.frame(es)
saveRDS(es,file = '~/scRNA_BC_metastases/Result/Score/MET500.score.es.rds')

score.es <- readRDS("~/scRNA_BC_metastases/Result/Score/MET500.score.es.rds")

score.es = as.data.frame(t(score.es))
score.es$Sample_id = rownames(score.es)

plot.score = left_join(score.es,mBC_metadata,by ='Sample_id')

plot.score$biopsy_tissue = if_else(plot.score$biopsy_tissue == 'bone_marrow','bone',
                                   if_else(plot.score$biopsy_tissue == 'lymph_node','lymph',plot.score$biopsy_tissue))

mycomparsion = list(c('lymph','brain'),c('bone','brain'))
plot.score = filter(plot.score,biopsy_tissue != 'soft_tissue')
plot.score = filter(plot.score,biopsy_tissue != 'pancreas')
plot.score = filter(plot.score,biopsy_tissue != 'breast')

p1 = ggplot(plot.score, aes(x = biopsy_tissue, y = P1,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p2 = ggplot(plot.score, aes(x = biopsy_tissue, y = P2,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p3 = ggplot(plot.score, aes(x = biopsy_tissue, y = P3,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p4 = ggplot(plot.score, aes(x = biopsy_tissue, y = P4,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p5 = ggplot(plot.score, aes(x = biopsy_tissue, y = P5,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p6 = ggplot(plot.score, aes(x = biopsy_tissue, y = P6,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p7 = ggplot(plot.score, aes(x = biopsy_tissue, y = P7,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p8 = ggplot(plot.score, aes(x = biopsy_tissue, y = P8,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p9 = ggplot(plot.score, aes(x = biopsy_tissue, y = P9,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  
p10 = ggplot(plot.score, aes(x = biopsy_tissue, y = P10,color = biopsy_tissue))+
  geom_boxplot(outlier.shape = 1)+theme_classic() + 
  stat_compare_means(method = 't.test',comparisons = mycomparsion)  

#library(patchwork)

(p1|p3|p4)/(p7|p9|p10)+plot_layout(guides = 'collect')


# TCGA surv--------------------------------------------------------------------
#top5gene <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/top5gene.rds")
TCGA.BRCA.htseq_fpkm.uq.tsv <- read.delim("~/scRNA_BC_metastases/Data/bulk_BC/TCGA/TCGA-BRCA.htseq_fpkm-uq.tsv.gz")
gencode.v22.annotation.gene <- read.delim("~/scRNA_BC_metastases/Data/bulk_BC/TCGA/gencode.v22.annotation.gene.probeMap")
metaProgram.gene <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/nmf/metaProgram.gene.rds")
metaProgram.gene.list = list()

for(i in 1:10){
  gene = filter(metaProgram.gene,plot.pro == paste0('P',i))$metaProgram.gene
  metaProgram.gene.list[[i]] = gene
  names(metaProgram.gene.list)[i] = paste0('P',i)
}

colnames(gencode.v22.annotation.gene)[1] = 'Ensembl_ID'

genetrans = gencode.v22.annotation.gene[,c(1,2)]

BRCA = left_join(TCGA.BRCA.htseq_fpkm.uq.tsv,genetrans,by = 'Ensembl_ID')

BRCA$Ensembl_ID = NULL
BRCA = BRCA[!duplicated(BRCA$gene),]
BRCA = BRCA[!is.na(BRCA$gene),]

rownames(BRCA) = BRCA$gene
BRCA$gene = NULL

es<-gsva(as.matrix(BRCA),metaProgram.gene.list,method='ssgsea')

es = as.data.frame(es)
saveRDS(es,file = 'TCGA.score.es.rds')

TCGA.score.es <- readRDS("~/scRNA_BC_metastases/Result/Score/TCGA.score.es.rds")
TCGA_metadata <- read.delim("~/scRNA_BC_metastases/Data/bulk_BC/TCGA/TCGA-BRCA.survival.tsv")
TCGA_metadata$tumor = if_else(substring(TCGA_metadata$sample,first = 14,last = 15) == 11,'normal','tumor')

TCGA.score.es = as.data.frame(t(TCGA.score.es))
TCGA.score.es$sample = rownames(TCGA.score.es)

TCGA_metadata$sample = gsub('-','.',TCGA_metadata$sample)

plot.score = left_join(TCGA.score.es,TCGA_metadata,by ='sample')

plot.score = na.omit(plot.score)

plot.score = filter(plot.score,tumor =='tumor')

# km plot

library(survival)
library(survminer)

# filter 4000 day survival
plot.score.filter = filter(plot.score,OS.time<=4000)
res.cut <- surv_cutpoint(plot.score.filter, time = "OS.time", event = "OS",
                         variables = c("P1","P2","P3",'P4','P5',
                                       'P6','P7','P8','P9','P10'))

#plot(res.cut, "P1", palette = "npg")

res.cat <- surv_categorize(res.cut)

#summary(res.cut)\
fit.list = list()
fit.list[[1]] = survfit(Surv(OS.time, OS)~P1,data = res.cat)
fit.list[[2]] = survfit(Surv(OS.time, OS)~P2,data = res.cat)
fit.list[[3]] = survfit(Surv(OS.time, OS)~P3,data = res.cat)
fit.list[[4]] = survfit(Surv(OS.time, OS)~P4,data = res.cat)
fit.list[[5]] = survfit(Surv(OS.time, OS)~P5,data = res.cat)
fit.list[[6]] = survfit(Surv(OS.time, OS)~P6,data = res.cat)
fit.list[[7]] = survfit(Surv(OS.time, OS)~P7,data = res.cat)
fit.list[[8]] = survfit(Surv(OS.time, OS)~P8,data = res.cat)
fit.list[[9]] = survfit(Surv(OS.time, OS)~P9,data = res.cat)
fit.list[[10]] = survfit(Surv(OS.time, OS)~P10,data = res.cat)

p1 = ggsurvplot(fit.list[[1]],res.cat,pval = T)
p2 = ggsurvplot(fit.list[[2]],res.cat,pval = T)
p3 = ggsurvplot(fit.list[[3]],res.cat,pval = T)
p4 = ggsurvplot(fit.list[[4]],res.cat,pval = T)
p5 = ggsurvplot(fit.list[[5]],res.cat,pval = T)
p6 = ggsurvplot(fit.list[[6]],res.cat,pval = T)
p7 = ggsurvplot(fit.list[[7]],res.cat,pval = T)
p8 = ggsurvplot(fit.list[[8]],res.cat,pval = T)
p9 = ggsurvplot(fit.list[[9]],res.cat,pval = T)
p10 = ggsurvplot(fit.list[[10]],res.cat,pval = T)

#fit.list = list()
#for (i in 1:10) {
#  fit.list[[i]] = survfit(Surv(OS.time, OS)~res.cat[,i+2],data = res.cat)
#  names(fit.list[[i]][["strata"]]) = c('High','Low')
#}
#ggsurvplot_list(fit.list,data = res.cat)
plot.list = list(p1,p6,p2,p7,p3,p8,p4,p9,p5,p10)
res = arrange_ggsurvplots(plot.list,ncol = 5,nrow = 2)

ggsave('~/scRNA_BC_metastases/Result/Score/TCGAsurv.png',res,width = 20,height = 6)
