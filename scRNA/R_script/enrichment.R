#enrichment
seu <- readRDS("~/Project/MultiOmics/data/skin/res/cellbender_anno_count.rds")

seu_sub = seu %>% subset(.,cl_major == 'Keratinocyte')
count = seu_sub@assays$RNA@counts

library(GSVA)
library(GSEABase)

ncores = 4

gmt_file = '~/Reference/gmt/KEGG_metabolism_nc.gmt'

geneSets <- getGmt(gmt_file) 
gsva_es <- gsva(as.matrix(count), min.sz = 3,geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) 
signature_exp<-as.matrix(gsva_es) %>% t()
saveRDS(signature_exp, file = 'meta_score.rds')

# Msigdb:C7 immune 


# metabolism
# scMetabolism

# run deg first
source("~/Project/MultiOmics/code/func/rungsea.R")

res = markers %>% run.gsea(.,gmt = '/home/zhepan/Reference/gmt/KEGG_metabolism_nc.gmt', 
                           clusters ='Macrophage-TREM2', method = 'gsea' )

df = res$res@result
df$group = if_else(df$NES>0, 'ARD','Normal')

library(GseaVis)

sigs = c(
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_IL2_STAT5_SIGNALING')

for(i in sigs){
  p = gseaNb(object = res$res,geneSetID = 'Propanoate metabolism',addPval = F
             ,subPlot = 2)
  ggsave(paste0(i,'_gsea.png'),p, width = 4.98, height = 4.26)
}

