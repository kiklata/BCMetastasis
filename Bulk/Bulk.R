fpkm = read.delim('MET500_geneExpression_M.mx.txt.gz')
probeMap_gencode.v23.annotation.gene <- read.delim("~/scRNA_BC_metastases/Data/bulk_BC/MET500/probeMap_gencode.v23.annotation.gene.probemap")
MET500_geneExpression_M.meta.plus <- read.delim("~/scRNA_BC_metastases/Data/bulk_BC/MET500/MET500_geneExpression_M.meta.plus.txt")
colnames(fpkm)[1] = 'id'
gene = probeMap_gencode.v23.annotation.gene[,c(1,2)]
fpkm = left_join(fpkm,gene,by ='id')
metadata = filter(MET500_geneExpression_M.meta.plus,cohort == 'BRCA')

fpkm = fpkm[,colnames(fpkm) %in% metadata$Sample_id]

name.col = colnames(fpkm)
meta.name = metadata$Sample_id

meta.name = gsub("-",'\\.',meta.name)
metadata$Sample_id = meta.name

exp = fpkm[,colnames(fpkm) %in% metadata$Sample_id]
exp$ENSEMBL = fpkm$id
exp$SYMBOL = fpkm$gene

saveRDS(exp,file = 'FPKM_mBC.rds')
saveRDS(metadata,file = 'mBC_metadata.rds')
