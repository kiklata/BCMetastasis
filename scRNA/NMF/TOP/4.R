library(tidyverse)

mat.files=dir("./res2/",pattern = "dt_0_02.txt$")

all.mat=data.frame()
for (fi in mat.files) {
  tmp.mat=read.table(paste0("./res2/",fi),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  tmp.mat=as.data.frame(t(tmp.mat))
  sampleid=str_replace(fi,"\\..*$","")
  colnames(tmp.mat)=paste(sampleid,colnames(tmp.mat),sep = "_")
  tmp.mat$gene=rownames(tmp.mat)
  
  if (sampleid == "HNSCC22") {
    all.mat=tmp.mat
  }else{
    all.mat=all.mat %>% full_join(tmp.mat,by="gene") #元素的并集进行合并
  }
}

# 对于某一个模块
signature.programs=c("HNSCC6_3","HNSCC22_4","HNSCC5_3")
signature.loading=all.mat[,c("gene",signature.programs)]

used.gene=c()
for (pi in signature.programs) {
  tmp.df=signature.loading[,c("gene",pi)]
  tmp.loading=tmp.df[,2]
  names(tmp.loading)=tmp.df[,1]
  
  tmp.loading=tmp.loading[!is.na(tmp.loading)]
  used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
}
used.gene=unique(used.gene)

signature.loading=signature.loading[signature.loading$gene %in% used.gene,]
rownames(signature.loading)=signature.loading$gene
signature.loading$gene=NULL
signature.loading[is.na(signature.loading)]<-0
signature.loading$total_loading=rowSums(signature.loading)
signature.loading$average_loading=signature.loading$total_loading / length(signature.programs)

signature.loading=signature.loading%>%arrange(desc(average_loading))
head(rownames(signature.loading),30)

# library(Seurat)
# 
# tmp.gene=c("TUBA1B", "HMGB2",  "TOP2A",  "UBE2C",  "NUSAP1", "MKI67",  "PBK",
#            "BIRC5",  "CDK1",   "AURKB",  "UBE2T",  "CKS1B",  "H2AFZ",  "TPX2",   
#            "TK1",    "CCNA2",  "GTSE1", "CEP55",  "KIF23",  "CENPF",  "CKS2",
#            "HMGB1",  "PTTG1",  "CDCA3",  "CDKN3",  "PRC1",   "NUF2",   "CCNB1",
#            "CCNB2",  "FOXM1")
# sum(tmp.gene %in% c(Seurat::cc.genes.updated.2019$s.genes,Seurat::cc.genes.updated.2019$g2m.genes))
