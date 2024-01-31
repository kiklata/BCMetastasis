nmf_count <- read.table("~/Project/scRNA_BC_metastases/Result/nmf/tumor/nmf_count.txt", quote="\"", comment.char="")
for (i in nmf_count$V1) { #nmf_count$V1
  count_p = paste0(i,'/count.txt')
  df = read.delim(count_p,sep = '\t')
  rp_gene = colnames(df) %>% grep('^RP[S/L]',.,value = T)
  df_new = df[,!colnames(df) %in% rp_gene]
  write.table(df_new, file = paste0(i,'/count_rm_ribo.txt'),quote = F,sep = '\t',row.names = T, col.names = T)
  print(i)
}
