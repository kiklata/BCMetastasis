

all_path = read.delim('/home/zhepan/Project/scRNA_BC_metastases/Result/nmf/tumor/nmf_cell_number.txt', sep = ' ')

all_nmf_gene_score = list()

for (i in 1:nrow(all_path)) {
  nmf_path = all_path[i, 'path']
  sample = all_path[i, 'sample']
  source = all_path[i, 'source']
  
  hvg = read.delim(
    paste0(nmf_path, '/cNMF/cNMF','.overdispersed_genes.txt'),
    sep = '\t',
    header = F
  )
  
  threshold = '0_02' # 2_0 0_02
  
  k_range = seq(4, 10)
  
  hvg_score_list = list()
  
  for (k in k_range) {
    gene_score = read.delim(file.path(
      nmf_path,
      paste0('/cNMF/cNMF','.gene_spectra_score.k_', k, '.dt_', threshold, '.txt')
    ), sep = '\t')
    gene_score$X = paste(sample,
                         'k',
                         k,
                         'threshold',
                         threshold,
                         'module',
                         gene_score$X,
                         sep = '_')
    gene_score = tibble::column_to_rownames(gene_score, 'X') %>% t() %>% as.data.frame()
    rownames(gene_score) = rownames(gene_score) %>% gsub('\\.', '-', .)
    hvg_score_list[[paste0(sample, '_', k)]] = gene_score[hvg$V1, ]
  }
  
  hvg_score_list = unname(hvg_score_list)
  hvg_score = do.call(cbind, hvg_score_list)
  
  colnames(hvg_score) = paste0(source, '_', colnames(hvg_score))
  
  all_nmf_gene_score[[paste0(source,'_',sample)]] = hvg_score
  print(sample)
}

saveRDS(all_nmf_gene_score, file = '/home/zhepan/Project/scRNA_BC_metastases/Result/nmf/tumor/cNMF_gene_score_0_02.rds')
