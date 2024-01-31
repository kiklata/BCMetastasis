nmf_count <- read.table("~/Project/scRNA_BC_metastases/Result/nmf/tumor/nmf_count.txt", quote="\"", comment.char="")

obj_list = list()

for (i in nmf_count$V1) {
  count_p = paste0(i,'/count_rm_ribo.txt')
  df = read.delim(count_p) %>% t()
  seu = CreateSeuratObject(df)
  names = i %>% gsub('/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/','',.) %>% gsub('/','_',.)
  source = names %>% stringr::str_split(.,'_') %>% unlist() %>% .[1]
  site = names %>% stringr::str_split(.,'_') %>% unlist() %>% .[2]
  seu$SampleID = names
  seu$Source = source
  seu$Site = site
  obj_list[[names]] = seu
}
saveRDS(obj_list, file = '/home/zhepan/Project/scRNA_BC_metastases/Result/nmf/tumor/mc_obj_list.rds')

seu = merge(obj_list[[1]], obj_list[2:60],add.cell.ids = names(obj_list))
seu = JoinLayers(seu)
saveRDS(seu, file = '/home/zhepan/Project/scRNA_BC_metastases/Result/nmf/tumor/mc_merge.rds')