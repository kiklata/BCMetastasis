library(infercnv)
options(scipen = 100)

workdir = '~/Project/MultiOmics/data/snRNA/Object/summary/test'

cellranger_filter <- readRDS("~/Project/MultiOmics/data/snRNA/Object/summary/cellranger_filter_count.rds")
anno.matrix <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/anno.tsv",sep = '\t',row.names = 1)

cellranger_filter = AddMetaData(cellranger_filter,anno.matrix)
samplelist = c('P1013S2','P1015S2','P1018S1')

for(i in samplelist){
  
seu = subset(cellranger_filter, SampleID == i)

epi_index = colnames(seu)[seu$manual_celltype_annotation == 'Epithelial']
normal_index = colnames(seu)[!(seu$manual_celltype_annotation %in% c('Epithelial','PVL','CAF'))] %>% sample(.,length(epi_index))

# ref: A single-cell and spatially resolved atlas of human breast cancers
count = as.matrix(subset(seu, cells = c(epi_index,normal_index))[['RNA']]@counts)

anno_file = anno.matrix[c(epi_index,normal_index),'manual_celltype_annotation'] %>% as.data.frame(.,row.name = c(epi_index,normal_index))
ref_group = setdiff(names(table(anno_file[1])),'Epithelial')
  
if(dir.exists(file.path(workdir,'cnv', i))==F){ 
  dir.create(file.path(workdir,'cnv', i), recursive = T )}

write.table(anno_file, file= file.path(workdir,'cnv', i,"infercnv_anno.txt"), col.names = F, row.names = T, quote=F, sep="\t")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = count,
                                    annotations_file= file.path(workdir,'cnv', i,"infercnv_anno.txt"),
                                    gene_order_file= "/home/zhepan/Reference/gencode_v32_gene_pos_gene_name.txt",
                                    ref_group_names= ref_group) 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=file.path(workdir,'cnv',i),  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   # cluster
                             denoise=T,
                             analysis_mode="subclusters",
                             tumor_subcluster_partition_method='random_trees',
                             HMM=T,write_phylo = TRUE,
                             num_threads = 8
)
}
#run.final.infercnv_obj <- readRDS("~/Project/MultiOmics/data/snRNA/Object/P1013S2/cnv/run.final.infercnv_obj.rds")
#cnv = run.final.infercnv_obj@expr.data

#scores=apply(cnv,2,function(x){ sum(x < 0.95 | x > 1.05)/length(x) })

# infercna based ------------------------------------------------------

# ref: An integrative model of cellular states, plasticity, and genetics for glioblastoma. 
#library(infercna) # kiklata/infercna

#df = data.frame(sig = cnaSignal(cnv,gene.quantile = 0.9),
#                cor = cnaCor(cnv,gene.quantile = 0.9))

#cnaScatterPlot(cnv, gene.quantile = 0.9, main = 'threshold: 0.9')

#Modes = infercna::findMalignant(cnv - 1, gene.quantile = .9, verbose = T)
#mal_index = Modes$malignant

#anno.matrix$infercnv = if_else((rownames(anno.matrix) %in% mal_index) & (anno.matrix$manual_celltype_annotation == 'Epithelial'),
#                               'malignant',anno.matrix$manual_celltype_annotation)

#saveRDS(anno.matrix, file = file.path(workdir,'cnv', "mal.anno.matrix.rds"))
