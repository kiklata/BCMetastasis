library(infercnv)
options(scipen = 100)

workdir = '~/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM'

obj <- readRDS("~/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM/BrM_decontx_qc.rds")
anno <- read.delim("~/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM/BrM_decontx_celltype.csv",sep = ',',row.names = 1)

obj = AddMetaData(obj,anno)
samplelist = c('BrM1','BrM2','BrM3','BrM4')

for(i in samplelist){
  
seu = subset(obj, SampleID == i)

epi_index = colnames(seu)[seu$cl_major == 'Epi']

normal_index = colnames(seu)[!(seu$cl_major %in% c('Epi','Fb'))]

# ref: A single-cell and spatially resolved atlas of human breast cancers
count = as.matrix(subset(seu, cells = c(epi_index,normal_index))[['RNA']]@counts)

anno_file = anno[c(epi_index,normal_index),'cl_major'] %>% as.data.frame(.,row.name = c(epi_index,normal_index))
ref_group = setdiff(names(table(anno_file[1])),'Epi')
  
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
                             HMM=T,
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
