all <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/all/all.rds")
all.anno = all@meta.data
epi = subset(all,celltype.compartment == 'epi/cancer')
sample.n = names(table(epi$sample))

setwd('~/scRNA_BC_metastases/Data/integrate_sc_BC/epi/inferCNV')

for (i in 1:length(sample.n)) {
  
  seu = subset(all,sample == sample.n[i])
  seu1 = subset(seu,celltype.compartment == 'epi/cancer')
  
    dir.create(sample.n[i])
    
    # anno file
    
    anno = seu@meta.data
    anno = anno[,c('cell.id','celltype.compartment')]
    anno$celltype.compartment = if_else(anno$celltype.compartment %in% c('epi/cancer'),'tumor',
                                        if_else(anno$celltype.compartment %in% c('immune'),'normal','filterd'))
    norm.cell = rownames(filter(anno,celltype.compartment =='normal'))
    
    norm = subset(seu, cells =  norm.cell)
    
    count = cbind(seu1[['RNA']]@counts,norm[['RNA']]@counts)
    
    anno = anno[anno$cell.id %in% colnames(count),]
    write.table(anno, 
                file= paste0(sample.n[i],"/anno.txt"), col.names = F, row.names = F, quote=F, sep="\t")
    
    
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix = count,
                                        annotations_file= paste0(sample.n[i],"/anno.txt"),
                                        gene_order_file= "~/scRNA_BC_metastases/Analysis/scRNA/EpiCell_Copy_number/gencode_v32_gene_pos_gene_name.txt",
                                        ref_group_names=c("normal")) 
    
    options(scipen = 100)
    try( infercnv::run(infercnv_obj,
                       cutoff=0.1, 
                       # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                       out_dir= paste0(sample.n[i]), no_plot = TRUE,
                       num_threads = 16,
                       denoise=TRUE,
                       HMM=T,
                       HMM_type = 'i6',HMM_transition_prob = 1e-06
    ))
    }
