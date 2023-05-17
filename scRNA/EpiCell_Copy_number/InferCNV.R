library(Seurat)
library(dplyr)
library(infercnv)

all <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/all/all.rds")
all.anno <- readRDS("~/scRNA_BC_metastases/Data/integrate_sc_BC/all/all.anno.rds")
all.anno = all@meta.data
epi = subset(all,celltype.compartment == 'epi/cancer')
# loop

sample.n = names(table(epi$sample))
sample.n = c('cell.pid2','primary.CID44991')
for (i in 1:length(sample.n)) {

  seu = subset(all,sample == sample.n[i])
  seu1 = subset(seu,celltype.compartment == 'epi/cancer')
  
 if(ncol(seu1)<500){
   dir.create(sample.n[i])
  
  # anno file
  
  anno = filter(all.anno,sample == sample.n[i])
  anno = anno[,c('cell.id','celltype.compartment')]
  anno$celltype.compartment = if_else(anno$celltype.compartment %in% c('epi/cancer'),'tumor','normal')
  norm.cell = rownames(filter(anno,celltype.compartment =='normal'))
  
  if(length(norm.cell)>1000){
    norm.cell = norm.cell[sample(length(norm.cell),1000)]
  }
  
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
                               out_dir= paste0(sample.n[i]), #no_plot = TRUE,
                               num_threads = 16,
                               denoise=TRUE,
                               HMM=T,
                               HMM_type = 'i6',HMM_transition_prob = 1e-06
                               ))

  }else if(ncol(seu1)>500){
    dir.create(sample.n[i])
    
    # anno file
    
    anno = filter(all.anno,sample == sample.n[i])
    anno = anno[,c('cell.id','celltype.compartment')]
    anno$celltype.compartment = if_else(anno$celltype.compartment %in% c('epi/cancer'),'tumor','normal')
    norm.cell = rownames(filter(anno,celltype.compartment =='normal'))
    
    if(length(norm.cell)>2000){
      norm.cell = norm.cell[sample(length(norm.cell),2000)]
    }
    
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
    try(infercnv::run(infercnv_obj,
                                 cutoff=0.1, 
                                 # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir= paste0(sample.n[i]), #no_plot = TRUE,
                                 num_threads = 16,
                                 denoise=TRUE,
                                 HMM=T,
                                 HMM_type = 'i6',HMM_transition_prob = 1e-06
    ))
  }
}

