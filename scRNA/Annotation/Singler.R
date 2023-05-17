source('obj_data.R')
# define singler ----------------------------------------------------------
library(Seurat)

run.singler = function(obj.list,sample.list,savepath){
  
  library(SingleR)
  library(celldex)
  library(BiocParallel)
  library(SingleCellExperiment)
  
  ref = celldex::MonacoImmuneData()
  singler.list = list()
  
  for (i in 1:length(obj.list)) {
    
    seu = obj.list[[i]]
    seu = CreateSeuratObject(seu[['RNA']]@counts)
    seu = as.SingleCellExperiment(seu)
    seu.singler = SingleR(test = seu, ref = ref, labels = ref$label.fine, fine.tune = T,
                          de.method="wilcox",de.n = 50, BPPARAM=MulticoreParam(16))
    
    print('SingleR predicted ')
    singler.list[sample.list[i]] = seu.singler
    
  }
  dir.create(savepath)
  saveRDS(singler.list,file = paste0(savepath,'/singler.list.rds'))
}


# run ---------------------------------------------------------------------

obj.list = c(atac, 
             EMBOJ.GSM4909308,EMBOJ.GSM4909310,EMBOJ.GSM4909312,EMBOJ.GSM4909314,
             EMBOJ.GSM4909316,EMBOJ.GSM4909318,EMBOJ.GSM4909321,
             oncogenesis.B2,oncogenesis.C2,oncogenesis.D2,oncogenesis.D3,oncogenesis.E2, 
             cell.pid1,cell.pid2,cell.pid3, 
             GSE143423,
             GSE202501,
             science.GSM4555888,science.GSM4555889,science.GSM4555891,
             bom1,bom2,
             CID3586,CID3838,CID3921,CID3941,CID3946,CID3948,CID3963,
             CID4040,CID4066,CID4067,CID4290A,CID4398,CID44041,CID4461,CID4463,CID4465,CID4471,CID4495,CID44971,CID44991,
             CID4513,CID4515,CID45171,CID4523,CID4530N,CID4535)

sample.list = c('atac', 
                'EMBOJ.GSM4909308','EMBOJ.GSM4909310','EMBOJ.GSM4909312','EMBOJ.GSM4909314',
                'EMBOJ.GSM4909316','EMBOJ.GSM4909318','EMBOJ.GSM4909321',
                'oncogenesis.B2','oncogenesis.C2','oncogenesis.D2','oncogenesis.D3','oncogenesis.E2', 
                'cell.pid1','cell.pid2','cell.pid3', 
                'GSE143423',
                'GSE202501',
                'science.GSM4555888','science.GSM4555889','science.GSM4555891',
                'bom1','bom2',
                'CID3586','CID3838','CID3921','CID3941','CID3946','CID3948','CID3963',
                'CID4040','CID4066','CID4067','CID4290A','CID4398','CID44041','CID4461','CID4463','CID4465','CID4471','CID4495','CID44971','CID44991',
                'CID4513','CID4515','CID45171','CID4523','CID4530N','CID4535')

savepath = '~/scRNA_BC_metastases/singler'

run.singler(obj.list = obj.list,sample.list = sample.list,savepath=savepath)

# tcell ---------------------

tcell.list = SplitObject(tcell,'sample')
sample.list = names(table(tcell$sample))
savepath = "~"
options(warn = -1)
run.singler(tcell.list,sample.list,savepath)
