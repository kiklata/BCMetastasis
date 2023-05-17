
# sccatch -----------------------------------------------------------------

run.scCATCH.study = function(object,study.name){
  
  library(scCATCH)
  library(harmony)
  library(Seurat)  
  
  primary.tissue = c('Breast','Mammary epithelium')
  blood.tissue = c('Blood', 'Blood vessel','Plasma')
  lymph.tissue = c('Lymph', 'Lymph node', 'Lymphoid tissue')
  brain.tissue = c('Brain')
  bone.tissue = c('Bone', 'Bone marrow')
  
  
  for (i in 1:8) {
    
    study = subset(object,study == study.name[i])
    DefaultAssay(study) = 'RNA'
    
    study = study %>% FindVariableFeatures() %>%
      NormalizeData() %>% ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars = c('sample'),
                                                                  assay.use='RNA', plot_convergence = TRUE,theta = 2.5, 
                                                                  kmeans_init_nstart=20, kmeans_init_iter_max=5000) %>% 
      FindNeighbors(reduction = "harmony",dims = 1:30) %>% 
      FindClusters() %>% RunUMAP(reduction = "harmony",dims = 1:30)
    
    study.catch = createscCATCH(study[['RNA']]@data,cluster = as.character(study$seurat_clusters))
    
    if (study.name[i] %in% c('natgene')) {
      study.catch <- findmarkergene(object = study.catch, marker = cellmatch,species = 'Human',
                                    tissue = c(primary.tissue,blood.tissue))
      study.catch <- findcelltype(object = study.catch)
      
    } else if (study.name[i] %in% c('atac','emboj','oncognesis')) {
      study.catch <- findmarkergene(object = study.catch, marker = cellmatch,species = 'Human',
                                    tissue = c(lymph.tissue,primary.tissue,blood.tissue))
      study.catch <- findcelltype(object = study.catch)
    } else if (study.name[i] %in% c('cell','GSE143423','GSE202501')) {
      study.catch <- findmarkergene(object = study.catch, marker = cellmatch,species = 'Human',
                                    tissue = c(brain.tissue,primary.tissue,blood.tissue))
      study.catch <- findcelltype(object = study.catch)
    } else if (study.name[i] %in% c('bone')) {
      study.catch <- findmarkergene(object = study.catch, marker = cellmatch,species = 'Human',
                                    tissue = c(bone.tissue,primary.tissue,blood.tissue))
      study.catch <- findcelltype(object = study.catch)
    }
    
    study.catch
  }}

