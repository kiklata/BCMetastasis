obj_create = function(type,
                      datapath,
                      savepath,
                      NeoRadRes = NULL,
                      SampleID = NULL,
                      Site = NULL,
                      Kit = NULL,
                      Age = NULL,
                      meta_ER = NULL,
                      meta_PR = NULL,
                      meta_HER = NULL,
                      meta_KI67 = NULL) {
  # only for test -----------------------------------------------------------
  #datapath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM'
  #savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BoM'
  #SampleID = 'BoM1'
  #Site = 'Bone'
  #Kit = 'sc'
  #Age = ''
  #ER = ''
  #PR = ''
  #HER = ''
  #KI67 = ''
  
  
  # load pkgs ---------------------------------------------------------------
  
  suppressMessages(library(Seurat))
  suppressMessages(library(dplyr))
  suppressMessages(library(scDblFinder))
  
  
  # load data ---------------------------------------------------------------
  
    cellranger_path = file.path(datapath, SampleID, 'filtered_feature_bc_matrix.h5')
    soupx_path = file.path(datapath, SampleID, 'soupx.rds')
    decontx_path = file.path(datapath, SampleID, 'decontx.rds')
    cellbender_path = file.path(
      datapath,
      SampleID,
      'cellbender_feature_bc_matrix_filtered.h5'
    )
    
    if (type == 'cellranger') {
      count = Read10X_h5(cellranger_path)
      seu = CreateSeuratObject(
        count,
        project = SampleID,
        min.cells = 3,
        min.features = 200
      )
    } else if (type == 'soupx') {
      seu = readRDS(soupx_path)
    } else if (type == 'cellbender') {
      count = scCustomize::Read_CellBender_h5_Mat(cellbender_path)
      seu = CreateSeuratObject(
        count,
        project = SampleID,
        min.cells = 3,
        min.features = 200
      )
    } else if (type == 'decontx') {
      seu = readRDS(decontx_path)
    }
    
    seu$SampleID = SampleID
    seu$Site = Site
    seu$Kit = Kit
    seu$Age = Age
    seu$meta_ER = meta_ER
    seu$meta_PR = meta_PR
    seu$meta_HER = meta_HER
    seu$meta_KI67 = meta_KI67

  
  # doublet detect ----------------------------------------------------------
  
  set.seed(42)
  sce = scDblFinder(as.SingleCellExperiment(seu))
  seu$scDblFinder.score = sce$scDblFinder.score
  seu$scDblFinder.class = sce$scDblFinder.class
  
  
  # some qc info ------------------------------------------------------------
  
  seu$percent_mt = PercentageFeatureSet(seu, pattern = '^MT-')
  seu$percent_hb = PercentageFeatureSet(seu, pattern = "^HB[^(P)]")
  seu$percent_rb = PercentageFeatureSet(seu, pattern = "^RP[SL]")
  
  cyclegenes = read.delim('/home/zhepan/Reference/regev_lab_cell_cycle_genes.txt')
  seu <- NormalizeData(seu)
  
  seu = CellCycleScoring(seu, s.features = cyclegenes[1:42, 1], g2m.features = cyclegenes[43:96, 1])
  
  #VlnPlot(seu,features = c('nCount_RNA','nFeature_RNA','percent_mt','S.Score'), group.by = 'scDblFinder.class', pt.size = 0)
  
  
  # convert to h5ad ---------------------------------------------------------
  source("/home/zhepan/Project/scRNA_BC_metastases/Analysis/scRNA/func/convertSeu5Format.R")
  
  if (type == 'cellranger') {
    h5names = 'cellranger_doublet.h5ad'
    rdsnames = 'cellranger_doublet.rds'
  } else if (type == 'soupx') {
    h5names = 'soupx_doublet.h5ad'
    rdsnames = 'soupx_doublet.rds'
  } else if (type == 'cellbender') {
    h5names = 'cellbender_doublet.h5ad'
    rdsnames = 'cellbender_doublet.rds'
  } else if (type == 'decontx') {
    h5names = 'decontx_doublet.h5ad'
    rdsnames = 'decontx_doublet.rds'
  } 
  
  seu = DietSeurat(seu, layers = 'counts')
  convertSeu5Format(seu, savepaths = file.path(savepath, SampleID, h5names))
  saveRDS(seu, file = file.path(savepath, SampleID, rdsnames))
  
}
