cyclegenes = read.delim('/home/zhepan/Reference/regev_lab_cell_cycle_genes.txt')

seu = Patient_E_doublet

savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/public/raw_sc_BC/BrM/science'
SampleID = 'BrM3'
  
seu$SampleID = SampleID

seu$percent_mt = PercentageFeatureSet(seu, pattern = '^MT-')
seu$percent_hb = PercentageFeatureSet(seu, pattern = "^HB[^(P)]")
seu$percent_rb = PercentageFeatureSet(seu, pattern = "^RP[SL]")

seu <- NormalizeData(seu)

seu = CellCycleScoring(seu, s.features = cyclegenes[1:42, 1], g2m.features = cyclegenes[43:96, 1])

#VlnPlot(seu,features = c('nCount_RNA','nFeature_RNA','percent_mt','S.Score'), group.by = 'scDblFinder.class', pt.size = 0)


# convert to h5ad ---------------------------------------------------------
source("/home/zhepan/Project/scRNA_BC_metastases/Analysis/scRNA/func/convertSeu5Format.R")

h5names = 'cellranger_doublet.h5ad'
rdsnames = 'cellranger_doublet.rds'

seu = DietSeurat(seu, layers = 'counts')
convertSeu5Format(seu, savepaths = file.path(savepath,SampleID, h5names))
saveRDS(seu, file = file.path(savepath,SampleID, rdsnames))
