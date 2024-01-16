cyclegenes = read.delim('/home/zhepan/Reference/regev_lab_cell_cycle_genes.txt')
source("/home/zhepan/Project/scRNA_BC_metastases/Analysis/scRNA/func/convertSeu5Format.R")
library(scDblFinder)

h5names = 'cellranger_doublet.h5ad'
rdsnames = 'cellranger_doublet.rds'
set.seed(42)

meta = read.delim('~/Project/scRNA_BC_metastases/Data/SingleCell/public/LymphM/NC/meta_all.csv.gz',sep = ',')
sample_list = c('LymphM1','LymphM2','LymphM3','LymphM4','LymphM5')

savepath = '/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/public/LymphM/NC'

for (SampleID in sample_list) {
  
count = Read10X(paste0(savepath,'/',SampleID))
#count = data.table::fread(paste0(savepath,'/',SampleID,'/expression_matrix.txt.gz')) %>% as.data.frame()
#count = tibble::column_to_rownames(count, 'V1')

seu = CreateSeuratObject(count)

  
seu$SampleID = SampleID

seu$percent_mt = PercentageFeatureSet(seu, pattern = '^MT-')
seu$percent_hb = PercentageFeatureSet(seu, pattern = "^HB[^(P)]")
seu$percent_rb = PercentageFeatureSet(seu, pattern = "^RP[SL]")

seu <- NormalizeData(seu)

seu = CellCycleScoring(seu, s.features = cyclegenes[1:42, 1], g2m.features = cyclegenes[43:96, 1])

sce = scDblFinder(as.SingleCellExperiment(seu))
seu$scDblFinder.class = sce$scDblFinder.class

# convert to h5ad ---------------------------------------------------------

seu = DietSeurat(seu, layers = 'counts')
convertSeu5Format(seu, savepaths = file.path(savepath,SampleID, h5names))
saveRDS(seu, file = file.path(savepath,SampleID, rdsnames))

}
