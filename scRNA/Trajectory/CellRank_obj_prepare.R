
site.lymph = c("atac","EMBOJ.GSM4909308","EMBOJ.GSM4909310","EMBOJ.GSM4909312",
               "EMBOJ.GSM4909314","EMBOJ.GSM4909316","EMBOJ.GSM4909318","EMBOJ.GSM4909321",
               "oncogenesis.B2","oncogenesis.C2","oncogenesis.D2","oncogenesis.D3","oncogenesis.E2")
site.bone = c( "BoM1","BoM2")
site.brain = c("cell.pid1","cell.pid2", "cell.pid3","GSE143423","GSE202501",
               "science.GSM4555888", "science.GSM4555889","science.GSM4555891" )
site.primary = c("primary.CID3586" ,   "primary.CID3838" ,   "primary.CID3921" ,   "primary.CID3941" ,  
                 "primary.CID3946" ,   "primary.CID3948" ,   "primary.CID3963" ,   "primary.CID4040" ,  
                 "primary.CID4066"  ,  "primary.CID4067" ,   "primary.CID4290A" ,  "primary.CID4398",   
                 "primary.CID44041" ,  "primary.CID4461" ,   "primary.CID4463",    "primary.CID4465" ,  
                 "primary.CID4471" ,   "primary.CID4495" ,   "primary.CID44971" ,  "primary.CID44991",  
                 "primary.CID4513" ,   "primary.CID4515" ,   "primary.CID45171",   "primary.CID4523" ,  
                 "primary.CID4530N" ,  "primary.CID4535")
# brain-------------------
library(Seurat)
library(harmony)

brain.list = list(cell.pid1,cell.pid2,cell.pid3,
                  GSE143423,
                  GSE202501,
                  science.GSM4555888,science.GSM4555889,science.GSM4555891)

brain.list = lapply(brain.list,
                    FUN = CellCycleScoring,
                    s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T )

brain.list = lapply(X = brain.list, FUN = SCTransform, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score","percent.mt"))

var_features_ind_sct = SelectIntegrationFeatures(brain.list, nfeatures = 3000)

brain.merge = merge(x = brain.list[[1]], y = brain.list[2:length(brain.list)])

VariableFeatures(brain.merge) = var_features_ind_sct

brain.merge = RunPCA(brain.merge)

brain.merge = RunHarmony(brain.merge,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                   kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

brain.merge = brain.merge %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

count = brain.merge@assays$SCT@counts
metad = brain.merge@meta.data
redu = brain.merge@reductions

new.brain = CreateSeuratObject(count,meta.data = metad)
new.brain@reductions = redu

saveRDS(brain.merge,file = '~/brain.rds')

library(SeuratDisk)
SaveH5Seurat(new.brain, filename = "~/brain.h5Seurat",assay = 'RNA')
Convert("~/brain.h5Seurat", dest = "h5ad")

# bone-------------------

BoM1 = CreateSeuratObject(BoM1@assays$RNA@counts,meta.data = BoM1@meta.data,min.cells = 3)
BoM2 = CreateSeuratObject(BoM2@assays$RNA@counts,meta.data = BoM2@meta.data,min.cells = 3)

bone.list = list(BoM1,BoM2)
bone.list = lapply(bone.list,
                    FUN = CellCycleScoring,
                    s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T )

bone.list = lapply(X = bone.list, FUN = SCTransform, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score","percent.mt"))

var_features_ind_sct = SelectIntegrationFeatures(bone.list, nfeatures = 3000)

bone.merge = merge(x = bone.list[[1]], y = bone.list[2:length(bone.list)])

VariableFeatures(bone.merge) = var_features_ind_sct

bone.merge = RunPCA(bone.merge)

bone.merge = RunHarmony(bone.merge,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                         kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

bone.merge = bone.merge %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

saveRDS(bone.merge,file = '~/bone.rds')

count = bone.merge@assays$SCT@counts
metad = bone.merge@meta.data
redu = bone.merge@reductions

new.bone = CreateSeuratObject(count,meta.data = metad)
new.bone@reductions = redu

library(SeuratDisk)
SaveH5Seurat(new.bone, filename = "~/bone.h5Seurat",assay = 'RNA')
Convert("~/bone.h5Seurat", dest = "h5ad")

# primary-------------------

primary.list = list(primary.CID3586, primary.CID3921, primary.CID3941,  
                  primary.CID3948,  primary.CID3963,  
                  primary.CID4066, primary.CID4067, primary.CID4290A, 
                  primary.CID44041, primary.CID4461, primary.CID4463, primary.CID4465,  
                  primary.CID4471, primary.CID4495, primary.CID44971, primary.CID44991,  
                  primary.CID4513, primary.CID4515, primary.CID45171,   primary.CID4523,  
                  primary.CID4530N, primary.CID4535)

primary.list = lapply(primary.list,
                   FUN = CellCycleScoring,
                   s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T )

primary.list = lapply(X = primary.list, FUN = SCTransform, 
                   method = "glmGamPoi",
                   vars.to.regress = c("S.Score", "G2M.Score","percent.mt"))

var_features_ind_sct = SelectIntegrationFeatures(primary.list, nfeatures = 3000)

primary.merge = merge(x = primary.list[[1]], y = primary.list[2:length(primary.list)])

VariableFeatures(primary.merge) = var_features_ind_sct

primary.merge = RunPCA(primary.merge)

primary.merge = RunHarmony(primary.merge,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                        kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

primary.merge = primary.merge %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

count = primary.merge@assays$SCT@counts
metad = primary.merge@meta.data
redu = primary.merge@reductions

new.primary = CreateSeuratObject(count,meta.data = metad)
new.primary@reductions = redu

library(SeuratDisk)
SaveH5Seurat(new.primary, filename = "~/primary.h5Seurat",assay = 'RNA')
Convert("~/primary.h5Seurat", dest = "h5ad")

# lymph-------------------

lymph.list = list(atac,EMBOJ.GSM4909308,EMBOJ.GSM4909310,EMBOJ.GSM4909312,
                  EMBOJ.GSM4909314,EMBOJ.GSM4909316,EMBOJ.GSM4909318,EMBOJ.GSM4909321,
                  oncogenesis.B2,oncogenesis.C2,oncogenesis.D2,oncogenesis.D3,oncogenesis.E2)

lymph.list = lapply(lymph.list,
                    FUN = CellCycleScoring,
                    s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T )

lymph.list = lapply(X = lymph.list, FUN = SCTransform, 
                    method = "glmGamPoi",
                    vars.to.regress = c("S.Score", "G2M.Score","percent.mt"))

var_features_ind_sct = SelectIntegrationFeatures(lymph.list, nfeatures = 3000)

lymph.merge = merge(x = lymph.list[[1]], y = lymph.list[2:length(lymph.list)])

VariableFeatures(lymph.merge) = var_features_ind_sct

lymph.merge = RunPCA(lymph.merge)

lymph.merge = RunHarmony(lymph.merge,group.by.vars = c('sample','study'),assay.use='SCT', plot_convergence = TRUE,theta = c(2.5,1.5), 
                         kmeans_init_nstart=20, kmeans_init_iter_max=5000)
print('harmonied')

lymph.merge = lymph.merge %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunUMAP(reduction = "harmony", dims = 1:30)

count = lymph.merge@assays$SCT@counts
metad = lymph.merge@meta.data
redu = lymph.merge@reductions

new.lymph = CreateSeuratObject(count,meta.data = metad)
new.lymph@reductions = redu

saveRDS(lymph.merge,file = '~/lymph.rds')

library(SeuratDisk)
SaveH5Seurat(new.lymph, filename = "~/lymph.h5Seurat",assay = 'RNA')
Convert("~/lymph.h5Seurat", dest = "h5ad")


# CD4 ---------------------------

cd4tcell.harmony <- readRDS("~/cd4tcell.harmony.rds")

count = cd4tcell.harmony@assays$RNA@counts
metad = cd4tcell.harmony@meta.data
redu = cd4tcell.harmony@reductions

new = CreateSeuratObject(count,meta.data = metad)
new@reductions = redu

library(SeuratDisk)
SaveH5Seurat(new, filename = "~/cd4tcell.harmony.h5Seurat",assay = 'RNA')
Convert("~/cd4tcell.harmony.h5Seurat", dest = "h5ad")


# cd8 ------------------------------------------------------------------
cd8tcell.harmony <- readRDS("~/cd8tcell.harmony.rds")
count = cd8tcell.harmony@assays$RNA@counts
metad = cd8tcell.harmony@meta.data
redu = cd8tcell.harmony@reductions

new = CreateSeuratObject(count,meta.data = metad)
new@reductions = redu

SaveH5Seurat(new, filename = "~/cd8tcell.harmony.h5Seurat",assay = 'RNA')
Convert("~/cd8tcell.harmony.h5Seurat", dest = "h5ad")


