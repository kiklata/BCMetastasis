tcell = c('CRTAM+ gamma-delta T cells','Cycling gamma-delta T cells','Cycling T cells','Follicular helper T cells','gamma-delta T cells',
          'MAIT cells','Memory CD4+ cytotoxic T cells','Regulatory T cells','T(agonist)','Tcm/Naive cytotoxic T cells','Tcm/Naive helper T cells',
          'Tem/Effector helper T cells','Tem/Effector helper T cells PD1+','Tem/Temra cytotoxic T cells','Tem/Trm cytotoxic T cells','Treg(diff)',
          'Trm cytotoxic T cells','Type 1 helper T cells','Type 17 helper T cells','CD8a/b(entry)','CD8a/a')

bcell = c('Age-associated B cells','B cells','Cycling B cells','Follicular B cells','Germinal center B cells','Large pre-B cells',
         'Memory B cells','Naive B cells','Plasma cells','Plasmablasts','Pro-B cells','Proliferative germinal center B cells',
         'Small pre-B cells','Transitional B cells')

mono = c('Alveolar macrophages','Classical monocytes','Erythrophagocytic macrophages','Intermediate macrophages',
         'Intestinal macrophages','Kupffer cells','Macrophages','Mono-mac','Monocyte precursor','Monocytes','Non-classical monocytes')


ilc = c('CD16- NK cells','CD16+ NK cells','Cycling NK cells','ILC','ILC precursor','ILC3','NK cells','NKT cells','Transitional NK')

early = c('CMP','Double-negative thymocytes','Double-positive thymocytes','Early erythroid','Early lymphoid/T lymphoid',
          'ELP','Endothelial cells','Epithelial cells','Erythrocytes','ETP','Fibroblasts','Hofbauer cells','HSC/MPP','Late erythroid',
          'Megakaryocytes/platelets','MEMP','Mid erythroid','MNP')

dc = c('Cycling DCs','DC','DC precursor','DC1','DC2','DC3','Migratory DCs','pDC','pDC precursor','Transitional DC')

myeloid = c('Granulocytes','Mast cells','Myelocytes','Neutrophil-myeloid progenitor')

immune.celltypist.low$celltype.major = 
  if_else(immune.celltypist.low$celltypist.low %in% tcell,'T cell',
  if_else(immune.celltypist.low$celltypist.low %in% bcell,'B cell',
  if_else(immune.celltypist.low$celltypist.low %in% mono,'Monocyte/Macrophage',
  if_else(immune.celltypist.low$celltypist.low %in% ilc, 'ILC',
  if_else(immune.celltypist.low$celltypist.low %in% early ,'Early-developed or Proliferative cell',
  if_else(immune.celltypist.low$celltypist.low %in% dc, 'DC',
  if_else(immune.celltypist.low$celltypist.low %in% myeloid,'Myelocyte','na')))))))

t.helper = c('Follicular helper T cells','Tcm/Naive helper T cells','Tem/Effector helper T cells','Tem/Effector helper T cells PD1+',
             'Type 1 helper T cells','Type 17 helper T cells')

t.cyto = c('Memory CD4+ cytotoxic T cells','Tcm/Naive cytotoxic T cells','Tem/Temra cytotoxic T cells','Tem/Trm cytotoxic T cells',
           'Trm cytotoxic T cells')

t.reg = c('Regulatory T cells','Treg(diff)')

t.gd = c('CRTAM+ gamma-delta T cells','Cycling gamma-delta T cells','gamma-delta T cells')

t.cycle = c('Cycling T cells','CD8a/b(entry)','CD8a/a')

b.pre = c('Large pre-B cells','Pro-B cells','Small pre-B cells','Transitional B cells','Cycling B cells')

b.mature = c('Age-associated B cells','B cells','Follicular B cells','Germinal center B cells',
             'Memory B cells','Naive B cells','Proliferative germinal center B cells')

b.plasm = c('Plasma cells','Plasmablasts')

monocyte = c('Classical monocytes','Mono-mac','Monocyte precursor','Monocytes','Non-classical monocytes')

macro = c('Alveolar macrophages','Erythrophagocytic macrophages','Intermediate macrophages','Intestinal macrophages',
          'Kupffer cells','Macrophages') 

nk = c('CD16- NK cells','CD16+ NK cells','Cycling NK cells','NK cells','NKT cells','Transitional NK')

ilc.tran = c('ILC','ILC precursor','ILC3')

dc.pre = c('Cycling DCs','DC precursor','Migratory DCs','pDC precursor','Transitional DC')

dc.mature = c('DC','DC1','DC2','DC3','pDC')

immune.celltypist.low$celltype.minor = 
  if_else(immune.celltypist.low$celltypist.low %in% t.cycle,'cycling of proliferative T cell',
  if_else(immune.celltypist.low$celltypist.low %in% t.helper,'helper T cell',
  if_else(immune.celltypist.low$celltypist.low %in% t.cyto, 'cytotoxic T cell',
  if_else(immune.celltypist.low$celltypist.low %in% t.reg, 'regulatory T cell',
  if_else(immune.celltypist.low$celltypist.low %in% t.gd, 'gamma-delta T cell',
  if_else(immune.celltypist.low$celltypist.low %in% b.pre, 'cycling of proliferative B cell',
  if_else(immune.celltypist.low$celltypist.low %in% b.mature, 'mature B cell',
  if_else(immune.celltypist.low$celltypist.low %in% b.plasm, 'Plasma cell',
  if_else(immune.celltypist.low$celltypist.low %in% monocyte ,'Monocyte',
  if_else(immune.celltypist.low$celltypist.low %in% macro, 'Macrophage',
  if_else(immune.celltypist.low$celltypist.low %in% nk,'NK cell',
  if_else(immune.celltypist.low$celltypist.low %in% ilc.tran, 'ILC',
  if_else(immune.celltypist.low$celltypist.low %in% dc.pre,'cycling of proliferative DC',
  if_else(immune.celltypist.low$celltypist.low %in% dc.mature,'mature DC',
  if_else(immune.celltypist.low$celltypist.low %in% early,'Early-developed or Proliferative cell',
          immune.celltypist.low$celltypist.low)))))))))))))))
  
immune.celltypist.low$celltype.subset = immune.celltypist.low$celltypist.low


early = c('CMP','Double-negative thymocytes','Double-positive thymocytes','Early erythroid','Early lymphoid/T lymphoid',
          'ELP','Endothelial cells','Epithelial cells','Erythrocytes','ETP','Fibroblasts','Hofbauer cells','HSC/MPP','Late erythroid',
          'Megakaryocytes/platelets','MEMP','Mid erythroid','MNP')

immune.celltypist.low$celltype.subset = if_else(immune.celltypist.low$celltypist.low %in% c('Endothelial cells','Epithelial cells','Fibroblasts'),'Proliferative cell',
                                                immune.celltypist.low$celltype.subset)

immune.anno = immune.celltypist.low[,c(1,3:5)]
