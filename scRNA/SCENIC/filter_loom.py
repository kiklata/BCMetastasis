
import os
import numpy as np
import scanpy as sc
import loompy as lp

def CreateLoom(h5adpath,loompath):

    EXP_H5AD_FNAME = h5adpath
    EXP_FILTER_LOOM_FNAME = loompath

    adata = sc.read_h5ad(EXP_H5AD_FNAME)
    row_attrs = {
        "Gene": np.array(adata.var_names) ,
    }
    col_attrs = {
        "CellID": np.array(adata.obs_names) ,
        "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
        "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
    }
    lp.create( EXP_FILTER_LOOM_FNAME, adata.X.transpose(), row_attrs, col_attrs)


TumorPath = '/home/shpc_100839/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/h5ad/filter_NG/'
SCENICpath = '/home/shpc_100839/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/SCENIC/'

file = 'tumor.h5ad'
expath = os.path.join(TumorPath,file)
savepath = os.path.join(SCENICpath,'{}.filter.loom'.format(file))
print(savepath)
CreateLoom(h5adpath = expath, loompath = savepath)
