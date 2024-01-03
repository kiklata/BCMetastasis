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


dataPath = '/home/zhepan/Project/MultiOmics/data/skin/res/'
RESULTS_FOLDERNAME = "/home/zhepan/Project/MultiOmics/data/skin/res/SCENIC"

SAMPLE_ID = 'cellbender_Keratinocyte_count'

EXP_FILTER_LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.loom'.format(SAMPLE_ID))
expath = os.path.join(dataPath,'{}.h5ad'.format(SAMPLE_ID))

#CreateLoom(h5adpath = expath, loompath = EXP_FILTER_LOOM_FNAME)

# step 1

AUXILLIARIES_FOLDERNAME = "/home/zhepan/Reference/SCENIC/"

RANKING_DBS_FNAMES = os.path.join(AUXILLIARIES_FOLDERNAME,'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather')
MOTIF_ANNOTATIONS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')
HUMAN_TFS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'allTFs_hg38.txt')

ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.adjacencies.tsv'.format(SAMPLE_ID))

os.system("pyscenic grn %s %s \
    -o %s \
    --num_workers 8 \
    --method grnboost2" %(EXP_FILTER_LOOM_FNAME,HUMAN_TFS_FNAME,ADJACENCIES_FNAME))


# step 2

DBS_PARAM = RANKING_DBS_FNAMES

MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.motifs.csv'.format(SAMPLE_ID))

os.system("pyscenic ctx %s %s \
    --annotations_fname %s \
    --expression_mtx_fname %s \
    --output %s \
    --num_workers 4" %(ADJACENCIES_FNAME,DBS_PARAM,MOTIF_ANNOTATIONS_FNAME,EXP_FILTER_LOOM_FNAME,MOTIFS_FNAME))


# step 4
AUCELL_LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.aucell.csv'.format(SAMPLE_ID))


os.system("pyscenic aucell %s %s \
    --output %s \
    --num_workers 8" %(EXP_FILTER_LOOM_FNAME,MOTIFS_FNAME,AUCELL_LOOM_FNAME))
