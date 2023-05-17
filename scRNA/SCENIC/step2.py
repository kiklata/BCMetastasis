import os

AUXILLIARIES_FOLDERNAME = "/home/shpc_100839/software/SCENICdata/"
RESULTS_FOLDERNAME = "/home/shpc_100839/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/SCENIC/"
RANKING_DBS_FNAMES = os.path.join(AUXILLIARIES_FOLDERNAME,'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather')
MOTIF_ANNOTATIONS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')


HUMAN_TFS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'allTFs_hg38.txt')

DBS_PARAM = RANKING_DBS_FNAMES

SAMPLE_ID = 'tumor'

EXP_FILTER_LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.h5ad.filter.loom'.format(SAMPLE_ID))
ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.adjacencies.tsv'.format(SAMPLE_ID))
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.motifs.csv'.format(SAMPLE_ID))

os.system("pyscenic ctx %s %s \
    --annotations_fname %s \
    --expression_mtx_fname %s \
    --output %s \
    --num_workers 4" %(ADJACENCIES_FNAME,DBS_PARAM,MOTIF_ANNOTATIONS_FNAME,EXP_FILTER_LOOM_FNAME,MOTIFS_FNAME))

