import os

AUXILLIARIES_FOLDERNAME = "/home/shpc_100839/software/SCENICdata/"
RESULTS_FOLDERNAME = "/home/shpc_100839/scRNA_BC_metastases/Data/integrate_sc_BC/tumor/SCENIC/"

HUMAN_TFS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'allTFs_hg38.txt')

alloom = os.listdir(RESULTS_FOLDERNAME)

SAMPLE_ID = 'tumor'

EXP_FILTER_LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.h5ad.filter.loom'.format(SAMPLE_ID))
ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.adjacencies.tsv'.format(SAMPLE_ID))

os.system("pyscenic grn %s %s \
    -o %s \
    --num_workers 8 \
    --method grnboost2" %(EXP_FILTER_LOOM_FNAME,HUMAN_TFS_FNAME,ADJACENCIES_FNAME))
