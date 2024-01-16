#!/home/zhepan/miniconda3/envs/cnmf_env/bin/python

import os
import numpy as np
from cnmf import cNMF
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path')
parser.add_argument('--site')
parser.add_argument('--sample_list', nargs='+')

opts = parser.parse_args()

np.random.seed(42)

path = opts.path
site = opts.site
sample_list = opts.sample_list

for sample in sample_list:
    numiter = 50  # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
    numhvgenes = 7000  # Number of over-dispersed genes to use for running the actual factorizations

    ## Results will be saved to [output_directory]/[run_name] which in this example is example_PBMC/cNMF/pbmc_cNMF
    output_directory = path + '/' + site + '/' + sample
    run_name = 'cNMF'

    ## Specify the Ks to use as a space separated list in this case "5 6 7 8 9 10"
    K = np.arange(3,10)

    ## To speed this up, you can run it for only K=7-8 with the option below
    #K = ' '.join([str(i) for i in range(7,9)])
    run_seed = 42 ## Specify a seed pseudorandom number generation for reproducibility

    ## Path to the filtered counts dataset we output previously
    countfn = path+ '/' +site+ '/' +sample+ '/' +'count.txt'
    ## Initialize the cnmf object that will be used to run analyses
    cnmf_obj = cNMF(output_dir=output_directory, name=run_name)
    ## Prepare the data, I.e. subset to 2000 high-variance genes, and variance normalize
    cnmf_obj.prepare(counts_fn=countfn, components=K, n_iter=numiter, seed=run_seed, num_highvar_genes=numhvgenes)
    cnmf_obj.factorize_multi_process(2)
    cnmf_obj.combine()
    cnmf_obj.k_selection_plot(close_fig=False)
    density_threshold = 2.00
    for selected_K in K:
        cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
    density_threshold = 0.02
    for selected_K in K:
        cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)