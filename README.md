``` r
> sessionInfo()
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] SCP_0.1.10                  survminer_0.4.9             survival_3.2-13            
 [4] ggpubr_0.4.0                clusterProfiler_4.4.4       GSEABase_1.58.0            
 [7] graph_1.74.0                annotate_1.74.0             XML_3.99-0.10              
[10] AnnotationDbi_1.58.0        GSVA_1.44.2                 BiocParallel_1.30.3        
[13] SCENIC_1.3.1                SeuratDisk_0.0.0.9020       patchwork_1.1.1            
[16] harmony_0.1.0               Rcpp_1.0.9                  SeuratWrappers_0.3.0       
[19] monocle3_1.2.9              SingleCellExperiment_1.18.0 SummarizedExperiment_1.26.1
[22] GenomicRanges_1.48.0        GenomeInfoDb_1.32.2         IRanges_2.30.0             
[25] S4Vectors_0.34.0            MatrixGenerics_1.8.1        matrixStats_0.62.0         
[28] Biobase_2.56.0              BiocGenerics_0.42.0         Seurat_4.1.1               
[31] ggsci_2.9                   viridis_0.6.2               viridisLite_0.4.0          
[34] scibet_1.0                  forcats_0.5.1               stringr_1.4.0              
[37] dplyr_1.0.9                 purrr_0.3.4                 readr_2.1.2                
[40] tidyr_1.2.0                 tibble_3.1.8                tidyverse_1.3.2            
[43] ggplot2_3.3.6               sp_1.5-0                    SeuratObject_4.1.0         

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                ica_1.0-3                 foreach_1.5.2             lmtest_0.9-40            
  [5] crayon_1.5.1              spatstat.core_2.4-4       MASS_7.3-58               rhdf5filters_1.8.0       
  [9] nlme_3.1-157              backports_1.4.1           reprex_2.0.1              GOSemSim_2.22.0          
 [13] rlang_1.0.4               XVector_0.36.0            ROCR_1.0-11               readxl_1.4.0             
 [17] irlba_2.3.5               nloptr_2.0.3              filelock_1.0.2            rjson_0.2.21             
 [21] bit64_4.0.5               glue_1.6.2                sctransform_0.3.3         parallel_4.2.1           
 [25] spatstat.sparse_2.1-1     DOSE_3.22.0               spatstat.geom_2.4-0       haven_2.5.0              
 [29] tidyselect_1.1.2          km.ci_0.5-6               fitdistrplus_1.1-8        zoo_1.8-10               
 [33] xtable_1.8-4              magrittr_2.0.3            cli_3.3.0                 zlibbioc_1.42.0          
 [37] rstudioapi_0.13           miniUI_0.1.1.1            parallelDist_0.2.6        rpart_4.1.16             
 [41] fastmatch_1.1-3           treeio_1.20.1             shiny_1.7.2               BiocSingular_1.12.0      
 [45] xfun_0.31                 clue_0.3-61               cluster_2.1.3             tidygraph_1.2.1          
 [49] KEGGREST_1.36.3           ggrepel_0.9.1             ape_5.6-2                 listenv_0.8.0            
 [53] Biostrings_2.64.0         png_0.1-7                 future_1.27.0             withr_2.5.0              
 [57] bitops_1.0-7              ggforce_0.3.3             plyr_1.8.7                cellranger_1.1.0         
 [61] RcppParallel_5.1.5        pillar_1.8.0              GlobalOptions_0.1.2       cachem_1.0.6             
 [65] fs_1.5.2                  hdf5r_1.3.5               GetoptLong_1.0.5          RcppML_0.3.7             
 [69] DelayedMatrixStats_1.18.0 vctrs_0.4.1               ellipsis_0.3.2            generics_0.1.3           
 [73] tools_4.2.1               munsell_0.5.0             tweenr_1.0.2              fgsea_1.22.0             
 [77] DelayedArray_0.22.0       fastmap_1.1.0             compiler_4.2.1            abind_1.4-5              
 [81] httpuv_1.6.5              plotly_4.10.0             rgeos_0.5-9               GenomeInfoDbData_1.2.8   
 [85] gridExtra_2.3             ggnewscale_0.4.7          lattice_0.20-45           deldir_1.0-6             
 [89] utf8_1.2.2                later_1.3.0               BiocFileCache_2.4.0       jsonlite_1.8.0           
 [93] princurve_2.1.6           scales_1.2.0              ScaledMatrix_1.4.0        tidytree_0.3.9           
 [97] pbapply_1.5-0             carData_3.0-5             sparseMatrixStats_1.8.0   lazyeval_0.2.2           
[101] promises_1.2.0.1          car_3.1-0                 doParallel_1.0.17         R.utils_2.11.0           
[105] goftest_1.2-3             spatstat.utils_2.3-1      reticulate_1.25           cowplot_1.1.1            
[109] Rtsne_0.16                downloader_0.4            uwot_0.1.11               igraph_1.3.4             
[113] proxyC_0.2.4              HDF5Array_1.24.1          htmltools_0.5.3           memoise_2.0.1            
[117] graphlayouts_0.8.0        digest_0.6.29             assertthat_0.2.1          mime_0.12                
[121] rappdirs_0.3.3            KMsurv_0.1-5              RSQLite_2.2.15            yulab.utils_0.0.5        
[125] future.apply_1.9.0        remotes_2.4.2             data.table_1.14.2         blob_1.2.3               
[129] R.oo_1.25.0               survMisc_0.5.6            splines_4.2.1             Rhdf5lib_1.18.2          
[133] googledrive_2.0.0         RCurl_1.98-1.8            broom_1.0.0               hms_1.1.1                
[137] modelr_0.1.8              rhdf5_2.40.0              colorspace_2.0-3          BiocManager_1.30.18      
[141] shape_1.4.6               aplot_0.1.4               RANN_2.6.1                circlize_0.4.15          
[145] enrichplot_1.16.1         fansi_1.0.3               tzdb_0.3.0                parallelly_1.32.1        
[149] R6_2.5.1                  grid_4.2.1                ggridges_0.5.3            lifecycle_1.0.1          
[153] curl_4.3.2                ggsignif_0.6.3            googlesheets4_1.0.0       minqa_1.2.4              
[157] leiden_0.4.2              DO.db_2.9                 Matrix_1.5-0              qvalue_2.28.0            
[161] RcppAnnoy_0.0.19          RColorBrewer_1.1-3        iterators_1.0.14          R.cache_0.15.0           
[165] htmlwidgets_1.5.4         beachmat_2.12.0           polyclip_1.10-0           biomaRt_2.52.0           
[169] shadowtext_0.1.2          gridGraphics_0.5-1        terra_1.6-7               rvest_1.0.2              
[173] ComplexHeatmap_2.13.2     mgcv_1.8-40               globals_0.15.1            spatstat.random_2.2-0    
[177] slingshot_2.4.0           progressr_0.10.1          codetools_0.2-18          lubridate_1.8.0          
[181] GO.db_3.15.0              prettyunits_1.1.1         dbplyr_2.2.1              R.methodsS3_1.8.2        
[185] gtable_0.3.0              DBI_1.1.3                 ggfun_0.0.6               tensor_1.5               
[189] httr_1.4.3                KernSmooth_2.23-20        stringi_1.7.8             progress_1.2.2           
[193] reshape2_1.4.4            farver_2.1.1              ggtree_3.4.1              xml2_1.3.3               
[197] boot_1.3-28               AUCell_1.18.1             lme4_1.1-30               ggplotify_0.1.0          
[201] scattermore_0.8           bit_4.0.4                 scatterpie_0.1.7          spatstat.data_2.2-0      
[205] ggraph_2.0.5              TrajectoryUtils_1.4.0     pkgconfig_2.0.3           gargle_1.2.0             
[209] rstatix_0.7.0             knitr_1.39   
```

```shell
root:~$ conda list -n Renv

# packages in environment at /home/shpc_100839/miniconda3/envs/Renv:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
anndata                   0.8.0                    pypi_0    pypi
annoy                     1.17.1                   pypi_0    pypi
arpack                    3.7.0                hdefa2d7_2    conda-forge
bypy                      1.8                      pypi_0    pypi
bzip2                     1.0.8                h7f98852_4    conda-forge
ca-certificates           2022.6.15            ha878542_0    conda-forge
certifi                   2022.6.15       py310hff52083_0    conda-forge
charset-normalizer        2.1.0                    pypi_0    pypi
cnmf                      1.3.4                    pypi_0    pypi
cycler                    0.11.0                   pypi_0    pypi
dill                      0.3.5.1                  pypi_0    pypi
fastcluster               1.2.6                    pypi_0    pypi
fbpca                     1.0                      pypi_0    pypi
fonttools                 4.37.2                   pypi_0    pypi
geosketch                 1.2                      pypi_0    pypi
glpk                      4.65              h9202a9a_1004    conda-forge
gmp                       6.2.1                h58526e2_0    conda-forge
h5py                      3.7.0                    pypi_0    pypi
icu                       70.1                 h27087fc_0    conda-forge
idna                      3.3                      pypi_0    pypi
igraph                    0.9.9                h026ac8f_0    conda-forge
intervaltree              2.1.0                    pypi_0    pypi
jags                      4.3.0             h236a147_1004    conda-forge
joblib                    1.1.0                    pypi_0    pypi
kiwisolver                1.4.4                    pypi_0    pypi
ld_impl_linux-64          2.36.1               hea4e1c9_2    conda-forge
leidenalg                 0.8.10          py310hd8f1fbe_0    conda-forge
libblas                   3.9.0           16_linux64_openblas    conda-forge
libcblas                  3.9.0           16_linux64_openblas    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc-ng                 12.1.0              h8d9b700_16    conda-forge
libgfortran-ng            12.1.0              h69a702a_16    conda-forge
libgfortran5              12.1.0              hdcd56e2_16    conda-forge
libgomp                   12.1.0              h8d9b700_16    conda-forge
libiconv                  1.16                 h516909a_0    conda-forge
liblapack                 3.9.0           16_linux64_openblas    conda-forge
libnsl                    2.0.0                h7f98852_0    conda-forge
libopenblas               0.3.21          pthreads_h78a6416_0    conda-forge
libsqlite                 3.39.2               h753d276_1    conda-forge
libstdcxx-ng              12.1.0              ha89aaad_16    conda-forge
libuuid                   2.32.1            h7f98852_1000    conda-forge
libxml2                   2.9.14               h22db469_4    conda-forge
libzlib                   1.2.12               h166bdaf_2    conda-forge
llvmlite                  0.39.0                   pypi_0    pypi
matplotlib                3.5.3                    pypi_0    pypi
metis                     5.1.0             h58526e2_1006    conda-forge
mpfr                      4.1.0                h9202a9a_1    conda-forge
multiprocess              0.70.13                  pypi_0    pypi
natsort                   8.2.0                    pypi_0    pypi
ncurses                   6.3                  h27087fc_1    conda-forge
networkx                  2.8.6                    pypi_0    pypi
numba                     0.56.0                   pypi_0    pypi
numpy                     1.22.4                   pypi_0    pypi
openssl                   3.0.5                h166bdaf_1    conda-forge
packaging                 21.3                     pypi_0    pypi
palettable                3.3.0                    pypi_0    pypi
pandas                    1.4.4                    pypi_0    pypi
patsy                     0.5.2                    pypi_0    pypi
pillow                    9.2.0                    pypi_0    pypi
pip                       22.2.2             pyhd8ed1ab_0    conda-forge
pkg-config                0.29.2            h36c2ea0_1008    conda-forge
pynndescent               0.5.7                    pypi_0    pypi
pyparsing                 3.0.9                    pypi_0    pypi
python                    3.10.5          ha86cf86_0_cpython    conda-forge
python-dateutil           2.8.2                    pypi_0    pypi
python-igraph             0.9.11          py310hb58d747_0    conda-forge
python_abi                3.10                    2_cp310    conda-forge
pytz                      2022.2.1                 pypi_0    pypi
pyyaml                    6.0                      pypi_0    pypi
readline                  8.1.2                h0f457ee_0    conda-forge
requests                  2.28.1                   pypi_0    pypi
requests-toolbelt         0.9.1                    pypi_0    pypi
scanoramact               1.2.0                    pypi_0    pypi
scanpy                    1.9.1                    pypi_0    pypi
scikit-learn              1.1.2                    pypi_0    pypi
scipy                     1.9.0                    pypi_0    pypi
seaborn                   0.12.0                   pypi_0    pypi
session-info              1.0.0                    pypi_0    pypi
setuptools                65.0.0          py310hff52083_0    conda-forge
six                       1.16.0                   pypi_0    pypi
some-package              0.1                      pypi_0    pypi
sortedcontainers          2.4.0                    pypi_0    pypi
sqlite                    3.39.2               h4ff8645_1    conda-forge
statsmodels               0.13.2                   pypi_0    pypi
stdlib-list               0.8.0                    pypi_0    pypi
suitesparse               5.10.1               h9e50725_1    conda-forge
tbb                       2021.5.0             h924138e_1    conda-forge
texttable                 1.6.4              pyhd8ed1ab_0    conda-forge
threadpoolctl             3.1.0                    pypi_0    pypi
tk                        8.6.12               h27826a3_0    conda-forge
tqdm                      4.64.0                   pypi_0    pypi
tzdata                    2022b                h191b570_0    conda-forge
umap-learn                0.5.3                    pypi_0    pypi
urllib3                   1.26.11                  pypi_0    pypi
wheel                     0.37.1             pyhd8ed1ab_0    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
```

```shell
root:~$ conda list -n pyscenic

# packages in environment at /home/shpc_100839/miniconda3/envs/pyscenic:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main    defaults
_openmp_mutex             5.1                       1_gnu    defaults
aiohttp                   3.8.1                    pypi_0    pypi
aiosignal                 1.2.0                    pypi_0    pypi
anndata                   0.8.0                    pypi_0    pypi
arboreto                  0.1.6                    pypi_0    pypi
asttokens                 2.0.5              pyhd3eb1b0_0    defaults
async-timeout             4.0.2                    pypi_0    pypi
attrs                     22.1.0                   pypi_0    pypi
backcall                  0.2.0              pyhd3eb1b0_0    defaults
blas                      1.0                         mkl    defaults
bokeh                     2.4.3                    pypi_0    pypi
boltons                   21.0.0                   pypi_0    pypi
bottleneck                1.3.5           py310ha9d4c09_0    defaults
brotli                    1.0.9                h5eee18b_7    defaults
brotli-bin                1.0.9                h5eee18b_7    defaults
bzip2                     1.0.8                h7b6447c_0    defaults
ca-certificates           2022.07.19           h06a4308_0    defaults
certifi                   2022.6.15       py310h06a4308_0    defaults
charset-normalizer        2.1.1                    pypi_0    pypi
click                     8.1.3                    pypi_0    pypi
cloudpickle               2.2.0                    pypi_0    pypi
ctxcore                   0.2.0                    pypi_0    pypi
cycler                    0.11.0             pyhd3eb1b0_0    defaults
cytoolz                   0.12.0                   pypi_0    pypi
dask                      2022.9.1                 pypi_0    pypi
dbus                      1.13.18              hb2f20db_0    defaults
debugpy                   1.5.1           py310h295c915_0    defaults
decorator                 5.1.1              pyhd3eb1b0_0    defaults
dill                      0.3.5.1                  pypi_0    pypi
distributed               2022.9.1                 pypi_0    pypi
entrypoints               0.4             py310h06a4308_0    defaults
executing                 0.8.3              pyhd3eb1b0_0    defaults
expat                     2.4.4                h295c915_0    defaults
fftw                      3.3.9                h27cfd23_1    defaults
fontconfig                2.13.1               h6c09931_0    defaults
fonttools                 4.25.0             pyhd3eb1b0_0    defaults
freetype                  2.11.0               h70c0345_0    defaults
frozendict                2.3.4                    pypi_0    pypi
frozenlist                1.3.1                    pypi_0    pypi
fsspec                    2022.8.2                 pypi_0    pypi
giflib                    5.2.1                h7b6447c_0    defaults
glib                      2.69.1               h4ff587b_1    defaults
gst-plugins-base          1.14.0               h8213a91_2    defaults
gstreamer                 1.14.0               h28cd5cc_2    defaults
h5py                      3.7.0                    pypi_0    pypi
heapdict                  1.0.1                    pypi_0    pypi
icu                       58.2                 he6710b0_3    defaults
idna                      3.4                      pypi_0    pypi
intel-openmp              2021.4.0          h06a4308_3561    defaults
interlap                  0.2.7                    pypi_0    pypi
ipykernel                 6.15.2          py310h06a4308_0    defaults
ipython                   8.4.0           py310h06a4308_0    defaults
jedi                      0.18.1          py310h06a4308_1    defaults
jinja2                    3.1.2                    pypi_0    pypi
joblib                    1.2.0                    pypi_0    pypi
jpeg                      9e                   h7f8727e_0    defaults
jupyter_client            7.3.5           py310h06a4308_0    defaults
jupyter_core              4.10.0          py310h06a4308_0    defaults
kiwisolver                1.4.2           py310h295c915_0    defaults
krb5                      1.19.2               hac12032_0    defaults
lcms2                     2.12                 h3be6417_0    defaults
ld_impl_linux-64          2.38                 h1181459_1    defaults
lerc                      3.0                  h295c915_0    defaults
libbrotlicommon           1.0.9                h5eee18b_7    defaults
libbrotlidec              1.0.9                h5eee18b_7    defaults
libbrotlienc              1.0.9                h5eee18b_7    defaults
libclang                  10.0.1          default_hb85057a_2    defaults
libdeflate                1.8                  h7f8727e_5    defaults
libedit                   3.1.20210910         h7f8727e_0    defaults
libevent                  2.1.12               h8f2d780_0    defaults
libffi                    3.3                  he6710b0_2    defaults
libgcc-ng                 11.2.0               h1234567_1    defaults
libgfortran-ng            11.2.0               h00389a5_1    defaults
libgfortran5              11.2.0               h1234567_1    defaults
libgomp                   11.2.0               h1234567_1    defaults
libllvm10                 10.0.1               hbcb73fb_5    defaults
libpng                    1.6.37               hbc83047_0    defaults
libpq                     12.9                 h16c4e8d_3    defaults
libsodium                 1.0.18               h7b6447c_0    defaults
libstdcxx-ng              11.2.0               h1234567_1    defaults
libtiff                   4.4.0                hecacb30_0    defaults
libuuid                   1.0.3                h7f8727e_2    defaults
libwebp                   1.2.2                h55f646e_0    defaults
libwebp-base              1.2.2                h7f8727e_0    defaults
libxcb                    1.15                 h7f8727e_0    defaults
libxkbcommon              1.0.1                hfa300c1_0    defaults
libxml2                   2.9.14               h74e7548_0    defaults
libxslt                   1.1.35               h4e12654_0    defaults
llvmlite                  0.39.1                   pypi_0    pypi
locket                    1.0.0                    pypi_0    pypi
loompy                    3.0.7                    pypi_0    pypi
lz4-c                     1.9.3                h295c915_1    defaults
markupsafe                2.1.1                    pypi_0    pypi
matplotlib                3.5.2           py310h06a4308_0    defaults
matplotlib-base           3.5.2           py310hf590b9c_0    defaults
matplotlib-inline         0.1.6           py310h06a4308_0    defaults
mkl                       2021.4.0           h06a4308_640    defaults
mkl-service               2.4.0           py310h7f8727e_0    defaults
mkl_fft                   1.3.1           py310hd6ae3a3_0    defaults
mkl_random                1.2.2           py310h00e6091_0    defaults
msgpack                   1.0.4                    pypi_0    pypi
multidict                 6.0.2                    pypi_0    pypi
multiprocessing-on-dill   3.5.0a4                  pypi_0    pypi
munkres                   1.1.4                      py_0    defaults
natsort                   8.2.0                    pypi_0    pypi
ncurses                   6.3                  h5eee18b_3    defaults
nest-asyncio              1.5.5           py310h06a4308_0    defaults
networkx                  2.8.6                    pypi_0    pypi
nspr                      4.33                 h295c915_0    defaults
nss                       3.74                 h0370c37_0    defaults
numba                     0.56.2                   pypi_0    pypi
numexpr                   2.8.3           py310hcea2de6_0    defaults
numpy                     1.22.4                   pypi_0    pypi
numpy-base                1.21.5          py310hcba007f_3    defaults
numpy-groupies            0.9.19                   pypi_0    pypi
openssl                   1.1.1q               h7f8727e_0    defaults
packaging                 21.3               pyhd3eb1b0_0    defaults
pandas                    1.4.4                    pypi_0    pypi
parso                     0.8.3              pyhd3eb1b0_0    defaults
partd                     1.3.0                    pypi_0    pypi
patsy                     0.5.2                    pypi_0    pypi
pcre                      8.45                 h295c915_0    defaults
pexpect                   4.8.0              pyhd3eb1b0_3    defaults
pickleshare               0.7.5           pyhd3eb1b0_1003    defaults
pillow                    9.2.0                    pypi_0    pypi
pip                       22.1.2          py310h06a4308_0    defaults
ply                       3.11            py310h06a4308_0    defaults
prompt-toolkit            3.0.20             pyhd3eb1b0_0    defaults
psutil                    5.9.2                    pypi_0    pypi
ptyprocess                0.7.0              pyhd3eb1b0_2    defaults
pure_eval                 0.2.2              pyhd3eb1b0_0    defaults
pyarrow                   9.0.0                    pypi_0    pypi
pygments                  2.11.2             pyhd3eb1b0_0    defaults
pynndescent               0.5.7                    pypi_0    pypi
pyparsing                 3.0.9           py310h06a4308_0    defaults
pyqt                      5.15.7          py310h6a678d5_1    defaults
pyqt5-sip                 12.11.0                  pypi_0    pypi
pyscenic                  0.12.0                   pypi_0    pypi
python                    3.10.4               h12debd9_0    defaults
python-dateutil           2.8.2              pyhd3eb1b0_0    defaults
pytz                      2022.2.1                 pypi_0    pypi
pyyaml                    6.0                      pypi_0    pypi
pyzmq                     23.2.0          py310h6a678d5_0    defaults
qt-main                   5.15.2               h327a75a_7    defaults
qt-webengine              5.15.9               hd2b0992_4    defaults
qtwebkit                  5.212                h4eab89a_4    defaults
readline                  8.1.2                h7f8727e_1    defaults
requests                  2.28.1                   pypi_0    pypi
scanpy                    1.9.1                    pypi_0    pypi
scikit-learn              1.1.2                    pypi_0    pypi
scipy                     1.9.1                    pypi_0    pypi
seaborn                   0.11.2             pyhd3eb1b0_0    defaults
session-info              1.0.0                    pypi_0    pypi
setuptools                59.8.0                   pypi_0    pypi
sip                       6.6.2           py310h6a678d5_0    defaults
six                       1.16.0             pyhd3eb1b0_1    defaults
sortedcontainers          2.4.0                    pypi_0    pypi
sqlite                    3.39.2               h5082296_0    defaults
stack_data                0.2.0              pyhd3eb1b0_0    defaults
statsmodels               0.13.2                   pypi_0    pypi
stdlib-list               0.8.0                    pypi_0    pypi
tblib                     1.7.0                    pypi_0    pypi
threadpoolctl             3.1.0                    pypi_0    pypi
tk                        8.6.12               h1ccaba5_0    defaults
toml                      0.10.2             pyhd3eb1b0_0    defaults
toolz                     0.12.0                   pypi_0    pypi
tornado                   6.1                      pypi_0    pypi
tqdm                      4.64.1                   pypi_0    pypi
traitlets                 5.1.1              pyhd3eb1b0_0    defaults
typing-extensions         4.3.0                    pypi_0    pypi
tzdata                    2022c                h04d1e81_0    defaults
umap-learn                0.5.3                    pypi_0    pypi
urllib3                   1.26.12                  pypi_0    pypi
wcwidth                   0.2.5              pyhd3eb1b0_0    defaults
wheel                     0.37.1             pyhd3eb1b0_0    defaults
xz                        5.2.5                h7f8727e_1    defaults
yarl                      1.8.1                    pypi_0    pypi
zeromq                    4.3.4                h2531618_0    defaults
zict                      2.2.0                    pypi_0    pypi
zlib                      1.2.12               h5eee18b_3    defaults
zstd                      1.5.2                ha4553b6_0    defaults
```

```shell
root:~$ conda list -n trajectory

# packages in environment at /home/shpc_100839/miniconda3/envs/trajectory:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                  2_kmp_llvm    conda-forge
alabaster                 0.7.12             pyhd3eb1b0_0    defaults
anndata                   0.8.0              pyhd8ed1ab_1    conda-forge
asttokens                 2.0.5              pyhd3eb1b0_0    defaults
attrs                     21.4.0             pyhd3eb1b0_0    defaults
babel                     2.9.1              pyhd3eb1b0_0    defaults
backcall                  0.2.0              pyhd3eb1b0_0    defaults
brotli                    1.0.9                h5eee18b_7    defaults
brotli-bin                1.0.9                h5eee18b_7    defaults
brotlipy                  0.7.0           py39h27cfd23_1003    defaults
c-ares                    1.18.1               h7f8727e_0    defaults
ca-certificates           2022.07.19           h06a4308_0    defaults
cached-property           1.5.2                      py_0    defaults
cellrank                  1.5.1              pyhdfd78af_0    bioconda
cellrank-krylov           1.5.1              pyh1e54042_0    bioconda
certifi                   2022.9.14        py39h06a4308_0    defaults
cffi                      1.15.1           py39h74dc2b5_0    defaults
charset-normalizer        2.0.4              pyhd3eb1b0_0    defaults
click                     8.0.4            py39h06a4308_0    defaults
colorama                  0.4.5            py39h06a4308_0    defaults
cryptography              37.0.1           py39h9ce1e76_0    defaults
cycler                    0.11.0             pyhd3eb1b0_0    defaults
dbus                      1.13.18              hb2f20db_0    defaults
debugpy                   1.5.1            py39h295c915_0    defaults
decorator                 5.1.1              pyhd3eb1b0_0    defaults
docrep                    0.3.2              pyh44b312d_0    conda-forge
docutils                  0.18.1           py39h06a4308_3    defaults
entrypoints               0.4              py39h06a4308_0    defaults
executing                 0.8.3              pyhd3eb1b0_0    defaults
expat                     2.4.9                h6a678d5_0    defaults
fftw                      3.3.10          mpi_openmpi_hdeb57f9_5    conda-forge
fontconfig                2.13.1               h6c09931_0    defaults
fonttools                 4.25.0             pyhd3eb1b0_0    defaults
freetype                  2.11.0               h70c0345_0    defaults
future                    0.18.2           py39h06a4308_1    defaults
giflib                    5.2.1                h7b6447c_0    defaults
glib                      2.69.1               h4ff587b_1    defaults
gmp                       6.2.1                h295c915_3    defaults
gst-plugins-base          1.14.0               h8213a91_2    defaults
gstreamer                 1.14.0               h28cd5cc_2    defaults
h5py                      3.7.0           nompi_py39hd51670d_101    conda-forge
hdf5                      1.12.2          mpi_openmpi_hb3f3608_0    conda-forge
hypre                     2.25.0          mpi_openmpi_ha709252_0    conda-forge
icu                       58.2                 he6710b0_3    defaults
idna                      3.3                pyhd3eb1b0_0    defaults
imagesize                 1.4.1            py39h06a4308_0    defaults
importlib-metadata        4.11.3           py39h06a4308_0    defaults
importlib_metadata        4.11.3               hd3eb1b0_0    defaults
iniconfig                 1.1.1              pyhd3eb1b0_0    defaults
ipykernel                 6.15.2           py39h06a4308_0    defaults
ipython                   8.4.0            py39h06a4308_0    defaults
jedi                      0.18.1           py39h06a4308_1    defaults
jinja2                    3.0.3              pyhd3eb1b0_0    defaults
joblib                    1.1.0              pyhd3eb1b0_0    defaults
jpeg                      9e                   h7f8727e_0    defaults
jupyter_client            7.3.5            py39h06a4308_0    defaults
jupyter_core              4.10.0           py39h06a4308_0    defaults
kiwisolver                1.4.2            py39h295c915_0    defaults
krb5                      1.19.2               hac12032_0    defaults
lcms2                     2.12                 h3be6417_0    defaults
ld_impl_linux-64          2.38                 h1181459_1    defaults
lerc                      3.0                  h295c915_0    defaults
libblas                   3.9.0           16_linux64_openblas    conda-forge
libbrotlicommon           1.0.9                h5eee18b_7    defaults
libbrotlidec              1.0.9                h5eee18b_7    defaults
libbrotlienc              1.0.9                h5eee18b_7    defaults
libcblas                  3.9.0           16_linux64_openblas    conda-forge
libclang                  10.0.1          default_hb85057a_2    defaults
libcurl                   7.84.0               h91b91d3_0    defaults
libdeflate                1.8                  h7f8727e_5    defaults
libedit                   3.1.20210910         h7f8727e_0    defaults
libev                     4.33                 h7f8727e_1    defaults
libevent                  2.1.12               h8f2d780_0    defaults
libffi                    3.3                  he6710b0_2    defaults
libgcc-ng                 12.1.0              h8d9b700_16    conda-forge
libgfortran-ng            11.2.0               h00389a5_1    defaults
libgfortran5              11.2.0               h1234567_1    defaults
liblapack                 3.9.0           16_linux64_openblas    conda-forge
libllvm10                 10.0.1               hbcb73fb_5    defaults
libllvm11                 11.1.0               h9e868ea_5    defaults
libnghttp2                1.46.0               hce63b2e_0    defaults
libopenblas               0.3.21          pthreads_h78a6416_3    conda-forge
libpng                    1.6.37               hbc83047_0    defaults
libpq                     12.9                 h16c4e8d_3    defaults
libsodium                 1.0.18               h7b6447c_0    defaults
libssh2                   1.10.0               h8f2d780_0    defaults
libstdcxx-ng              12.1.0              ha89aaad_16    conda-forge
libtiff                   4.4.0                hecacb30_0    defaults
libuuid                   1.0.3                h7f8727e_2    defaults
libwebp                   1.2.2                h55f646e_0    defaults
libwebp-base              1.2.2                h7f8727e_0    defaults
libxcb                    1.15                 h7f8727e_0    defaults
libxkbcommon              1.0.1                hfa300c1_0    defaults
libxml2                   2.9.14               h74e7548_0    defaults
libxslt                   1.1.35               h4e12654_0    defaults
libzlib                   1.2.12               h166bdaf_3    conda-forge
llvm-openmp               14.0.6               h9e868ea_0    defaults
llvmlite                  0.38.1           py39h7d9a04d_0    conda-forge
loompy                    3.0.6                      py_0    conda-forge
lz4-c                     1.9.3                h295c915_1    defaults
markupsafe                2.1.1            py39h7f8727e_0    defaults
matplotlib                3.5.2            py39h06a4308_0    defaults
matplotlib-base           3.5.2            py39hf590b9c_0    defaults
matplotlib-inline         0.1.6            py39h06a4308_0    defaults
metis                     5.1.0                hf484d3e_4    defaults
mpfr                      4.1.0                h9202a9a_1    conda-forge
mpi                       1.0                     openmpi    defaults
mpi4py                    3.1.3            py39h5418507_2    conda-forge
mumps-include             5.2.1               ha770c72_11    conda-forge
mumps-mpi                 5.2.1               hfb3545b_11    conda-forge
munkres                   1.1.4                      py_0    defaults
natsort                   7.1.1              pyhd3eb1b0_0    defaults
ncurses                   6.3                  h5eee18b_3    defaults
nest-asyncio              1.5.5            py39h06a4308_0    defaults
networkx                  2.8.4            py39h06a4308_0    defaults
nspr                      4.33                 h295c915_0    defaults
nss                       3.74                 h0370c37_0    defaults
numba                     0.55.2           py39h66db6d7_0    conda-forge
numpy                     1.22.4           py39hc58783e_0    conda-forge
numpy_groupies            0.9.19             pyhd8ed1ab_0    conda-forge
openmpi                   4.1.4              ha1ae619_101    conda-forge
openssl                   1.1.1q               h7f8727e_0    defaults
packaging                 21.3               pyhd3eb1b0_0    defaults
pandas                    1.5.0            py39h4661b88_0    conda-forge
parmetis                  4.0.3             he9a3056_1005    conda-forge
parso                     0.8.3              pyhd3eb1b0_0    defaults
patsy                     0.5.2            py39h06a4308_1    defaults
pcre                      8.45                 h295c915_0    defaults
petsc                     3.17.4          real_h4502189_101    conda-forge
petsc4py                  3.17.4          real_he586954_100    conda-forge
pexpect                   4.8.0              pyhd3eb1b0_3    defaults
pickleshare               0.7.5           pyhd3eb1b0_1003    defaults
pillow                    9.2.0            py39hace64e9_1    defaults
pip                       22.1.2           py39h06a4308_0    defaults
pluggy                    1.0.0            py39h06a4308_1    defaults
ply                       3.11             py39h06a4308_0    defaults
progressbar2              3.37.1           py39h06a4308_0    defaults
prompt-toolkit            3.0.20             pyhd3eb1b0_0    defaults
psutil                    5.9.0            py39h5eee18b_0    defaults
ptscotch                  6.0.9                h0a9c416_2    conda-forge
ptyprocess                0.7.0              pyhd3eb1b0_2    defaults
pure_eval                 0.2.2              pyhd3eb1b0_0    defaults
py                        1.11.0             pyhd3eb1b0_0    defaults
pycparser                 2.21               pyhd3eb1b0_0    defaults
pygam                     0.8.0                      py_0    conda-forge
pygments                  2.11.2             pyhd3eb1b0_0    defaults
pygpcca                   1.0.3            py39hf3d152e_1    conda-forge
pynndescent               0.5.4              pyhd3eb1b0_0    defaults
pyopenssl                 22.0.0             pyhd3eb1b0_0    defaults
pyparsing                 3.0.9            py39h06a4308_0    defaults
pyqt                      5.15.7           py39h6a678d5_1    defaults
pyqt5-sip                 12.11.0          py39h6a678d5_1    defaults
pysocks                   1.7.1            py39h06a4308_0    defaults
pytest                    7.1.2            py39h06a4308_0    defaults
pytest-runner             5.3.1              pyhd3eb1b0_0    defaults
python                    3.9.13               haa1d7c7_1    defaults
python-dateutil           2.8.2              pyhd3eb1b0_0    defaults
python-utils              3.3.3            py39h06a4308_0    defaults
python_abi                3.9                      2_cp39    conda-forge
pytz                      2022.1           py39h06a4308_0    defaults
pyzmq                     23.2.0           py39h6a678d5_0    defaults
qt-main                   5.15.2               h327a75a_7    defaults
qt-webengine              5.15.9               hd2b0992_4    defaults
qtwebkit                  5.212                h4eab89a_4    defaults
readline                  8.1.2                h7f8727e_1    defaults
requests                  2.28.1           py39h06a4308_0    defaults
scalapack                 2.2.0                h67de57e_1    conda-forge
scanpy                    1.9.1              pyhd8ed1ab_0    conda-forge
scikit-learn              1.1.1            py39h6a678d5_0    defaults
scipy                     1.9.1            py39h8ba3f38_0    conda-forge
scotch                    6.0.9                hb2e6521_2    conda-forge
scvelo                    0.2.4              pyhdfd78af_0    bioconda
seaborn                   0.11.2             pyhd3eb1b0_0    defaults
session-info              1.0.0              pyhd8ed1ab_0    conda-forge
setuptools                63.4.1           py39h06a4308_0    defaults
sip                       6.6.2            py39h6a678d5_0    defaults
six                       1.16.0             pyhd3eb1b0_1    defaults
slepc                     3.17.2          real_h624bf36_100    conda-forge
slepc4py                  3.17.2          real_hb9e9d30_100    conda-forge
snowballstemmer           2.2.0              pyhd3eb1b0_0    defaults
sphinx                    5.0.2            py39h06a4308_0    defaults
sphinxcontrib-applehelp   1.0.2              pyhd3eb1b0_0    defaults
sphinxcontrib-devhelp     1.0.2              pyhd3eb1b0_0    defaults
sphinxcontrib-htmlhelp    2.0.0              pyhd3eb1b0_0    defaults
sphinxcontrib-jsmath      1.0.1              pyhd3eb1b0_0    defaults
sphinxcontrib-qthelp      1.0.3              pyhd3eb1b0_0    defaults
sphinxcontrib-serializinghtml 1.1.5              pyhd3eb1b0_0    defaults
sqlite                    3.39.2               h5082296_0    defaults
stack_data                0.2.0              pyhd3eb1b0_0    defaults
statsmodels               0.13.2           py39h7f8727e_0    defaults
stdlib-list               0.7.0                      py_2    conda-forge
suitesparse               5.10.1               h9e50725_1    conda-forge
superlu                   5.2.2                h00795ac_0    conda-forge
superlu_dist              7.2.0                h34f6f4d_0    conda-forge
tbb                       2021.5.0             hd09550d_0    defaults
threadpoolctl             2.2.0              pyh0d69192_0    defaults
tk                        8.6.12               h1ccaba5_0    defaults
toml                      0.10.2             pyhd3eb1b0_0    defaults
tomli                     2.0.1            py39h06a4308_0    defaults
tornado                   6.2              py39h5eee18b_0    defaults
tqdm                      4.64.0           py39h06a4308_0    defaults
traitlets                 5.1.1              pyhd3eb1b0_0    defaults
typing_extensions         4.3.0            py39h06a4308_0    defaults
tzdata                    2022c                h04d1e81_0    defaults
umap-learn                0.5.3            py39hf3d152e_0    conda-forge
urllib3                   1.26.11          py39h06a4308_0    defaults
wcwidth                   0.2.5              pyhd3eb1b0_0    defaults
wheel                     0.37.1             pyhd3eb1b0_0    defaults
wrapt                     1.14.1           py39h5eee18b_0    defaults
xz                        5.2.6                h5eee18b_0    defaults
yaml                      0.2.5                h7b6447c_0    defaults
zeromq                    4.3.4                h2531618_0    defaults
zipp                      3.8.0            py39h06a4308_0    defaults
zlib                      1.2.12               h5eee18b_3    defaults
zstd                      1.5.2                ha4553b6_0    defaults
```
