## use cellbender env

## google colab failed due to limited memory (13GB) of colab, run it in local. see https://github.com/broadinstitute/CellBender/issues/266

## just download all file required in local PC, run in cellbender env, succeed                                               
conda run -n cellbender cellbender remove-background --input raw_feature_bc_matrix.h5 \
                                                     --output cellbender_feature_bc_matrix.h5 \
                                                     --cuda

## P1013S2 UMAP display a poor perform, might not use it