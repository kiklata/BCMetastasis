## use cellbender env

datapath='/home/zhepan/Project/MultiOmics/data/snRNA/P1013S2/Result/2.Summary/P1013S2'

## run in google colab to obtain checkpoint file, first upload raw h5 file to google drive
#conda run -n cellbender cellbender remove-background --input drive/MyDrive/raw_feature_bc_matrix.h5 \
#                                                     --output drive/MyDrive/cellbender_feature_bc_matrix.h5 \
#                                                     --cpu-threads 16 --cuda

## run in local due to limited memory (13GB) of colab, see https://github.com/broadinstitute/CellBender/issues/266

# local hash df7718350c remote hash 74c288e45e

#tar -zxvf ckpt.tar.gz
#rename 's/^74c288e45e/df7718350c/' 74c288e45e*
#tar -czvf newckpt.tar.gz df7718350c*
#rm -rf df7718350c*

## the code is corrupted for now, a GPU is needed

#conda run -n cellbender cellbender remove-background --input $datapath/raw_feature_bc_matrix.h5 \
#                                                     --output $datapath/cellbender_feature_bc_matrix.h5 \
#                                                     --cpu-threads 16 --checkpoint newckpt.tar.gz

## just download all file required in local PC, run in cellbender env, succeed                                               
conda run -n cellbender cellbender remove-background --input raw_feature_bc_matrix.h5 \
                                                     --output cellbender_feature_bc_matrix.h5 \
                                                     --cuda --estimator_multiple_cpu	

## P1013S2 UMAP display a poor perform, might not use it