## use cellbender env

sample=/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/BrM/BrM4

conda run -n cellbender cellbender remove-background --input $sample/raw_feature_bc_matrix.h5 \
                                                     --output $sample/cellbender.h5 \
                                                     --checkpoint-mins 