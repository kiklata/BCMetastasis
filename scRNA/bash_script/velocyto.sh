## use velocyto env
sample=P1015S2
genes='/home/zhepan/Reference/Cellranger/genes/genes.gtf'
rmsk='/home/zhepan/Reference/GRCh38/hg38_rmsk.gtf'
cellrangerfold='/home/zhepan/Project/MultiOmics/data/snRNA/Result/'${sample}'/cellranger/'${sample}

velocyto run10x -m $rmsk -@ 8 -vvv $cellrangerfold $genes
