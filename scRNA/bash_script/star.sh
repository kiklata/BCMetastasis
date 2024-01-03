#!/bin/bash
# STAR env
set -e # Exit immediately if a command exits with a non-zero status

sample=P1015S2
datapath=/home/zhepan/Project/MultiOmics/data/snRNA/Result/$sample/cellranger/$sample

refer=/home/zhepan/Reference/STAR
whitelist=/home/zhepan/Reference/STAR/3M-february-2018.txt

bam=$datapath/outs/possorted_genome_bam.bam
STAR --runThreadN 8 \
     --genomeDir $refer \
     --soloType CB_UMI_Simple \
     --readFilesType SAM SE --readFilesIn $bam  \
     --readFilesCommand samtools view -F 0x100 \
     --soloInputSAMattrBarcodeSeq CR UR \
     --soloInputSAMattrBarcodeQual CY UY \
     --soloCBwhitelist $whitelist \
     --soloCellFilter EmptyDrops_CR \
     --soloFeatures Gene SJ Velocyto \
     --outFileNamePrefix $datapath/star \
     --outSAMtype None
