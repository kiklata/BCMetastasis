#!/bin/bash
# STAR env
set -e # Exit immediately if a command exits with a non-zero status

refer=/home/zhepan/Reference/STAR_index
whitelist=/home/zhepan/Reference/STAR_index/3M-february-2018.txt

for sample in $(ls /home/zhepan/BCLM/count | grep \\-1);
do
  bam=/home/zhepan/BCLM/count/$sample/outs/*.bam
  STAR --runThreadN 8 \
       --genomeDir $refer \
       --soloType CB_UMI_Simple \
       --readFilesType SAM SE --readFilesIn $bam  \
       --readFilesCommand samtools view -F 0x100 \
       --soloInputSAMattrBarcodeSeq CR UR \
       --soloInputSAMattrBarcodeQual CY UY \
       --soloCBwhitelist $whitelist \
       --soloCellFilter EmptyDrops_CR \
       --soloFeatures SJ Velocyto\
       --outFileNamePrefix /home/zhepan/BCLM/splicing/$sample/ \
       --outSAMtype None
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done
