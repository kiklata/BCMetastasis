#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

refer=/home/zhepan/Reference/Cellranger
whitelist=/home/zhepan/Reference/STAR_index/3M-february-2018

for sample in $(ls /home/zhepan/BCLM/fastq | grep \\-1);
do
  bam=/home/zhepan/BCLM/count/$sample/outs/*.bam
  STAR --runThreadN 8 \
       --genomeDir $refer \
       --soloType CB_UMI_Simple \
       --readFilesIn $bam \
       --readFilesCommand samtools view -F 0x100 \
       --readFilesType SAM PE \
       --soloInputSAMattrBarcodeSeq CR UR \
       --soloInputSAMattrBarcodeQual CY UY \
       --soloCBwhitelist $whitelist \
       --soloFeatures Gene SJ
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done
