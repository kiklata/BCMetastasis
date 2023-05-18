#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

project_p = /home/zhepan/BCLM
refer=/home/zhepan/Reference/STAR_index
whitelist=/home/zhepan/Reference/STAR_index/3M-february-2018.txt

for sample in $(ls $project_p/fastq | grep \\-1);
do
  bam=$project_p/count/$sample/outs/*.bam
  sam=$project_p/count/$sample/outs/tmp.sam
  samtools view -F0x100 $bam > $sam
  STAR --runThreadN 8 \
       --genomeDir $refer \
       --soloType CB_UMI_Simple \
       --readFilesIn $sam \
       --readFilesType SAM PE \
       --soloInputSAMattrBarcodeSeq CR UR \
       --soloInputSAMattrBarcodeQual CY UY \
       --soloCBwhitelist $whitelist \
       --soloFeatures Gene SJ \
       --outFileNamePrefix $project_p/splicing/$sample
  rm -rf $sam
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done
