#!/bin/bash
# velocyto env
set -e # Exit immediately if a command exits with a non-zero status

gtf=/home/zhepan/Reference/Cellranger/genes/genes.gtf
rmsk=/home/zhepan/Reference/GRCh38/mm10_rmsk.gtf

for sample in $(ls /home/zhepan/BCLM/count| grep \\-1);
do
  velocyto run10x -v -m $rmsk \
    --samtools-threads 8 --samtools-memory 15827 \
    /home/zhepan/BCLM/count/$sample \
    $gtf 
  rm -rf /home/zhepan/BCLM/count/$sample/outs/cellsorted_*.bam
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done

