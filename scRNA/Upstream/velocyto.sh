#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

gtf=/home/zhepan/Reference/Cellranger/genes.gtf
rmsk=/home/zhepan/Reference/GRCh38/mm10_rmsk.gtf

for sample in $(ls /home/zhepan/BCLM/count| grep \\-1);
do
  velocyto run10x -m $rmsk \
    -@ 8 \
    -o /home/zhepan/BCLM/velocyto/$sample \
    /home/zhepan/BCLM/count/$sample \
    $gtf 
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done

