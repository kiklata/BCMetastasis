#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

project_p = /home/zhepan/BCLM

gtf=/home/zhepan/Reference/Cellranger/genes.gtf
rmsk=/home/zhepan/Reference/GRCh38/mm10_rmsk.gtf

for sample in $(ls $project_p/count| grep \\-1);
do
  velocyto run10x -m $rmsk \
    -@ 8 \
    -o $project_p/velocyto/$sample \
    $project_p/count/$sample \
    $gtf 
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done

