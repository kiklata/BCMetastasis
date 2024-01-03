#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

cd /home/zhepan/BCLM/count

ref=/home/zhepan/Reference/Cellranger

for sample in `ls /home/zhepan/BCLM/fastq`;
do  
  fastq=/home/zhepan/BCLM/fastq/$sample
  cellranger count --id $sample \
                   --localcores=12 \
                   --localmem=124 \
                   --transcriptome=$ref \
                   --include-introns=false \
                   --nosecondary \
                   --fastqs=$fastq 
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done

