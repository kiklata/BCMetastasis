#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

project_p = /home/zhepan/BCLM

for sample in $(ls $project_p/fastq | grep \\-1);
do
  mixcr analyze 10x-5gex-full-length \
    -s hsa \
    -t 8 \
    $project_p/fastq/$sample/*.fastq.gz \
    $project_p/mixcr/$sample/report
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done
