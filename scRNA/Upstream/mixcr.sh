#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

for sample in $(ls /home/zhepan/BCLM/fastq | grep \\-1);
do
  mixcr analyze 10x-5gex-full-length \
    -s hsa \
    -t 8 \
    ~/BCLM/fastq/$sample/*.fastq.gz \
    ~/BCLM/mixcr/$sample/report
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done
