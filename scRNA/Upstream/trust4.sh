#!/bin/bash
# ImmunoRepertoire env
set -e # Exit immediately if a command exits with a non-zero status

vdj=/home/zhepan/Reference/Trust4/hg38_bcrtcr.fa
ref=/home/zhepan/Reference/Trust4/human_IMGT_C.fa

for sample in $(ls /home/zhepan/BCLM/count | grep \\-2);
do
  bam=/home/zhepan/BCLM/count/$sample/outs/*.bam
  res=/home/zhepan/BCLM/trust4/$sample
  mkdir $res
  run-trust4 -f $vdj --ref $ref -b $bam --barcode CB \
             -t 8 \
             --od $res -o Trust
done