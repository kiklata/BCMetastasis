#!/bin/bash
# CSP env

set -e # Exit immediately if a command exits with a non-zero status

project=/home/zhepan/BCLM
REGION_VCF=/home/zhepan/Reference/SNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz


for sample in $(ls /home/zhepan/BCLM/count | grep \\-1);

do
  path=count/$sample/outs
  
  bam=$project/$path/*.bam
  barcode=$project/$path/filtered_feature_bc_matrix/barcodes.tsv.gz
  out_dir=$project/snp/$sample
  mkdir $out_dir
  
  cellsnp-lite -s $bam -b $barcode -O $out_dir -R $REGION_VCF -p 8 --minMAF 0.1 --minCOUNT 20 --gzip
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=CSP
