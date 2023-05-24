#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

Project=BCLM
workdir=/home/zhepan

ref=$workdir/Reference/Cellranger

vdj=$workdir/Reference/Trust4/hg38_bcrtcr.fa
ref_trust=$workdir/Reference/Trust4/human_IMGT_C.fa

gtf=$workdir/Reference/Cellranger/genes/genes.gtf
rmsk=$workdir/Reference/GRCh38/mm10_rmsk.gtf

refer=$workdir/Reference/STAR_index
whitelist=$workdir/Reference/STAR_index/3M-february-2018.txt

for sample in `ls $workdir/$Project/fastq`;
do  
  fastq=$workdir/$Project/fastq/$sample
  cellranger count --id $workdir/$Project/count/$sample \
                   --localcores=12 \
                   --localmem=124 \
                   --transcriptome=$ref \
                   --include-introns=false \
                   --nosecondary \
                   --fastqs=$fastq 
                   
  mixcr analyze 10x-5gex-full-length \
    -s hsa \
    -t 8 \
    ~/$Project/fastq/$sample/*.fastq.gz \
    ~/$Project/mixcr/$sample/report
    
  bam=$workdir/$Project/count/$sample/outs/*.bam
  res=$workdir/$Project/trust4/$sample
  
  run-trust4 -f $vdj --ref $ref_trust -b $bam --barcode CB \
             -t 8 \
             --od $res
             
  velocyto run10x -v -m $rmsk \
    --samtools-threads 8 --samtools-memory 65536 \
    $workdir/$Project/count/$sample \
    $gtf 
    
  rm -rf $workdir/$Project/count/$sample/outs/cellsorted_*.bam

  STAR --runThreadN 8 \
       --genomeDir $refer \
       --soloType CB_UMI_Simple \
       --readFilesType SAM SE --readFilesIn $bam  \
       --readFilesCommand samtools view -F 0x100 \
       --soloInputSAMattrBarcodeSeq CR UR \
       --soloInputSAMattrBarcodeQual CY UY \
       --soloCBwhitelist $whitelist \
       --soloCellFilter EmptyDrops_CR \
       --soloFeatures Gene SJ Velocyto\
       --outFileNamePrefix $workdir/$Project/splicing/$sample \
       --outSAMtype None
done

curl https://sctapi.ftqq.com/SCT167047TINLpe9UwpLd9LmqT04NylUuh.send?title=Done
