## in velocyto env

set -e # Exit immediately if a command exits with a non-zero status

samples=P1015S2

datapath=/home/zhepan/Project/MultiOmics/data/snRNA/$samples/Rawdata/

ref=/home/zhepan/Reference/Cellranger

cd /home/zhepan/Project/MultiOmics/data/snRNA/Result/$samples/cellranger

fastq=$datapath/$samples

cellranger count --id=${samples} \
                 --sample='P1015S2-1','P1015S2-2','P1015S2-3','P1015S2-4','P1015S2-5','P1015S2-6' \
                 --localcores=16 \
                 --localmem=180 \
                 --transcriptome=$ref \
                 --nosecondary --fastqs=$fastq

