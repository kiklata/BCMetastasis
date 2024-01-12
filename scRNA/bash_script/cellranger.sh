## in velocyto env

set -e # Exit immediately if a command exits with a non-zero status

for i in $(cat ~/bofile);
do
  samples=$i
  datapath=/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/public/raw_sc_BC/meta_Bone/$samples

  ref=/home/zhepan/Reference/Cellranger

  cd /home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/public/raw_sc_BC/meta_Bone/$samples

  fastq=$datapath

  cellranger count --id=${samples} \
                  --localcores=16 \
                  --localmem=180 \
                  --transcriptome=$ref \
                  --nosecondary --fastqs=$fastq
done