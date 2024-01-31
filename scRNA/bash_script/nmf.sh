for i in $(cat /home/zhepan/Project/scRNA_BC_metastases/Result/nmf/tumor/nmf_count.txt);
do
  rm -rf $i/cNMF;
  python /home/zhepan/nmf_count.py --path $i --count count.txt --name cNMF;
done