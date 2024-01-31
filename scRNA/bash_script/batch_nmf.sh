inhouse=/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/inHouse/
pub=/home/zhepan/Project/scRNA_BC_metastases/Data/SingleCell/public/

script_path=/home/zhepan/Project/scRNA_BC_metastases/Analysis/scRNA/py_script

#python $script_path/nmf.py --path $pub --site Primary --sample_list CID3921 CID3941 CID3948 CID3963 CID4066 CID4067 CID4290A CID4461 CID4463 CID4465 CID4471 CID4495 CID44971 CID44991 CID4513 CID4515 CID45171 CID4523 CID4530N CID4535

#python $script_path/nmf.py --path $pub --site BoM --sample_list BoM1 BoM2

#python $script_path/nmf.py --path $pub/BrM --site science --sample_list BrM1 BrM2 BrM3 BrM4
#python $script_path/nmf.py --path $pub/BrM --site cell --sample_list BrM1 BrM2 BrM3
#python $script_path/nmf.py --path $pub --site BrM --sample_list GSE143423
#python $script_path/nmf.py --path $pub --site BrM --sample_list GSE202501

#python $script_path/nmf.py --path $pub/LymphM --site EMBOJ --sample_list LymphM1 LymphM2 LymphM3 LymphM4 LymphM5 LymphM6 LymphM7
python $script_path/nmf.py --path $pub/LymphM --site NC --sample_list LymphM1 LymphM5 LymphM6 LymphM7 LymphM8
#python $script_path/nmf.py --path $pub/LymphM --site oncogenesis --sample_list LymphM1 LymphM2 LymphM3 LymphM4 LymphM5


#python $script_path/nmf.py --path $inhouse --site BoM --sample_list BoM1 BoM3
#python $script_path/nmf.py --path $inhouse --site BrM --sample_list BrM1 BrM2 BrM3
#python $script_path/nmf.py --path $inhouse --site LiverM --sample_list LiverM1 LiverM2 LiverM3 LiverM4 LiverM5 LiverM6
#python $script_path/nmf.py --path $inhouse --site LungM --sample_list LungM1