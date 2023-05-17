library(reticulate)

### 调用子环境的python
use_condaenv(condaenv = "cnmf_env", required = T,conda = "/home/hsy/miniconda3/bin/conda")
py_config() #如果显示cnmf_env环境里面的python就OK

### 第一步 ################
#表达矩阵文件命名尽量和我的示例文件一致(HNSCC17.count.txt)，保证第一个点(.)前面是样本名称/ID。这些文件都保存在某个文件夹里面，比如示例中的"count_data"
source("1.R")
step1(dir_input = "count_data",dir_output = "res1",k=3:5,iteration = 50) #这里为了演示方便，取值都比较小

# ###
# 到这儿先停一停，step1会在dir_output下面为每一个样本生成一个文件夹。
# 此外还有一个k_selection文件夹以及一个k_selection.txt文件。
# k_selection文件夹里面有每一个样本选k值的图片（点击图片可以在浏览器打开），以此为依据，将k_selection.txt补充完整（可以直接在Rstudio打开）。
# k_selection.txt有两列，以逗号分割，第二列是空的，需要你填上去。
# ###


### 第二步 ################
source("2.R")
step2(dir_input = "res1",dir_output = "res2",dir_count = "count_data",usage_filter = 0.03,top_gene = 30,cor_min = 0,cor_max = 0.6)
#过滤阈值，10X 0.01 smart-seq 0.03，这是我用过的阈值，不一定适用于所有情况
#cor_min可以设为0，尽管真实范围会小于0；因为原文用的是0.6，所以我这里也限制最大相关系数为0.6，主要是为了图明显一些。试过10X数据，相关系数比较高(dropout导致的假象)，所以没有限制
#top_gene按照20 30 50来设计的
#有一个color参数，可以给样本涂色，格式是命名后的字符串向量；函数里面有一个默认的配色，也还OK。

### 实际分析中的一点点经验 ################
#1. 相关性聚类这一步需要多做几次，所以原始画图的数据我都导出了，方便你自己改图
#根据每个program的大致功能，如果发现另一种功能的program聚到某一种meta模块里面，这时可以将乱入的program删掉，再做一次相关性聚类。
#如果我们认定应该属于同一个meta模块的program分散在两个地方，可以试试调整聚类方法(参数clustering_method)，或者像"_1""_2"这样定义
#2. 对top基因做注释，尽量选择和研究背景有关联的基因集，比如研究肿瘤细胞异质性时，选择和肿瘤相关的基因集