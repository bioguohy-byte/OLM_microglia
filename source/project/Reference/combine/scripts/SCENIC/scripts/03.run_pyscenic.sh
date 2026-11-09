# 不同物种的数据库不一样,此处为human 
#!/bin/bash
dir=${MG_PROJECT_ROOT}/Reference/combine/scripts/SCENIC/database #设置自己的工作目录
tfs=$dir/allTFs_hg38.txt
feather=`find  ${MG_PROJECT_ROOT}/Reference/combine/scripts/SCENIC/database  -name "hg38*.rankings.feather"`
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
# make sure these database files are complete and error-free

input_loom=${MG_PROJECT_ROOT}/Reference/combine/scripts/SCENIC/data/tmp_obj.loom
ls $tfs  $feather  $tbl  

#conda activate 激活conda环境
source ${MG_CONDA_ACTIVATE}
conda activate pyscenic

# pyscenic 的3个步骤之 grn
pyscenic grn \
--num_workers 20 \
--output ${MG_PROJECT_ROOT}/Reference/combine/scripts/SCENIC/results/adj.sample.tsv \
--method grnboost2 \
$input_loom  \
$tfs 

#pyscenic 的3个步骤之 cistarget
pyscenic ctx \
${MG_PROJECT_ROOT}/Reference/combine/scripts/SCENIC/results/adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output ${MG_PROJECT_ROOT}/Reference/combine/scripts/SCENIC/results/reg.csv \
--num_workers 20  \
--mask_dropouts

#pyscenic 的3个步骤之 AUCell
pyscenic aucell \
$input_loom \
${MG_PROJECT_ROOT}/Reference/combine/scripts/SCENIC/results/reg.csv \
--output ${MG_PROJECT_ROOT}/Reference/combine/scripts/SCENIC/results/out_SCENIC.loom \
--num_workers 20 
