# 先在terminal中激活环境
# source ${MG_CONDA_ACTIVATE}
# conda activate pyscenic

import os, sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv(__import__("os").path.join(__import__("os").environ["MG_PROJECT_ROOT"], "Reference/combine/scripts/SCENIC/data/tmp_obj.csv"));
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create(__import__("os").path.join(__import__("os").environ["MG_PROJECT_ROOT"], "Reference/combine/scripts/SCENIC/data/tmp_obj.loom"),x.X.transpose(),row_attrs,col_attrs);

