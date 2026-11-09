# 不同物种的数据库不一样,此处为human 
#!/bin/bash

# ==========================================
# 项目：pySCENIC 基因通路调节网络分析 (TFA)
# ==========================================

# 1. 核心路径定义 (只需修改这部分) ----------------------
# 数据库基础路径
DB_DIR="${MG_PROJECT_ROOT}/scripts/SCENIC/database"
# 项目基础路径
PROJ_DIR="${MG_PROJECT_ROOT}/scripts/SCENIC"
# 结果保存路径
RESULT_DIR="$PROJ_DIR/results"

# 2. 数据库文件定义 -----------------------------------
TFS="$DB_DIR/allTFs_hg38.txt"
FEATHER=$(find $DB_DIR -name "hg38*.rankings.feather" | head -n 1) # 动态获取第一个匹配的feather文件
TBL="$DB_DIR/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

# 3. 输入文件定义 -------------------------------------
INPUT_LOOM="$PROJ_DIR/data/tmp_obj.loom"

# 4. 环境准备 -----------------------------------------
# 创建结果目录（如果不存在）
mkdir -p $RESULT_DIR

# 验证文件是否存在
ls $TFS $FEATHER $TBL $INPUT_LOOM || { echo "错误：核心数据库文件缺失"; exit 1; }

# 激活环境
source ${MG_CONDA_ACTIVATE}
conda activate pyscenic

# 5. pySCENIC 三步走流程 -------------------------------

# Step 1: GRN (Gene Regulatory Network) - 推断转录因子与靶基因的共表达模块
echo ">>> Starting Step 1: GRN Inference..."
pyscenic grn \
    --num_workers 20 \
    --output "$RESULT_DIR/adj.sample.tsv" \
    --method grnboost2 \
    "$INPUT_LOOM" \
    "$TFS"

# Step 2: CTX (Cistarget) - 利用基元富集分析精简模块，识别 Regulons
echo ">>> Starting Step 2: Cistarget (Pruning)..."
pyscenic ctx \
    "$RESULT_DIR/adj.sample.tsv" "$FEATHER" \
    --annotations_fname "$TBL" \
    --expression_mtx_fname "$INPUT_LOOM" \
    --mode "dask_multiprocessing" \
    --output "$RESULT_DIR/reg.csv" \
    --num_workers 20 \
    --mask_dropouts

# Step 3: AUCell - 在单细胞水平评估 Regulon 的活性评分
echo ">>> Starting Step 3: AUCell (Scoring)..."
pyscenic aucell \
    "$INPUT_LOOM" \
    "$RESULT_DIR/reg.csv" \
    --output "$RESULT_DIR/out_SCENIC.loom" \
    --num_workers 20

echo ">>> pySCENIC Pipeline Completed Successfully!"