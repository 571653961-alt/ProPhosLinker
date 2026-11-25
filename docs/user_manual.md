# ProPhosLinker 用户使用手册

> **版本**: v2.0 | **最后更新**: 2026-11-25

---

## 目录

1. [概述](#1-概述)
2. [安装](#2-安装)
   - [2.1 Docker 安装（推荐）](#21-docker-安装推荐)
   - [2.2 手动安装](#22-手动安装)
3. [输入数据格式](#3-输入数据格式)
   - [3.1 蛋白质组学数据](#31-蛋白质组学数据)
   - [3.2 磷酸化蛋白质组学数据](#32-磷酸化蛋白质组学数据)
   - [3.3 样本分组信息](#33-样本分组信息)
   - [3.4 蛋白质-磷酸化位点映射](#34-蛋白质-磷酸化位点映射)
4. [快速入门](#4-快速入门)
5. [命令行参数](#5-命令行参数)
6. [配置文件 (config.yaml)](#6-配置文件-configyaml)
   - [6.1 数据预处理](#61-数据预处理)
   - [6.2 模式分析](#62-模式分析)
   - [6.3 差异分析](#63-差异分析)
   - [6.4 功能分析](#64-功能分析)
7. [分析流程详解](#7-分析流程详解)
   - [步骤 0: 数据预处理](#步骤-0-数据预处理)
   - [步骤 1: 模式分析](#步骤-1-模式分析)
   - [步骤 2: 差异分析](#步骤-2-差异分析)
   - [步骤 3: 功能分析](#步骤-3-功能分析)
8. [输出文件说明](#8-输出文件说明)
9. [Case Study: 胰腺导管腺癌 (PDAC)](#9-case-study-胰腺导管腺癌-pdac)
10. [常见问题](#10-常见问题)
11. [引用](#11-引用)

---

## 1. 概述

ProPhosLinker 是一个整合蛋白质组学和磷酸化蛋白质组学的生物信息学分析工具。通过构建蛋白质-磷酸化相互作用网络，结合多组学数据，实现从数据预处理到功能解读的全流程自动化分析。

**核心功能：**

| 模块 | 功能 |
|------|------|
| 数据预处理 | 缺失值过滤、KNN插补、中位数归一化 |
| 模式分析 | Procrustes分析、NMF分子分型、WGCNA共表达网络 |
| 差异分析 | limma差异表达、磷酸化率分型、Bootstrap稳定差异网络 |
| 功能分析 | GO/KEGG富集、知识图谱功能网络、Hub蛋白网络 |

**适用场景：** 肿瘤 vs 正常组织的蛋白质组+磷酸化蛋白质组配对分析。

---

## 2. 安装

### 2.1 Docker 安装（推荐）

Docker 容器内置 Python 3.10 + R 4.1 + Neo4j 5.26，开箱即用。

```bash
# 拉取镜像
docker pull ppy1222/prophoslinker:v2

# 运行分析
docker run -it --rm \
  -v /path/to/your/data:/data \
  -v /path/to/output:/results \
  ppy1222/prophoslinker:v2 \
  --pro_file /data/protein_abundance.tsv \
  --phos_file /data/phosphoprotein_abundance.tsv \
  --sample_group /data/compare_groups.tsv \
  --mapping_file /data/protein_phosphoproSite.tsv \
  --group_comparing "T:N"
```

或使用 docker-compose：

```yaml
# docker-compose.yml
services:
  prophoslinker:
    image: ppy1222/prophoslinker:v2
    volumes:
      - ./casedata:/data
      - ./results:/results
    command: >
      --pro_file /data/protein_abundance.tsv
      --phos_file /data/phosphoprotein_abundance.tsv
      --sample_group /data/compare_groups.tsv
      --mapping_file /data/protein_phosphoproSite.tsv
      --group_comparing "T:N"
```

### 2.2 手动安装

**系统要求：**
- Python 3.10+
- R 4.1.2+
- Neo4j Community Edition 5.26+
- Java 17+

**安装步骤：**

```bash
# 1. 克隆仓库
git clone https://github.com/mzlab.research/ProPhosLinker.git
cd ProPhosLinker

# 2. 安装 Python 依赖
pip install -r requirements.txt

# 3. 安装 R 包
Rscript install_r_packages.R

# 4. 安装包本身
pip install -e .

# 5. 启动 Neo4j（需先安装）
neo4j start
```

---

## 3. 输入数据格式

所有输入文件均为 **TSV 格式**（制表符分隔）。

### 3.1 蛋白质组学数据

| 参数 | `--pro_file` |
|------|-------------|

| 列 | 内容 | 示例 |
|----|------|------|
| 第1列 | 基因名称 | `A1BG` |
| 第2列起 | 各样本表达量 | `254240956` |

> 缺失值标记为 `NA`。首行列名为 `Gene` 接各样本名。

**示例（前3行）：**

```
Gene	N_C3L_03123	N_C3N_03884	T_C3L_03635	T_C3N_00518
A1BG	254240956	207337793	430055999	530313237
A2M	494757222	468169804	551030445	890765355
```

### 3.2 磷酸化蛋白质组学数据

| 参数 | `--phos_file` |
|------|--------------|

| 列 | 内容 | 示例 |
|----|------|------|
| 第1列 | 磷酸化位点ID | `NP_000005.2_S928` |
| 第2列起 | 各样本表达量 | 数值 |

> 位点ID格式：`{RefSeq ID}_{磷酸化残基}{位点位置}`（如 `S`=丝氨酸, `T`=苏氨酸, `Y`=酪氨酸）。

**示例（前3行）：**

```
Index	N_C3L_03123	N_C3N_03884	T_C3L_03635
NP_000005.2_S928	1500000	1200000	2800000
NP_000011.2_S160	850000	920000	1100000
```

### 3.3 样本分组信息

| 参数 | `--sample_group` |
|------|-----------------|

**必须包含以下4列，列名不可修改：**

| 列名 | 说明 |
|------|------|
| `Protein_sample` | 蛋白质组学样本名 |
| `Phosphosite_sample` | 磷酸化组学样本名 |
| `Cor_sample` | 关联分析用样本名 |
| `group` | 分组标签（如 `T`=实验组, `N`=对照组） |

**示例：**

```
Protein_sample	Phosphosite_sample	Cor_sample	group
T_C3L_03394	T_C3L_03394	T_C3L_03394	T
N_C3L_03123	N_C3L_03123	N_C3L_03123	N
```

### 3.4 蛋白质-磷酸化位点映射

| 参数 | `--mapping_file` |
|------|-----------------|

**两列，无表头修改（注意第二列名为 `Index` 而非 `Phosphosite`）：**

| 列 | 内容 |
|----|------|
| 第1列 (`Gene`) | 基因名 |
| 第2列 (`Index`) | 磷酸化位点ID |

> 一个基因可对应多个磷酸化位点（一对多）。

**示例：**

```
Gene	Index
A2M	NP_000005.2_S928
ACVRL1	NP_000011.2_S160
ACVRL1	NP_000011.2_S161
```

---

## 4. 快速入门

以下使用项目自带的示例数据（`casedata/` 目录）演示完整流程。

### 最小化运行（默认参数）

```bash
ProPhosLinker \
  --pro_file casedata/protein_abundance.tsv \
  --phos_file casedata/phosphoprotein_abundance.tsv \
  --sample_group casedata/compare_groups.tsv \
  --mapping_file casedata/protein_phosphoproSite.tsv \
  --group_comparing "T:N"
```

输出将保存在 `./results/`。

### 调整差异阈值

```bash
ProPhosLinker \
  --pro_file casedata/protein_abundance.tsv \
  --phos_file casedata/phosphoprotein_abundance.tsv \
  --sample_group casedata/compare_groups.tsv \
  --mapping_file casedata/protein_phosphoproSite.tsv \
  --group_comparing "T:N" \
  --pro_FC 1.5 \
  --phos_FC 1.5 \
  --pro_diff_q_val 0.01 \
  --phos_diff_q_val 0.01
```

### 指定标识符类型

```bash
ProPhosLinker \
  ... \
  --identified_type UNIPROT       # 或 SYMBOL (默认)
```

### 只运行特定步骤

```bash
# 仅运行差异分析和功能分析（跳过预处理和模式分析）
ProPhosLinker ... --steps "2,3"
```

### 使用外部 Neo4j

```bash
ProPhosLinker \
  ... \
  --neo4j_path /opt/neo4j \
  --username neo4j \
  --password your_password
```

---

## 5. 命令行参数

### 必选参数

| 参数 | 说明 |
|------|------|
| `--pro_file` | 蛋白质组学数据文件路径 (TSV) |
| `--phos_file` | 磷酸化蛋白质组学数据文件路径 (TSV) |
| `--sample_group` | 样本分组文件路径 (TSV) |
| `--mapping_file` | 蛋白质-磷酸化位点映射文件路径 (TSV) |

### 差异分析参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--group_comparing` | `"T:N"` | 比较组，格式 `"实验组:对照组"` |
| `--identified_type` | `"SYMBOL"` | 基因标识符类型：`UNIPROT` 或 `SYMBOL` |
| `--pro_FC` | `2.0` | 蛋白质差异倍数阈值（线性，非 log2） |
| `--phos_FC` | `2.0` | 磷酸化差异倍数阈值（线性） |
| `--pro_diff_q_val` | `0.05` | 蛋白质差异显著性 q 值阈值 |
| `--phos_diff_q_val` | `0.05` | 磷酸化差异显著性 q 值阈值 |
| `--network_FC` | `2.0` | 网络分析差异倍数阈值 |
| `--disease` | `""` | 疾病名称（用于知识图谱检索） |

### 组学命名

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--omics1_name` | `"Protein"` | 组学1（蛋白质组）名称 |
| `--omics2_name` | `"PhosProtein"` | 组学2（磷酸化组）名称 |

### 系统参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--config` | 无 | YAML 配置文件路径（覆盖 CLI 参数） |
| `--steps` | `"all"` | 执行步骤：`0,1,2,3` 或 `all` |
| `--neo4j_path` | 无 | Neo4j 安装路径 |
| `--username` | `"neo4j"` | Neo4j 用户名 |
| `--password` | `"neo4j"` | Neo4j 密码 |
| `--version` | - | 显示版本号 |

---

## 6. 配置文件 (config.yaml)

通过 `--config /path/to/config.yaml` 加载。YAML 中的值覆盖命令行默认值。完整配置文件位于 `ProPhosLinker/config.yaml`。

### 6.1 数据预处理

```yaml
data_processing:
  pro_miss_value_ratio: 0.5        # 蛋白质缺失值过滤比例 (0-1)
  phos_miss_value_ratio: 0.5       # 磷酸化缺失值过滤比例 (0-1)
  pro_imputation_method: knn       # 插补方法: zero/min/mean/median/knn/gmm/MLE
  phos_imputation_method: knn
  knn_k: 5                         # KNN 近邻数
  pro_normalization_method: median # 归一化: median/quantile/PQN/none
  phos_normalization_method: median
```

### 6.2 模式分析

```yaml
pa:
  dim_rd_mtd: PCA                  # 降维方法

nmf:
  pro_filter_nums: 3000            # NMF 蛋白质过滤数
  phos_filter_nums: 3000           # NMF 磷酸化过滤数

wgcna:
  pro_filter_num: 5000             # WGCNA 蛋白质过滤数
  phos_filter_num: 5000            # WGCNA 磷酸化过滤数
  pro_SoftPower: auto              # 软阈值（auto=自动选择）
  phos_SoftPower: auto
  cor_method: spearman             # 相关性方法
  networkType: signed              # 网络类型
  minModuleSize: 30                # 最小模块大小
  mergingThresh: 0.25              # 模块合并阈值

module_cor_:
  threshold: 0.5                   # 模块-性状相关性阈值
  p_adj: 0.05                      # 显著性阈值
  method: spearman                 # 相关方法
```

### 6.3 差异分析

```yaml
phosrate_subtyping:
  robust_quantiles: [0.25, 0.50, 0.75]  # 分位数切点
  overlap_coefficients: 0.8             # 重叠系数阈值
  FDR: 0.05                             # 错误发现率
  spline_df: 3                          # 样条自由度
  grid_points: 100                      # 网格点数

differential_network:
  filter_num: 3000                      # 过滤特征数
  FC_threshold: 1.0                     # 差异倍数阈值 (log2)
  P_threshold: 0.05                     # 显著性阈值
  nBoots: 100                           # Bootstrap 次数
  stability: 0.8                        # 边稳定性阈值
  cor_method: spearman                  # 相关性方法
  R_threshold: 0.4                      # 相关系数阈值
```

### 6.4 功能分析

```yaml
functional_enrichment:
  GO_showCategory: 6                    # GO 显示类别数
  KEGG_showCategory: 18                 # KEGG 显示类别数
  pvalueCutoff: 0.05                    # 富集显著性阈值

kgbased_network:
  top_n: 6                              # 知识图谱网络 Top N
  visualization_num: 20                 # 可视化节点数
  max_phosphoSite: 5                    # 每蛋白最大磷酸化位点数
  node_filtering: betweenness           # 节点过滤方式
  layout: kk                            # 网络布局算法
  comm_detection: infomap               # 社区检测算法

hub_protein_network:
  top_n: 5                              # Hub 蛋白数量
  layout: kk                            # 布局算法
```

---

## 7. 分析流程详解

### 步骤 0: 数据预处理

**参数**: `--steps "0"` 或包含 `0`

**流程:**

```
原始数据 → 缺失值过滤 → 异常值处理 → KNN插补 → 归一化 → 保存
```

**关键操作：**

1. **缺失值过滤**：删除缺失比例超过 `miss_value_ratio`（默认 50%）的特征
2. **异常值处理**：识别并处理离群值
3. **KNN 插补**：基于 k 近邻（默认 k=5）填补剩余缺失值
4. **中位数归一化**：各样本缩放到全局中位数

**输出**: `0.Data_Preprocessing/`

### 步骤 1: 模式分析

**参数**: `--steps "1"` 或包含 `1`

**子步骤：**

**1.1 Procrustes 分析** — 评估蛋白质组与磷酸化蛋白质组的整体一致性。
- 输入：归一化后的两个组学矩阵
- 输出：Procrustes 变换图、统计量表

**1.2 NMF 分子分型** — 基于非负矩阵分解的样本聚类。
- 输入：Top N 变异最大的特征（默认 3000）
- 输出：一致性聚类图、样本分型表

**1.3 WGCNA** — 加权基因共表达网络分析。
- 输入：Top N 变异最大的特征（默认 5000）
- 输出：模块-性状热图、模块成员表、基因显著性

### 步骤 2: 差异分析

**参数**: `--steps "2"` 或包含 `2`

**子步骤：**

**2.1 组学差异分析** — limma 线性模型差异表达分析。
- 输入：归一化矩阵 + 分组信息
- 输出：火山图 + 全差异结果表 (`prodiff.tsv`, `phosdiff.tsv`)

**2.2 磷酸化率分型** — 以亲本蛋白质表达为条件，分析磷酸化位点响应模式。
- 算法流程：min-max overlap → 条件/全局路由 → 样条拟合 → limma → 一致性评估 → k-means 聚类
- 支持并行计算
- 输出：Sankey 图 + 热图 + 聚类报告

**2.3 稳定差异网络** — Bootstrap 重采样构建条件相关网络，比较差异边。
- 方法：Bootstrap (默认 100 次) + 相关性稳定性评估 + 条件网络 + 差异比较
- 输出：稳定性图 + 差异网络图 + 边表

### 步骤 3: 功能分析

**参数**: `--steps "3"` 或包含 `3`

**子步骤：**

**3.1 功能富集分析** — GO 和 KEGG 通路富集。
- 输入：差异表达基因列表
- 支持 UNIPROT / SYMBOL / ENSEMBL / ENTREZ 多种 ID 类型
- 输出：点图 + 富集结果表

**3.2 知识图谱功能网络** — 基于 Neo4j 的蛋白质相互作用网络。
- 步骤：Neo4j 查询 → 网络构建 → 社区检测 → 可视化
- 支持 8 种过滤模式、6 种社区检测算法、9 种中心性排序

**3.3 Hub 蛋白质功能网络** — 关键蛋白质的子网络分析。
- 复合中心性：degree (40%) + betweenness (30%) + closeness (20%) + eigenvector (10%)
- 每个 Hub 蛋白独立分析

---

## 8. 输出文件说明

### 输出目录结构

```
results/
├── 0.Data_Preprocessing/
│   ├── Protein_preprocessed.tsv          # 处理后蛋白质表达矩阵
│   ├── PhosProtein_preprocessed.tsv      # 处理后磷酸化表达矩阵
│   ├── sample_list.tsv                   # 样本列表
│   ├── mapping_list.tsv                  # 蛋白质-位点映射
│   ├── metadata.tsv                      # 样本元数据
│   ├── prodiff.tsv                       # 蛋白质差异结果
│   ├── phosdiff.tsv                      # 磷酸化差异结果
│   └── input_overview.txt                # 数据概览报告
│
├── 1.Pattern_Analysis/
│   ├── 1.1Procrustes/
│   │   ├── procrustes_plot.png/pdf       # Procrustes 图
│   │   └── procrustes_stats.tsv          # 统计量
│   ├── 1.2Molecular_Subtyping/
│   │   ├── nmf_consensus.png/pdf         # NMF 一致性聚类
│   │   └── nmf_membership.tsv            # 样本分型
│   └── 1.3WGCNA/
│       ├── module_trait_heatmap.png/pdf  # 模块-性状热图
│       ├── module_membership.tsv         # 模块成员
│       └── gene_significance.tsv         # 基因显著性
│
├── 2.Differential_Analysis/
│   ├── 2.1Omics_Differential_Analysis/
│   │   ├── volcano_plot.png/pdf          # 火山图
│   │   ├── quadrant_plot.png/pdf         # 四象限图
│   │   └── diff_results.tsv             # 全差异表
│   ├── 2.2PhosRate_Subtyping/
│   │   ├── sankey_plot.png/pdf           # Sankey 图
│   │   ├── heatmap.png/pdf               # 表达热图
│   │   └── cluster_report.tsv            # 聚类报告
│   └── 2.3Stable_Differential_Network/
│       ├── stability_plot.png/pdf        # 稳定性评估
│       ├── diff_network.png/pdf          # 差异网络图
│       └── network_edges.tsv             # 边表
│
└── 3.Functional_Analysis/
    ├── 3.1Omics_Function_Enrichment/
    │   ├── GO_dotplot.png/pdf            # GO 富集点图
    │   ├── KEGG_dotplot.png/pdf          # KEGG 富集点图
    │   └── enrichment_results.tsv        # 富集结果表
    ├── 3.2KGbased_Functional_Network/
    │   ├── Node_info.tsv                 # 节点信息 (13列)
    │   ├── Edge_info.tsv                 # 边信息
    │   ├── community_modules.tsv         # 社区模块
    │   └── network_visualization.png/pdf # 网络可视化
    └── 3.3Hub_Protein_Functional_Network/
        └── {hub_gene}/
            └── hub_network.png/pdf       # Hub 蛋白子网络
```

### 关键输出文件详解

**差异结果表 (`prodiff.tsv` / `phosdiff.tsv`)**

| 列 | 说明 |
|----|------|
| `logFC` | log2 差异倍数 |
| `AveExpr` | 平均表达量 |
| `t` | t 统计量 |
| `P.Value` | 原始 p 值 |
| `adj.P.Val` | 校正后 q 值 |
| `B` | B 统计量 |

**节点信息表 (`Node_info.tsv`)**

包含 13 列：基因名、差异状态（上调/下调/无差异）、差异倍数、显著性、磷酸化位点信息等。

---

## 9. Case Study: 胰腺导管腺癌 (PDAC)

本节以胰腺导管腺癌（Pancreatic Ductal Adenocarcinoma, PDAC）队列为例，展示 ProPhosLinker 从原始数据到功能解读的完整分析流程。该案例数据存放于项目 `casedata/` 目录下。

### 9.1 数据集概览

| 数据项 | 描述 |
|--------|------|
| **疾病** | 胰腺导管腺癌 (PDAC) |
| **样本量** | 215 例：140 肿瘤 (T) vs 75 正常 (N) |
| **蛋白质组** | 8,754 种蛋白质（11,646 行，含基因名） |
| **磷酸化组** | 16,651 个磷酸化位点，分布于 4,768 个磷酸化蛋白 |
| **蛋白质-磷酸化位点映射** | 51,399 条对应关系 |
| **临床信息** | 140 例患者表型数据（年龄、性别、肿瘤大小、分期等） |

**输入文件：**

| 文件 | 路径 | 大小 | 说明 |
|------|------|------|------|
| 蛋白丰度矩阵 | `casedata/protein_abundance.tsv` | 33.8 MB | 基因 × 样本，log2 强度值 |
| 磷酸化位点丰度矩阵 | `casedata/phosphoprotein_abundance.tsv` | 89.6 MB | 磷酸化位点 × 样本 |
| 样本分组 | `casedata/compare_groups.tsv` | 8.4 KB | 215 行，标注 T/N 分组 |
| 蛋白-位点映射 | `casedata/protein_phosphoproSite.tsv` | 1.3 MB | 51,399 对映射关系 |
| 临床信息 | `casedata/clinical_table_140.tsv` | 65 KB | 140 例患者的临床表型 |

### 9.2 运行命令

```bash
prophoslinker \
  --pro_file "casedata/protein_abundance.tsv" \
  --phos_file "casedata/phosphoprotein_abundance.tsv" \
  --sample_group "casedata/compare_groups.tsv" \
  --mapping_file "casedata/protein_phosphoproSite.tsv" \
  --outdir "results_new_V3" \
  --group_comparing "T:N" \
  --password "12345678" \
  --metadata_file "casedata/clinical_table_140.tsv" \
  --disease "pancreatic ductal adenocarcinoma"
```

**参数说明：**

| 参数 | 值 | 含义 |
|------|-----|------|
| `--group_comparing "T:N"` | Tumor vs Normal | 以 T（肿瘤）为实验组，N（正常）为对照组 |
| `--disease` | pancreatic ductal adenocarcinoma | 用于知识图谱构建的疾病上下文 |
| `--metadata_file` | clinical_table_140.tsv | 携带表型数据以支持 WGCNA 模块-性状关联分析 |

分析在单机完成，完整运行约需 2-4 小时（取决于硬件配置）。

### 9.3 各步骤产出

#### 步骤 0：数据预处理

**关键输出：**

| 文件 | 内容 |
|------|------|
| `input_overview.txt` | 完整的输入数据摘要报告 |
| `mapping_list.tsv` | 51,399 对蛋白质-磷酸化位点映射表 |
| `protein_preprocessed.tsv` | 清洗后的蛋白表达矩阵 |
| `phosprotein_preprocessed.tsv` | 清洗后的磷酸化位点表达矩阵 |
| `metadata.tsv` | 格式化的临床元数据 |
| `sample_list.tsv` | 样本分组清单（T=140, N=75） |

**数据质量：** 8,754 种蛋白与 4,768 个磷酸化蛋白中共有 3,998 个重叠蛋白，平均每个磷酸化蛋白含 3.49 个磷酸化位点。

#### 步骤 1：模式分析

**1.1 Procrustes 分析** — 评估蛋白质组与磷酸化组的数据一致性，生成 Procrustes 旋转图 (`procrustes.png`)。

**1.2 分子亚型** — 基于多组学数据的分子分型：
- `molecular_subtyping_sankey.png` — 桑基图展示蛋白亚型 ↔ 磷酸化亚型的对应关系
- `molecular_subtyping_sankey_by_group.png` — 按 T/N 分组展示亚型流转

**1.3 WGCNA 加权基因共表达网络**：
- 蛋白水平：`Protein_soft_threshold.png`（软阈值选择）、`Protein_Sample_Dendrogram_and_Metadata_Heatmap.png`（样本聚类 + 表型热图）
- 磷酸化水平：同上，对应 `PhosProtein_*` 文件
- 模块-性状关联：`module_metadata_correlation_pheatmap.png`（多组学模块与临床性状的相关性热图）

#### 步骤 2：差异分析

**2.1 组学差异分析** — 蛋白水平与磷酸化水平的差异表达四象限图 (`omics_diff_quadrant_plot.png`)：
- 右上象限：蛋白上调 + 磷酸化上调（协同激活）
- 左上/右下象限：蛋白与磷酸化方向不一致（翻译后调控主导）

**2.2 磷酸化速率亚型分析** — `ratio_clusters.png` 展示磷酸化速率聚类模式，`site_routing.png` 展示位点归类路径。16 个代表性位点在 `case_sites/` 下生成独立差异图（如 `YAP1_S164.png`、`CTNNB1_S675.png` 等），展示单个磷酸化位点的肿瘤 vs 正常差异模式。

**2.3 稳定差异网络**：
- 自举法稳定性检验：`Bootnet_T.png`、`Bootnet_N.png`（T 和 N 各自的网络稳定性）
- 整体网络 GO/KEGG 富集：突出 PDAC 相关通路
- 网络聚类：`network_by_top_modules.png` 展示整体模块结构
- 子网络分析：cluster 1/3/4/5/6/8/9 每个子网络均有独立可视化（包含差异网络图、配体-受体网络、转录因子调控网络）及功能富集分析

#### 步骤 3：功能分析

**3.1 组学功能富集**：
- `Protein_PhosProtein_GO_enrichment.png` — 整合蛋白+磷酸化差异基因的 GO 富集气泡图
- `Protein_PhosProtein_KEGG_enrichment.png` — KEGG 通路富集结果

**3.2 知识图谱功能网络** — 基于 Neo4j 图数据库，整合 PhosphoSitePlus 和 TRRUST 调控关系，以 PDAC 为疾病上下文构建的功能网络。

社区检测识别出多个功能模块（module 3-9），每个模块包含：
- `pro_phosphopro_network_ALL.png` — 模块全景网络图（含蛋白节点、磷酸化节点、调控边）
- `pro_phosphopro_network_LR.png` — 配体-受体相互作用子网络
- `pro_phosphopro_network_TF.png` — 转录因子调控子网络
- `Protein_PhosProtein_GO_enrichment.png` — 模块特异的 GO 功能注释
- `Protein_PhosProtein_KEGG_enrichment.png` — 模块特异的通路富集

### 9.4 结果解读要点

本次 PDAC 分析揭示了以下生物学发现方向：

1. **翻译后调控主导**：差异分析四象限图中，有大量蛋白不变但磷酸化水平显著变化的位点（左上/右下象限），提示 PDAC 中翻译后修饰层面（而非蛋白表达层面）的调控占据重要地位。

2. **关键磷酸化位点**：`case_sites/` 下的 YAP1、CTNNB1、COPA 等位点已知在胰腺癌中发挥重要作用，磷酸化速率亚型分析揭示了它们在不同样本中的调控模式差异。

3. **功能模块**：知识图谱网络聚类识别出多个功能模块，各模块富集于特定的 KEGG 通路和 GO 术语，可对应论文正文 2.6 节中的功能解读。

4. **临床关联**：WGCNA 模块-性状关联热图展示了特定共表达模块与临床特征（肿瘤大小、侵袭程度等）的显著性关联。

---

## 10. 常见问题

### Q: 运行报错 "Neo4j connection refused"

**A:** Neo4j 未启动或启动未完成。Docker 容器会自动等待最多 60 秒。手动安装时需先运行 `neo4j start`。

### Q: 内存不足 (out of memory)

**A:** 增大 Docker 内存限制：`docker run --memory 16g ...`。WGCNA 和差异网络步骤内存消耗较大，建议 ≥ 8GB。

### Q: R 包安装失败

**A:** 在中国大陆使用清华镜像（`install_r_packages.R` 已默认配置）。如仍有问题，尝试手动安装失败的具体包：

```r
install.packages("包名", repos = "https://cloud.r-project.org")
```

### Q: 如何只重新运行某个步骤？

**A:** 使用 `--steps` 参数。例如只运行步骤 2 和 3：

```bash
ProPhosLinker ... --steps "2,3"
```

已完成的步骤不会重复计算，除非删除对应输出目录。

### Q: 输入文件中的样本名是否需要完全匹配？

**A:** 是的。`sample_group` 中的样本名必须与 `pro_file` 和 `phos_file` 中的列名完全一致。程序会验证并报告不匹配的样本。

### Q: 支持哪些物种？

**A:** 默认使用人类数据库（`org.Hs.eg.db`）。如需其他物种，需修改配置文件或源码。

### Q: 如何自定义输出目录？

**A:** 通过 `--outdir` 参数或在 `config.yaml` 中设置 `outdir`。

---

## 11. 引用

如果 ProPhosLinker 对您的研究有帮助，请引用：

> *论文引用信息将在正式发表后更新。*

**GitHub:** [https://github.com/mzlab.research/ProPhosLinker](https://github.com/mzlab.research/ProPhosLinker)

**Docker Hub:** `ppy1222/prophoslinker`

**许可证:** MIT License

---

*如有问题或建议，请提交 GitHub Issue 或联系维护团队。*