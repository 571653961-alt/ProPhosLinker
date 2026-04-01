# ProPhosLinker

ProPhosLinker is an integrative analysis toolkit designed to provide an in-depth understanding of proteomics and phosphoproteomics data. It facilitates the exploration of complex biological processes by offering a reproducible pipeline for:

- **Differential expression analysis** of proteins and phosphosites.
- **Phosphorylation rate inference** (site-level vs protein-level) through advanced quantile modelling.
- **Knowledge graph-based functional network construction** (Neo4j-enabled), leveraging heterogeneous data sources (e.g., PhosphoSitePlus, TRRUST).
- **Subtype concordance and clustering** using advanced methods (Mfuzz, WGCNA, NMF).

---

## 🚀 Key Features

- **Comprehensive Analysis**: Unified pipeline for DE analysis, enrichment, and network construction.
- **Innovative Modelling**: Quantile modelling to assess site-specific phosphorylation dynamics.
- **Knowledge Graph Integration**: Integrated support for Neo4j-based functional network biology.
- **Flexible & Reproducible**: Fully controllable via CLI or YAML configuration.

---

## 📦 Installation

### 🐳 Method 1: Docker (Highly Recommended)

This is the easiest way to run ProPhosLinker. The pre-built image contains all necessary R and Python dependencies, as well as a pre-configured environment.

#### 1. Pull the Image
```bash
docker pull ppy1222/prophoslinker:v1
```

#### 2. Run the Container
Run the container and mount your local directory to save analysis results (replace `/your/local/path` with your actual folder):
```bash
docker run -it --rm \
  -p 7474:7474 -p 7687:7687 \
  -v /your/local/path:/app/results \
  --name prophoslinker_test \
  ppy1222/prophoslinker:v1
```

#### 3. Run analysis inside Docker
Once inside the container terminal, you can run the pipeline immediately:
```bash
ProPhosLinker \
  --pro_file "casedata/protein_abundance.tsv" \
  --phos_file "casedata/phosphoprotein_abundance.tsv" \
  --sample_group "casedata/compare_groups.tsv" \
  --mapping_file "casedata/protein_phosphoproSite.tsv" \
  --metadata_file "casedata/clinical_table_140.tsv"
```

---

### 🛠️ Method 2: Manual Installation

#### 1. Install Neo4j (>= 5.x)
1. Download from [official website](https://neo4j.com/download/).
2. Import the knowledge graph:
```bash
neo4j-admin database load neo4j --from-path=./ProPhosLinker/database --overwrite-destination=true
```

#### 2. Install R & Required Packages
Run this in your R console (ensure R >= 4.0):
```r
# Set mirrors
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
  "AnnotationDbi", "GO.db", "preprocessCore", "impute", "limma", 
  "Mfuzz", "GOSemSim", "DOSE", "qvalue", "ggtree", "enrichplot", 
  "clusterProfiler", "org.Hs.eg.db", "NMF", "WGCNA"
), ask = FALSE, update = FALSE)

# Install CRAN packages
cran_pkgs <- c(
  "optparse", "readr", "stringr", "dplyr", "vegan", "ggrepel", "ggplot2", 
  "tidyverse", "doParallel", "patchwork", "pheatmap", "plyr", "viridis", 
  "flashClust", "ggsankeyfier", "colorspace", "reshape2", "tibble", 
  "igraph", "ggraph", "tidygraph", "tidyr", "ggforce", "ggpubr", "Hmisc", 
  "bootnet", "graphlayouts", "scatterpie", "ggsci", "ggnewscale", 
  "svglite", "ggiraph", "units", "gdtools", "systemfonts"
)
install.packages(cran_pkgs, dependencies = c("Depends", "Imports", "LinkingTo"))
```

#### 3. Install Python Dependencies
```bash
pip install -r requirements.txt
pip install -e .
```

---

## 📌 Supported Pipeline Steps (`--steps`)

| Step name | Description |
| :--- | :--- |
| `data_preprocessing` | QC, normalization, and missing value imputation. |
| `pattern_analysis` | Multiscale integration (Procrustes, WGCNA, NMF subtyping). |
| `differential_analysis`| Statistical association, DPR modeling, and DCN inference. |
| `functional_analysis` | Neo4j-based functional network integration. |

---

## 📄 License

Released under the terms of the `LICENSE` file.