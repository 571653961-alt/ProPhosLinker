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

## 🎨 Conditional-response curve figures (step 2.2)

`scripts/2.2phosRate_subtyping.R` draws one representative site per response pattern
(`conditional_<GENE>_<INDEX>.png`). The layout is fixed — bold title on top left,
legend below it on the left, three interval background panels, joint-model confidence
bands, and a bold grey interval-statistics caption at the bottom left — while every
visual constant is exposed as a parameter.

Two equivalent entry points:

1. **`config.yaml`** → `phosrate_subtyping` → `plot` (used by `main.py` / the CLI pipeline).
2. **R command line** → `--plot-*` options (defaults below; `Rscript scripts/2.2phosRate_subtyping.R --help`).

| `config.yaml` (`plot:`) | R option | Default | Meaning |
| :--- | :--- | :--- | :--- |
| `width` / `height` / `dpi` | `--plot-width` / `--plot-height` / `--plot-dpi` | `8` / `6` / `300` | Figure size in inches and resolution |
| `base_size` | `--plot-base-size` | `13` | `theme_classic()` base size |
| `title_size` | `--plot-title-size` | `22` | Bold title, left aligned to the panel |
| `axis_title_size` | `--plot-axis-title-size` | `15` | Bold axis titles |
| `axis_text_size` | `--plot-axis-text-size` | `13` | Bold tick labels |
| `legend_size` | `--plot-legend-size` | `14` | Legend on top, left aligned (`N` / `T`) |
| `caption_size` | `--plot-caption-size` | `11` | `Low/Middle/High` statistics caption (plain, not bold) |
| `caption_color` | `--plot-caption-color` | `grey30` | Caption text color |
| `line_width` | `--plot-line-width` | `2.0` | Fitted-curve line width |
| `point_size` / `point_alpha` | `--plot-point-size` / `--plot-point-alpha` | `1.8` / `0.55` | Scatter points |
| `ribbon_alpha` | `--plot-ribbon-alpha` | `0.15` | Confidence band around each fitted curve |
| `spine_width` / `tick_width` / `tick_length` | `--plot-spine-width` / `--plot-tick-width` / `--plot-tick-length` | `1.2` / `1.2` / `4` | Axis lines (left + bottom only) and ticks |
| `dash_color` / `dash_width` | `--plot-dash-color` / `--plot-dash-width` | `grey45` / `1.1` | Low/middle/high interval boundaries |
| `panel_low` / `panel_mid` / `panel_high` | `--plot-panel-low` / `--plot-panel-mid` / `--plot-panel-high` | `#F7F9FC` / `#FFF9EF` / `#F9F5FA` | Interval background panels |
| `group_labels` | `--plot-group-labels` | `N,T` | Legend labels for Control, Case |
| `n_representative` | `--plot-n-representative` | `5` | Sites per global class (UP/DOWN/NS) for representative plots; requires `--global-phosprotein-diff` |

Group colors are shared with the routing figures: `--color-control` (default `#012da7`)
and `--color-case` (default `#ff7f00`). When the pipeline is used, the same values are
passed through `DifferentialConfig` (`plot_*` fields in `config/analysis_config.py`).

---

## 📄 License

Released under the terms of the `LICENSE` file.