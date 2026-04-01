# ProPhosLinker

ProPhosLinker is an integrative analysis toolkit designed to provide an in-depth understanding of proteomics and phosphoproteomics data. It facilitates the exploration of complex biological processes by offering a reproducible pipeline for:

- **Differential expression analysis** of proteins and phosphosites to uncover key molecular differences
- **Phosphorylation rate inference** (site-level vs protein-level) through advanced quantile modelling
- **Knowledge graph-based functional network construction** (Neo4j-enabled), leveraging heterogeneous data sources (e.g., PhosphoSitePlus, TRRUST)
- **Subtype concordance and clustering** using advanced methods (Mfuzz, WGCNA, NMF)

---

## 🚀 Key Features

- **Comprehensive Analysis**: Perform differential expression analysis, functional enrichment, and network analysis, all in a unified pipeline.
- **Innovative Modelling**: Quantile modelling to assess phosphorylation rate variations incorporating biological context.
- **Knowledge Graph Integration**: Construct and analyze dynamic functional networks using Neo4j.
- **Flexible & Reproducible**: Easily run specific pipeline steps or adjust parameters via CLI and YAML configuration.

---

## 📦 Installation & Usage

To use **ProPhosLinker**, you can choose the highly recommended Docker method (which avoids dependency issues) or install everything manually.

### 🐳 Method 1: Docker (Highly Recommended)

Using Docker is the easiest way to run ProPhosLinker as it bundles Python, R, all required packages, and Neo4j into a single environment.

#### 1. Build the Docker Image
Navigate to the root directory of the project and run:
```bash
docker build -t prophoslinker:v1 .
```

#### 2. Run and Enter the Container
To run the container interactively and access the results on your local machine, use volume mounting (`-v`):
```bash
docker run -it --rm \
  -p 7474:7474 -p 7687:7687 \
  -v /path/to/your/local/results:/app/results \
  --name prophoslinker_test \
  prophoslinker:v1
```
> 💡 **Note**: This command forwards Neo4j ports so you can view the knowledge graph at `http://localhost:7474` in your browser.

#### 3. Run the Analysis (Inside the Container)
Once you are inside the running container's terminal, you can execute the full analysis pipeline directly using the provided case data:
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

If you prefer to set up the environment on your native machine, follow these steps in order:

#### 1. Install Neo4j (>= 5.x)
ProPhosLinker relies on Neo4j for functional network analysis. 
1. Download Neo4j from the [official website](https://neo4j.com/download/).
2. Place the `ProPhosLinker_KG.dump` file in a directory accessible by Neo4j.
3. Run the following command to import the knowledge graph:
```bash
neo4j-admin database load neo4j --from-path=./ProPhosLinker/database --overwrite-destination=true
```

#### 2. Install R & Required Packages (R >= 4.0)
Open your R console and run the following script to install CRAN and Bioconductor packages with optimized mirrors:

```r
# --- 1. Base Tools & Environment ---
options(repos = c(CRAN = "[https://mirrors.tuna.tsinghua.edu.cn/CRAN/](https://mirrors.tuna.tsinghua.edu.cn/CRAN/)"))
options(BioC_mirror = "[https://mirrors.tuna.tsinghua.edu.cn/bioconductor](https://mirrors.tuna.tsinghua.edu.cn/bioconductor)")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# --- 2. Bioconductor Packages ---
bioc_packages <- c(
  "AnnotationDbi", "GO.db", "preprocessCore", "impute", 
  "limma", "Mfuzz", "GOSemSim", "DOSE", "qvalue",
  "ggtree", "enrichplot", "clusterProfiler", "org.Hs.eg.db"
)
BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

# --- 3. CRAN Packages ---
cran_packages <- c(
  "optparse", "readr", "stringr", "dplyr", "vegan", "ggrepel", "ggplot2", 
  "NMF", "tidyverse", "doParallel", "WGCNA", "patchwork", "pheatmap", 
  "plyr", "viridis", "flashClust", "colorspace", "reshape2", "tibble", 
  "igraph", "ggraph", "tidygraph", "tidyr", "ggforce", "ggpubr", 
  "Hmisc", "bootnet", "graphlayouts", "scatterpie", "ggsci", 
  "ggnewscale", "svglite", "ggiraph", "remotes", "ggsankeyfier",
  "units", "gdtools", "systemfonts"
)

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = c("Depends", "Imports", "LinkingTo"), upgrade = "never")
  }
}
```

#### 3. Install Python Dependencies
Finally, install the Python side of the toolkit:
```bash
pip install -r requirements.txt
pip install -e .
```

---

## 🛠️ Quick Start (CLI - Manual Mode)

If you installed ProPhosLinker manually, you can run it with additional parameter configurations:

### Full pipeline execution example
```bash
ProPhosLinker \
  --pro_file "casedata/protein_abundance.tsv" \
  --phos_file "casedata/phosphoprotein_abundance.tsv" \
  --sample_group "casedata/compare_groups.tsv" \
  --mapping_file "casedata/protein_phosphoproSite.tsv" \
  --metadata_file "casedata/clinical_table_140.tsv" \
  --config "ProPhosLinker/config.yaml" \
  --group_comparing "T:N" \
  --outdir "." \
  --pro_FC 2 \
  --phos_FC 2 \
  --network_FC 2 \
  --password "neo4j_password"
```

### Run only specific steps
```bash
ProPhosLinker \
  --pro_file protein_abundance.tsv \
  --phos_file protein_phosphoproSite.tsv \
  --sample_group compare_groups.tsv \
  --mapping_file mapping.tsv \
  --group_comparing "Case:Control" \
  --steps "data_preprocessing,differential_analysis"
```

---

## 📌 Supported Pipeline Steps (`--steps`)

| Step name | Description |
| :--- | :--- |
| `data_preprocessing` | Quality control, median centering normalization, and missing value imputation. |
| `pattern_analysis` | Multiscale integration via Procrustes analysis, WGCNA, and NMF-based subtyping. |
| `differential_analysis`| Abundance-dependent phosphorylation dynamics (DPR) and differential correlation networks (DCN). |
| `functional_analysis` | Mapping signals to Neo4j knowledge graphs and running network biology functional integration. |

> 💡 **Note**: The pipeline always runs in its internal logical order regardless of the order you specify in `--steps`.

---

## 📄 License

This project is released under the terms of the `LICENSE` file.
