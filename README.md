# ProPhosLinker

ProPhosLinker is an integrative analysis toolkit designed to provide an in-depth understanding of proteomics and phosphoproteomics data. It facilitates the exploration of complex biological processes by offering a reproducible pipeline for:

- **Differential expression analysis** of proteins and phosphosites to uncover key molecular differences
- **Phosphorylation rate inference** (site-level vs protein-level) through advanced quantile modelling, emphasizing the nuanced relationship between protein abundance and phosphorylation dynamics
- **Knowledge graph-based functional network construction** (Neo4j-enabled), leveraging heterogeneous data sources such as PhosphoSitePlus, TRRUST, CellChatDB, and CKG for robust protein-site functional mapping
- **Subtype concordance and clustering** using advanced methods (Mfuzz, WGCNA, NMF), to identify molecular subtypes and their phosphorylation profiles


---

## 🚀 Key Features

- **Comprehensive analysis**: Perform differential expression analysis, functional enrichment, and network analysis, all in a unified pipeline
- **Innovative modelling**: Quantile modelling to assess phosphorylation rate variations at the site and protein level, incorporating biological context for precision
- **Advanced functional network analysis**: Integrate heterogeneous biological knowledge for constructing and analyzing dynamic functional networks using Neo4j
- **Flexible and reproducible**: Easily run specific pipeline steps, adjust parameters via CLI or YAML configuration, and ensure reproducibility with version-controlled workflows


---

## 📦 Installation

To use **ProPhosLinker**, you need to install dependencies in the following order: **Neo4j**, **R packages**, and **Python dependencies**.

### 1) Install Neo4j (required for functional network analysis)

**ProPhosLinker** relies on **Neo4j 4.x+** for functional network analysis. Ensure Neo4j is installed and accessible via `neo4j` or `cypher-shell`.

#### Download Neo4j:

- Download Neo4j from the [official website](https://neo4j.com/download/). Choose **Neo4j 5.x+** for compatibility with ProPhosLinker.

#### Import the Knowledge Graph:

1. After installing Neo4j, place the `ProPhosLinker_KG.dump` file in a directory accessible by Neo4j.
2. Run the following command to import the knowledge graph into Neo4j:

```bash
neo4j-admin load --database=neo4j --from=path/to/ProPhosLinker_KG.dump --force
````

Make sure to replace `path/to/ProPhosLinker_KG.dump` with the actual path of the dump file.

---

### 2) Install R and required R packages

ProPhosLinker relies on several R scripts for certain analysis steps. Please ensure **R (>= 4.0)** is installed, and then install the necessary R packages.

Run the following command in your R console to install the required packages:

```r
install.packages(c(
  "optparse", "readr", "stringr", "dplyr", "vegan", "ggrepel", "ggplot2", "NMF", 
  "tidyverse", "doParallel", "WGCNA", "patchwork", "pheatmap", "plyr", "viridis", 
  "grid", "flashClust", "ggsankeyfier", "limma", "statmod", "colorspace", 
  "Mfuzz", "reshape2", "tibble", "igraph", "ggraph", "tidygraph", "tidyr", 
  "ggforce", "ggpubr", "clusterProfiler", "org.Hs.eg.db", "enrichplot", 
  "hmisc", "bootnet", "graphlayouts", "scatterpie", "ggsci", "ggnewscale", 
  "svglite", "ggiraph"
))
```

Once the packages are installed, proceed to the next step.

---

### 3) Install Python dependencies

Now that the R environment is set up, you can install the required Python packages for ProPhosLinker.

1. Install the Python dependencies via the following command:

```bash
pip install -r requirements.txt
```

2. **Editable installation** (recommended for development):

```bash
pip install -e .
```

This will allow you to modify the code locally and have those changes immediately reflected in your environment.

---

By following these steps, you will have everything set up to run **ProPhosLinker** and use it for analysis.

---

## 🛠️ Quick Start (CLI)

### Basic example (full pipeline)

```bash
ProPhosLinker \
  --pro_file "casedata/protein_abundance.tsv" \
  --phos_file "casedata/phosphoprotein_abundance.tsv" \
  --sample_group "casedata/compare_groups.tsv" \
  --mapping_file "casedata/protein_phosphoproSite.tsv" \
  --metadata_file "casedata/clinical_table_140.tsv"
  --config "ProPhosLinker/config.yaml" \
  --group_comparing "T:N" \
  --outdir "." \
  --pro_FC 2 \
  --phos_FC 2 \
  --network_FC 2 \
  --password "neo4j_password"
```

### Run only a specific step

```bash
ProPhosLinker \
  --pro_file protein_abundance.tsv \
  --phos_file protein_phosphoproSite.tsv \
  --sample_group compare_groups.tsv \
  --mapping_file mapping.tsv \
  --config "ProPhosLinker/config.yaml" \
  --group_comparing "Case:Control" \
  --outdir "." \
  --pro_FC 1.5 \
  --phos_FC 1.5 \
  --network_FC 1.5 \
  --steps functional_analysis
```

### Run multiple steps

```bash
ProPhosLinker \
  --pro_file protein_abundance.tsv \
  --phos_file protein_phosphoproSite.tsv \
  --sample_group compare_groups.tsv \
  --mapping_file mapping.tsv \
  --config "ProPhosLinker/config.yaml" \
  --group_comparing "Case:Control" \
  --outdir "." \
  --pro_FC 1.8 \
  --phos_FC 1.8 \
  --network_FC 1.8 \
  --steps "data_preprocessing,differential_analysis,functional_analysis"
```

---

## 📌 Supported Pipeline Steps (`--steps`)
| Step name               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `data_preprocessing`    | **Data preprocessing** ensures robust statistical inference through quality control, performs median centering normalization, and completes missing value imputation. Prepares the data for downstream analysis..                                                                                                                                                                                                        |
| `pattern_analysis`      | **Multiscale integration** through quantification of global similarity between omics datasets (Procrustes analysis) and weighted gene co-expression network analysis (WGCNA), identifying key molecular modules supporting cross-omics associations. Visualize the consistency or disconnection between proteomics and phosphoproteomics subtypes (NMF-based subtyping).                                                                                                                                                                             |
| `differential_analysis` | **Statistical association analysis** designed to uncover complex differential regulation features. This module not only covers standard differential expression but also introduces differential correlation network (DCN) inference to identify correlation shifts among differentially expressed molecules. Additionally, **abundance-dependent phosphorylation dynamics modeling** (DPR) is used to precisely analyze the nonlinear modification patterns constrained by protein expression levels. |
| `functional_analysis`   | **Functional and knowledge-driven integration**, mapping multidimensional statistical signals to heterogeneous biological networks integrated with prior knowledge. This step uses Neo4j for knowledge graph-based functional network construction.                                                                                                                                                                                                                                                    |


> Note: The pipeline always runs in the internal logical order regardless of the order you specify in `--steps`.

---

## 🧩 Project Structure

```text
E:/项目/ProPhosLinker/
├── LICENSE
├── MANIFEST.in
├── pyproject.toml
├── README.md
├── requirements.txt
├── setup.py
│
├── casedata/                     # example input files
│   ├── clinical_table_140.tsv
│   ├── compare_groups.tsv
│   ├── compare_groups_old_V.tsv
│   ├── phosphoprotein_abundance.tsv
│   ├── phosprotein_diff.tsv
│   ├── protein_abundance.tsv
│   ├── protein_diff.tsv
│   └── protein_phosphoproSite.tsv
│
└── ProPhosLinker/                # main Python package
    ├── cli.py                    # command line interface entry point
    ├── config.yaml               # default configuration
    ├── main.py                   # main orchestration module
    ├── __init__.py
    │
    ├── analysis/                 # analysis modules (differential, network, etc.)
    │   ├── data_preprocessing.py
    │   ├── differential_analysis.py
    │   ├── functional_analysis.py
    │   ├── functionnal_interaction_network.py
    │   ├── hub_protein_visualization.py
    │   ├── pattern_analysis.py
    │   ├── phos_rate_subtyping.py
    │   ├── pro_phosphopro_network_analysis.py
    │   ├── pro_phosphopro_network_contructor.py
    │   ├── pro_phosphopro_network_visualization.py
    │   └── __init__.py
    │
    ├── config/                   # configuration templates / defaults
    │   ├── analysis_config.py
    │   ├── basic_config.py
    │   ├── color_config.py
    │   ├── neo4j_config.py
    │   ├── preprocessing_config.py
    │   ├── result_dir_config.py
    │   ├── script_config.py
    │   └── __init__.py
    │
    ├── database/                 # Neo4j-related helpers / loaders
    │   ├── neo4j_manager.py
    │   ├── __init__.py
    │   └── dataload/
    │
    ├── scripts/                  # R scripts called by pipeline
    │   ├── 1.1procrustes.R
    │   ├── 1.2clusterbyNMF.R
    │   ├── 1.3WGCNA.R
    │   ├── 1.pattern_analysis.R
    │   ├── 2.1omics_differential_analysis.R
    │   ├── 2.2phosRate_quantile_subtyping.R
    │   ├── 2.3differential_network.R
    │   ├── 3.1functional_enrichment.R
    │   ├── 3.2functional_enrichment_function.R
    │   ├── 3.2functionnal_interaction_network_community_visualisation.R
    │   ├── 3.2functionnal_interaction_network_visualization.R
    │   ├── 3.2functionnal_interaction_network_visualization_CN.R
    │   ├── 3.3network_hubgen_visualization.R
    │   └── 2.3Stable_Differential_Network/
    │       ├── .Rhistory
    │       ├── differential_network.R
    │       ├── differential_subnetwork_plot.R
    │       ├── diff_net_community_detection_plot.R
    │       ├── functional_enrichment_function_delete.R
    │       ├── network_show.R
    │       ├── network_show_delete.R
    │       ├── pipline_save.R
    │       ├── run_cluster.R
    │       ├── run_color.R
    │       ├── run_conditional_network.R
    │       ├── run_corStability.R
    │       ├── run_diff.R
    │       ├── run_diffsubnet_enrichment.R
    │       ├── run_diff_enrichment.R
    │       ├── run_diff_network.R
    │       ├── run_enrichment.R
    │       ├── run_mediation.R
    │       ├── run_predata.R
    │       ├── run_prenetwork.R
    │       └── run_samplelist.R
    │
    └── tools/                    # helper scripts/tools
        ├── neo4j_import.py
        └── __init__.py
```

---

## 💡 Tips & Notes

- **If you skip a step** (e.g., `differential_analysis`), be sure the downstream step has required intermediate outputs.
- **Neo4j must be running** before using `functional_analysis`.
- Use `--config <path>` to customize parameters (see `ProPhosLinker/config.yaml`).

---

## 📄 License

This project is released under the terms of the `LICENSE` file.
