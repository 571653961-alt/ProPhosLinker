# Case Study: Pancreatic Ductal Adenocarcinoma (PDAC)

This directory contains a **curated subset** of results from running ProPhosLinker on a PDAC phosphoproteomics dataset (215 samples: 140 tumor, 75 normal). The full results (320 MB) are too large for the repository; the files here are representative and sufficient to verify the pipeline output.

## Data Source

Input data files: `../../casedata/` (119 MB, 5 files)
- `protein_abundance.tsv` — Protein-level expression matrix (8,754 proteins × 215 samples)
- `phosphoprotein_abundance.tsv` — Phosphosite-level expression matrix (16,651 sites × 215 samples)
- `compare_groups.tsv` — Group comparison design (Tumor vs Normal)
- `protein_phosphoproSite.tsv` — Protein-phosphosite mapping table (51,399 pairs)
- `clinical_table_140.tsv` — Clinical metadata for 140 tumor samples

## Results by Pipeline Step

### 0. Data Preprocessing
| File | Description |
|------|-------------|
| `input_overview.txt` | Input data summary report |
| `mapping_list.tsv` | Protein-phosphosite mapping (51,399 pairs) |
| `metadata.tsv` | Formatted clinical metadata |
| `sample_list.tsv` | Sample group assignment (T=140, N=75) |

### 1. Pattern Analysis
| File | Description |
|------|-------------|
| `molecular_subtyping_clusters.tsv` | Cluster assignments for molecular subtypes |
| `molecular_subtyping_sankey_by_group.png` | Sankey diagram: protein ↔ phospho subtypes by group |
| `WGCNA_softpower.tsv` | WGCNA soft-thresholding power analysis |
| `module_correlation_filtering.tsv` | Module-trait correlation results |
| `module_correlation_pheatmap.png` | Module-trait correlation heatmap |

### 2. Differential Analysis
| File | Description |
|------|-------------|
| `T_vs_N_Protein.tsv` | Protein-level differential expression (Tumor vs Normal) |
| `T_vs_N_PhosProtein.tsv` | Phosphosite-level differential expression |
| `T_vs_N_combined.tsv` | Combined protein + phosphosite differential results |
| `cor_value_results.tsv` | Protein-phosphosite correlation values |
| `headline.tsv` | Phosphorylation rate subtyping headline statistics |
| `site_routing.tsv` | Site-level routing classification |
| `network_by_top_modules.png` | Differential network clustering by top modules |

### 3. Functional Analysis
| File | Description |
|------|-------------|
| `enrichment/Protein_PhosProtein_GO_enrichment.tsv` | GO enrichment for integrated protein + phosphosite |
| `enrichment/Protein_PhosProtein_KEGG_enrichment.tsv` | KEGG pathway enrichment |
| `Node_info.tsv` | Knowledge graph node information |
| `Edge_info.tsv` | Knowledge graph edge information (PPI + KS + LR + TF) |
| `community_detection/pro_phosphopro_network_ALL.png` | Community detection module 3 network (representative) |

## Reproduction

To reproduce the full analysis:

```bash
prophoslinker \
  --pro_file "casedata/protein_abundance.tsv" \
  --phos_file "casedata/phosphoprotein_abundance.tsv" \
  --sample_group "casedata/compare_groups.tsv" \
  --mapping_file "casedata/protein_phosphoproSite.tsv" \
  --outdir "results" \
  --group_comparing "T:N" \
  --password "12345678" \
  --metadata_file "casedata/clinical_table_140.tsv" \
  --disease "pancreatic ductal adenocarcinoma"
```

> **Note:** The `--password` parameter above is for Neo4j authentication. Use a strong password in production.

See the [User Manual](../user_manual.md#9-case-study) for a detailed walkthrough of each result.