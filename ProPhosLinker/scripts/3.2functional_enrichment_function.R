#!/usr/bin/env Rscript

#===============================================================================
# Functional Enrichment Analysis for Multi-omics Data
#===============================================================================
# This script performs comprehensive functional enrichment analysis comparing
# proteomics and phosphoproteomics datasets. It includes:
#   - GO (Gene Ontology) enrichment analysis with comparison plots
#   - KEGG pathway enrichment analysis with comparison plots
#   - Support for both online and local KEGG databases
#   - Customizable color schemes and visualization parameters
#===============================================================================

# Load required libraries with suppressed startup messages
options(warn = -1)
suppressWarnings(suppressPackageStartupMessages(library(clusterProfiler)))  
suppressWarnings(suppressPackageStartupMessages(library(org.Hs.eg.db)))     
suppressWarnings(suppressPackageStartupMessages(library(enrichplot)))       
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))         
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))           

# Increase timeout for downloading KEGG data (if needed)
options(timeout = 1000)

#===============================================================================
# DATA PREPARATION FUNCTION
#===============================================================================

#' Prepare data for enrichment analysis from differential expression results
#'
#' This function extracts lists of differentially expressed proteins and
#' phosphoproteins based on fold change and adjusted p-value thresholds.
#' It also maps phosphosites back to their corresponding proteins.
#'
#' @param pro_diff_path Path to proteomics differential expression results (TSV)
#' @param phos_diff_path Path to phosphoproteomics differential expression results (TSV)
#' @param phos_pro_path Path to protein-phosphosite mapping file (TSV)
#' @param outdir Output directory (default: "./")
#' @param log2FC Log2 fold change threshold (default: 1.2)
#' @param p_adj Adjusted p-value threshold (default: 0.05)
#' @return List containing two vectors: protein_list and phosphoprotein_list
enrichment_predata <- function(pro_diff_path, phos_diff_path, phos_pro_path, 
                               outdir = "./", log2FC = 1.2, p_adj = 0.05) {
  
  # Load input files
  pro_diff <- read.csv(pro_diff_path, sep = '\t')
  phos_diff <- read.csv(phos_diff_path, sep = '\t')
  phos_pro <- read.csv(phos_pro_path, sep = '\t')
  
  # Filter differentially expressed proteins based on thresholds
  pro_diff_filter_list <- pro_diff[
    ((pro_diff$logFC >= log2FC | pro_diff$logFC <= -log2FC) & 
       pro_diff$adj.P.Val < p_adj), 
  ]$Protein
  
  # Filter differentially expressed phosphosites
  phos_diff_filter_list <- phos_diff[
    ((phos_diff$logFC >= log2FC | phos_diff$logFC <= -log2FC) & 
       phos_diff$adj.P.Val < p_adj), 
  ]$Protein
  
  # Map phosphosites back to their corresponding proteins
  phos_diff_filter_list <- phos_pro[
    phos_pro$Site %in% phos_diff_filter_list, 
  ]$Protein
  
  return(list(
    "protein_list" = pro_diff_filter_list,
    "phosphoprotein_list" = phos_diff_filter_list
  ))
}

#===============================================================================
# GO ENRICHMENT COMPARISON FUNCTION
#===============================================================================

#' Perform GO enrichment analysis and save results for two omics datasets
#'
#' This function performs comparative GO enrichment analysis on two lists of genes from different omics datasets,
#' saves the results to a TSV file, and generates a dot plot visualizing the enrichment results.
#'
#' @param omics1_list List of gene IDs for the first omics dataset.
#' @param omics2_list List of gene IDs for the second omics dataset.
#' @param omics1_name Name of the first omics dataset (used in output filenames and legends).
#' @param omics2_name Name of the second omics dataset (used in output filenames and legends).
#' @param outdir Directory where the output files will be saved.
#' @param pvalueCutoff P-value cutoff for significance (default: 0.05).
#' @param showCategory Number of top categories to display per ontology (default: 15).
#' @param enrich_fromType Key type for gene annotation (default: 'UNIPROT').
#' @param color_gradient_low Color for low -log10(p-adjust) values (default: "#175663").
#' @param color_gradient_high Color for high -log10(p-adjust) values (default: "#90362d").

two_omicses_GO_enrichment <- function(omics1_list, omics2_list,
                                      omics1_name, omics2_name,
                                      outdir = "./",
                                      pvalueCutoff = 0.05, showCategory = 15,
                                      enrich_fromType = 'UNIPROT',
                                      color_gradient_low = "#175663",
                                      color_gradient_high = "#90362d",
                                      title_size = 16,
                                      axis_text_size = 10,
                                      axis_title_size = 10,
                                      legend_title_size = 10,
                                      legend_text_size = 8,
                                      strip_text_size = 12) {
  
  # Create a named list of gene clusters for comparison
  compare_lists <- list(omics1_list, omics2_list)
  names(compare_lists) <- c(omics1_name, omics2_name)
  
  # Perform comparative GO enrichment analysis
  go_compare <- compareCluster(
    geneClusters = compare_lists,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = enrich_fromType,
    ont = "ALL",  # Include all three GO ontologies: BP, MF, CC
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = 0.2
  )
  
  # Extract and save enrichment results
  go_compare_data <- as.data.frame(go_compare)
  
  write.table(
    go_compare_data,
    file = paste0(outdir, "/", omics1_name, "_", omics2_name, "_GO_enrichment.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  #---------------------------------------------------------------------------
  # Data preprocessing for visualization
  #---------------------------------------------------------------------------
  
  go_compare_clean <- go_compare_data %>%
    # Calculate numeric GeneRatio and BgRatio for sorting
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"),
                             function(x) as.numeric(x[1]) / as.numeric(x[2])),
      BgRatio_num = sapply(strsplit(BgRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2])),
      FoldEnrichment = GeneRatio_num / BgRatio_num,
      log_p_adjust = -log10(p.adjust)  # Add -log10(p.adjust) column
    ) %>%
    # Group by ontology and sort by adjusted p-value
    group_by(ONTOLOGY) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    # Keep top categories per ontology
    slice_head(n = showCategory) %>%
    ungroup() %>%
    # Set factor levels for proper plotting order
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC"))
    ) %>%
    arrange(desc(GeneRatio_num), .by_group = TRUE) %>%
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC")),
      GeneRatio_factor = as.factor(GeneRatio_num)  # Optional factor version
    )
  if(is.null(go_compare_clean) || nrow(go_compare_clean) == 0) {
    message("No GO enrichment results to visualize.")
    return(NULL)
  }
  
  #---------------------------------------------------------------------------
  # Create GO enrichment dot plot
  #---------------------------------------------------------------------------
  go_compare_GO <- ggplot(go_compare_clean, aes(x = Cluster, y = Description)) +
    geom_point(aes(size = Count, color = log_p_adjust), alpha = 0.8) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient(low = color_gradient_low, high = color_gradient_high,
                         name = "-log10(Adjusted P-value)") +
    scale_size_continuous(range = c(2, 8),
                          name = paste0(omics1_name, "/", omics2_name, " Count"),
                          guide = guide_legend(override.aes = list(color = "grey60"))) +
    labs(x = "Protein Ratio",
         y = "",
         title = paste0("GO Enrichment Comparison: ", omics1_name, " vs ", omics2_name)) +
    
    theme_bw() +
    theme(
      # Title settings 
      plot.title = element_text(hjust = 0.5, face = "bold", size = title_size,
                                margin = margin(b = 20)),
      
      # Facet strip settings 
      strip.text = element_text(size = strip_text_size, face = "bold",
                                margin = margin(t = 5, r = 5, b = 5, l = 5)),
      strip.background = element_rect(fill = "lightgray"),
      
      # Axis settings 
      axis.text.x = element_text(size = axis_text_size, color = "grey20", face = "bold"),
      axis.text.y = element_text(size = axis_text_size, color = "grey20", face = "bold"),
      axis.title.x = element_text(size = axis_title_size, face = "bold",
                                  margin = margin(t = 15)),
      
      # Legend settings 
      legend.title = element_text(size = legend_title_size, face = "bold"),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(0.8, "cm"),
      
      # Grid and margins
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    guides(
      fill = guide_colorbar(barheight = unit(4, "cm")),
      shape = guide_legend(
        override.aes = list(
          fill = "grey80",
          size = 6,
          stroke = 0.5
        )
      )
    )
  
  # Save plot to file
  ggsave(paste0(outdir, "/", omics1_name, "_", omics2_name, "_GO_enrichment.png"),
         go_compare_GO, width = 12, height = 8, dpi = 300)
}



#' Perform KEGG enrichment analysis and save results for two omics datasets
#'
#' This function performs comparative KEGG enrichment analysis on two lists of genes from different omics datasets,
#' saves the results to a TSV file, and generates a dot plot visualizing the enrichment results.
#'
#' @param omics1_list List of gene IDs for the first omics dataset.
#' @param omics2_list List of gene IDs for the second omics dataset.
#' @param omics1_name Name of the first omics dataset (used in output filenames and legends).
#' @param omics2_name Name of the second omics dataset (used in output filenames and legends).
#' @param outdir Directory where the output files will be saved.
#' @param pvalueCutoff P-value cutoff for significance (default: 0.05).
#' @param showCategory Number of top categories to display (default: 15).
#' @param enrich_fromType Key type for gene annotation (default: 'UNIPROT').
#' @param color_gradient_low Color for low -log10(p-adjust) values (default: "#175663").
#' @param color_gradient_high Color for high -log10(p-adjust) values (default: "#90362d").

two_omicses_KEGG_enrichment <- function(omics1_list, omics2_list,
                                        omics1_name, omics2_name,
                                        outdir = "./",
                                        pvalueCutoff = 0.05, showCategory = 15,
                                        enrich_fromType = 'UNIPROT',
                                        color_gradient_low = "#175663",
                                        color_gradient_high = "#90362d",
                                        title_size = 16,
                                        axis_text_size = 10,
                                        axis_title_size = 10,
                                        legend_title_size = 10,
                                        legend_text_size = 8) {
  
  # Create named list of gene clusters
  compare_lists <- list(omics1_list, omics2_list)
  names(compare_lists) <- c(omics1_name, omics2_name)
  
  # Handle different ID types for KEGG enrichment
  if (enrich_fromType == 'UNIPROT') {
    # Direct KEGG enrichment with UniProt IDs
    kk_compare <- compareCluster(
      geneClusters = compare_lists,
      fun = "enrichKEGG",
      organism = "hsa",
      keyType = "uniprot",
      pvalueCutoff = pvalueCutoff
    )
  } else {
    # Convert SYMBOL to ENTREZID for KEGG enrichment
    omics1_list_entrez <- bitr(
      omics1_list,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )$ENTREZID
    
    omics2_list_entrez <- bitr(
      omics2_list,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )$ENTREZID
    
    # Update comparison lists with ENTREZ IDs
    compare_lists <- list(omics1_list_entrez, omics2_list_entrez)
    names(compare_lists) <- c(omics1_name, omics2_name)
    
    kk_compare <- compareCluster(
      geneClusters = compare_lists,
      fun = "enrichKEGG",
      organism = "hsa",
      pvalueCutoff = pvalueCutoff
    )
  }
  
  # Extract and save enrichment results
  kk_compare_data <- as.data.frame(kk_compare)
  write.table(
    kk_compare_data,
    file = paste0(outdir, "/", omics1_name, "_", omics2_name, "_KEGG_enrichment.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  #---------------------------------------------------------------------------
  # Data preprocessing for visualization
  #---------------------------------------------------------------------------

  if(is.null(kk_compare_data) || nrow(kk_compare_data) == 0) {
    message("No KEGG enrichment results to visualize.")
    return(NULL)
  }
  
  kk_compare_data_clean <- kk_compare_data %>%
    # Calculate numeric GeneRatio
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"),
                             function(x) as.numeric(x[1]) / as.numeric(x[2])),
      BgRatio_num = sapply(strsplit(BgRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2])),
      FoldEnrichment = GeneRatio_num / BgRatio_num,
      log_p_adjust = -log10(p.adjust)  # Add -log10(p.adjust) column
    ) %>%
    # Sort and select top categories
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = showCategory) %>%
    arrange(desc(Count), .by_group = TRUE) %>%
    # Set factor levels for proper ordering
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      GeneRatio_factor = as.factor(GeneRatio_num)
    )
 
  #---------------------------------------------------------------------------
  # Create KEGG enrichment dot plot
  #---------------------------------------------------------------------------
  p_kk_compare <- ggplot(kk_compare_data_clean, aes(x = Cluster, y = Description)) +
    
    # Add points with size = gene count, color = -log10(adjusted p-value)
    geom_point(aes(size = Count, color = log_p_adjust), alpha = 0.8) +
    
    # Color gradient for -log10(p-values)
    scale_color_gradient(low = color_gradient_low, high = color_gradient_high,
                         name = "-log10(Adjusted P-value)") +
    
    # Size scale for points
    scale_size_continuous(range = c(2, 8),
                          name = paste0(omics1_name, " / ", omics2_name, " Count"),
                          guide = guide_legend(
                            override.aes = list(color = "grey60")
                          )) +
    
    # Labels and title
    labs(x = "Protein Ratio",
         y = "",
         title = paste0("KEGG Enrichment Comparison: ", omics1_name, " vs ", omics2_name)) +
    
    # Theme customization - 使用传入的参数
    theme_bw() +
    theme(
      # Title settings
      plot.title = element_text(hjust = 0.5, face = "bold", size = title_size,
                                margin = margin(b = 20)),
      
      # Axis settings
      axis.text.x = element_text(size = axis_text_size, color = "grey20", face = "bold"),
      axis.text.y = element_text(size = axis_text_size, color = "grey20", face = "bold"),
      axis.title.x = element_text(size = axis_title_size, face = "bold",
                                  margin = margin(t = 15)),
      
      # Legend settings
      legend.title = element_text(size = legend_title_size, face = "bold"),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(0.8, "cm"),
      
      # Grid and margins
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    
    # Legend customization
    guides(
      fill = guide_colorbar(barheight = unit(4, "cm")),
      shape = guide_legend(
        override.aes = list(
          fill = "grey80",
          size = 6,
          stroke = 0.5
        )
      )
    )
  
  # Save plot to file
  ggsave(paste0(outdir, "/", omics1_name, "_", omics2_name, "_KEGG_enrichment.png"),
         p_kk_compare, width = 12, height = 8, dpi = 300)
}



#===============================================================================
# MAIN ENRICHMENT PIPELINE FUNCTION
#===============================================================================

#' Complete enrichment analysis pipeline for multi-omics comparison
#'
#' This function orchestrates the entire enrichment analysis workflow:
#'   - Creates output directory if needed
#'   - Performs GO enrichment comparison (if enough genes in both lists)
#'   - Performs KEGG enrichment comparison (if enough genes in both lists)
#'
#' @param omics1_list Vector of feature IDs for first omics dataset
#' @param omics2_list Vector of feature IDs for second omics dataset
#' @param outdir Output directory
#' @param omics1_name Name for first omics dataset (default: "Proteomics")
#' @param omics2_name Name for second omics dataset (default: "Phosphoproteomics")
#' @param pvalueCutoff P-value cutoff (default: 0.05)
#' @param GO_showCategory Number of top GO terms to show (default: 6)
#' @param KEGG_showCategory Number of top KEGG pathways to show (default: 15)
#' @param enrich_fromType ID type ('UNIPROT' or 'SYMBOL', default: 'UNIPROT')
#' @param color_gradient_low Low end color (default: "#175663")
#' @param color_gradient_high High end color (default: "#90362d")
omics_enrichment_list <- function(omics1_list, omics2_list, outdir = "./",
                                  omics1_name = "Proteomics", 
                                  omics2_name = "Phosphoproteomics",
                                  pvalueCutoff = 0.05, 
                                  GO_showCategory = 6, 
                                  KEGG_showCategory = 15,
                                  enrich_fromType = 'UNIPROT',
                                  color_gradient_low = "#175663",
                                  color_gradient_high = "#90362d",
                                  title_size = 16,           
                                  axis_text_size = 10,       
                                  axis_title_size = 10,      
                                  legend_title_size = 10,    
                                  legend_text_size = 8,      
                                  strip_text_size = 12) {    
  

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # GO Enrichment
  if (length(omics1_list) >= 1 & length(omics2_list) >= 1) {
    two_omicses_GO_enrichment(
      omics1_list = omics1_list,
      omics2_list = omics2_list,
      omics1_name = omics1_name,
      omics2_name = omics2_name,
      outdir = outdir,
      pvalueCutoff = 0.05,
      showCategory = GO_showCategory,
      enrich_fromType = enrich_fromType,
      color_gradient_low = color_gradient_low,
      color_gradient_high = color_gradient_high,
      title_size = title_size,
      axis_text_size = axis_text_size,
      axis_title_size = axis_title_size,
      legend_title_size = legend_title_size,
      legend_text_size = legend_text_size,
      strip_text_size = strip_text_size
    )
  }
  
  # KEGG Enrichment
  if (length(omics1_list) >= 2 & length(omics2_list) >= 2) {
    two_omicses_KEGG_enrichment(
      omics1_list = omics1_list,
      omics2_list = omics2_list,
      omics1_name = omics1_name,
      omics2_name = omics2_name,
      outdir = outdir,
      pvalueCutoff = 0.05,
      showCategory = KEGG_showCategory,
      enrich_fromType = enrich_fromType,
      color_gradient_low = color_gradient_low,
      color_gradient_high = color_gradient_high,
      title_size = title_size,
      axis_text_size = axis_text_size,
      axis_title_size = axis_title_size,
      legend_title_size = legend_title_size,
      legend_text_size = legend_text_size
    )
  }
}