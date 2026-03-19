#!/usr/bin/env Rscript

#===============================================================================
# Functional Enrichment Analysis for Multi-omics Data
#===============================================================================
# This script performs comparative functional enrichment analysis (GO and KEGG)
# between proteomics and phosphoproteomics datasets. It identifies significantly
# enriched biological processes, molecular functions, cellular components, and
# KEGG pathways, and generates comparative dot plots.
#
# Key features:
#   - Filter differentially expressed proteins and phosphosites
#   - Map phosphosites back to their corresponding proteins
#   - Perform GO enrichment comparison (BP, MF, CC)
#   - Perform KEGG pathway enrichment comparison
#   - Generate publication-quality dot plots with customizable colors
#===============================================================================

options(warn = -1)
suppressWarnings(library(optparse))

#===============================================================================
# COMMAND LINE ARGUMENT DEFINITION
#===============================================================================

option_list <- list(
  #---------------------------------------------------------------------------
  # Input Parameters
  #---------------------------------------------------------------------------
  make_option(c("--pro_diff"), type = "character", metavar = "FILE",
              help = "Differential protein expression file (TSV format)"),
  make_option(c("--phos_diff"), type = "character", metavar = "FILE",
              help = "Differential phosphoprotein expression file (TSV format)"),
  make_option(c("--phos_pro"), type = "character", metavar = "FILE",
              help = "Phosphosite to protein mapping file (TSV format)"),
  
  #---------------------------------------------------------------------------
  # Output Parameters
  #---------------------------------------------------------------------------
  make_option(c("--outdir"), type = "character", default = getwd(), metavar = "DIR",
              help = "Output directory for results [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Analysis Parameters
  #---------------------------------------------------------------------------
  make_option(c("--log2FC"), type = "numeric", default = 1.2, metavar = "NUM",
              help = "Log2 fold change cutoff for differential expression [default: %default]"),
  make_option(c("--diff_p_adj"), type = "numeric", default = 0.05, metavar = "NUM",
              help = "Adjusted p-value cutoff for differential expression significance [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Omics Names
  #---------------------------------------------------------------------------
  make_option(c("--omics1_name"), type = "character", default = "Proteomics", metavar = "STR",
              help = "Name for first omics dataset (proteomics) [default: %default]"),
  make_option(c("--omics2_name"), type = "character", default = "Phosphoproteomics", metavar = "STR",
              help = "Name for second omics dataset (phosphoproteomics) [default: %default]"),
  make_option(c("--enrich_fromType"), type = "character", default = "UNIPROT", metavar = "STR",
              help = "Input ID type for enrichment analysis (e.g., UNIPROT, SYMBOL, ENSEMBL, ENTREZ) [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Enrichment Parameters
  #---------------------------------------------------------------------------
  make_option(c("--GO_showCategory"), type = "integer", default = 6, metavar = "INT",
              help = "Number of GO terms to display per ontology [default: %default]"),
  make_option(c("--KEGG_showCategory"), type = "integer", default = 15, metavar = "INT",
              help = "Number of KEGG pathways to display [default: %default]"),
  make_option(c("--pvalueCutoff"), type = "numeric", default = 0.05, metavar = "NUM",
              help = "P-value cutoff for enrichment analysis [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Color Parameters
  #---------------------------------------------------------------------------
  make_option(c("--color_gradient_low"), type = "character", default = "#175663", metavar = "character",
              help = "Color for low end of enrichment gradient [default: %default]"),
  make_option(c("--color_gradient_high"), type = "character", default = "#90362d", metavar = "character",
              help = "Color for high end of enrichment gradient [default: %default]"),
  
  #---------------------------------------------------------------------------
  # General Parameters
  #---------------------------------------------------------------------------
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE, metavar = "FLAG",
              help = "Print detailed output messages [Default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                           description = "Perform comparative enrichment analysis (GO and KEGG) between protein and phosphoprotein datasets.")
opt <- parse_args(opt_parser)

#===============================================================================
# PARAMETER VALIDATION FUNCTION
#===============================================================================

#' Validate and clean input parameters
#'
#' This function performs comprehensive validation of all command-line parameters:
#'   - Cleans and normalizes file paths
#'   - Checks existence and non-emptiness of input files
#'   - Creates output directory if it doesn't exist
#'   - Validates numeric parameters within specified ranges
#'
#' @param opt List of parsed command-line arguments
#' @return Validated and cleaned parameter list
parameter_validation <- function(opt) {
  
  #---------------------------------------------------------------------------
  # Helper Functions
  #---------------------------------------------------------------------------
  
  #' Clean and normalize file paths
  #' @param path Raw path string
  #' @return Cleaned path with forward slashes
  clean_path <- function(path) {
    # Trim leading and trailing whitespace and quotes
    path <- gsub("^['\" ]+|['\" ]+$", "", path)
    # Replace backslashes with forward slashes for cross-platform compatibility
    path <- gsub("\\\\", "/", path)
    return(path)
  }
  
  #' Check if input file exists and is not empty
  #' @param file_path Path to file
  #' @param file_name Name of file for error messages
  check_input_file <- function(file_path, file_name) {
    if (!file.exists(file_path)) {
      stop("❌ Error: The file '", file_name, "' does not exist at the specified path: ", file_path)
    }
    if (file.size(file_path) == 0) {
      stop("❌ Error: The file '", file_name, "' is empty at the specified path: ", file_path)
    }
  }
  
  #' Create or validate output directory
  #' @param result_dir Output directory path
  #' @return Path to output directory
  check_output_dir <- function(result_dir) {
    if (!dir.exists(result_dir)) {
      # Recursively create directory and all parent directories
      dir.create(result_dir, recursive = TRUE)
      print(paste("✅ Output directory created:", result_dir))
    } else {
      print(paste("✅ Output directory already exists:", result_dir))
    }
    return(result_dir)
  }
  
  #' Validate categorical parameter against allowed values
  #' @param value Parameter value to check
  #' @param allowed_values Vector of allowed values
  #' @param option_name Name of option for error messages
  validate_option_choices <- function(value, allowed_values, option_name) {
    if (!value %in% allowed_values) {
      stop("❌ Error: Invalid", option_name, "value '", value, "'. Please use one of the following options: ",
           paste(allowed_values, collapse = ", "))
    }
  }
  
  #' Validate numeric parameter within specified range
  #' @param value Numeric value to check
  #' @param min_val Minimum allowed value
  #' @param max_val Maximum allowed value
  #' @param option_name Name of option for error messages
  validate_numeric_range <- function(value, min_val, max_val, option_name) {
    if (is.null(value)) {
      return()  # Skip validation if the value is NULL
    }
    if (value < min_val || value > max_val) {
      stop("❌ Error: ", option_name, " must be within the range [", min_val, ", ", max_val, "].")
    }
  }
  
  #---------------------------------------------------------------------------
  # Path Cleaning
  #---------------------------------------------------------------------------
  opt$pro_diff <- clean_path(opt$pro_diff)
  opt$phos_diff <- clean_path(opt$phos_diff)
  opt$phos_pro <- clean_path(opt$phos_pro)
  opt$outdir <- clean_path(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Input File Validation
  #---------------------------------------------------------------------------
  check_input_file(opt$pro_diff, "pro_diff")
  check_input_file(opt$phos_diff, "phos_diff")
  check_input_file(opt$phos_pro, "phos_pro")
  
  # Create output directory
  opt$outdir <- check_output_dir(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Numeric Parameter Validation
  #---------------------------------------------------------------------------
  validate_numeric_range(opt$log2FC, 0, 10, "log2FC")
  validate_numeric_range(opt$diff_p_adj, 0, 1, "diff_p_adj")
  validate_numeric_range(opt$GO_showCategory, 1, 20, "GO_showCategory")
  validate_numeric_range(opt$KEGG_showCategory, 1, 20, "KEGG_showCategory")
  validate_numeric_range(opt$pvalueCutoff, 0, 1, "pvalueCutoff")
  
  #---------------------------------------------------------------------------
  # Verbose Output
  #---------------------------------------------------------------------------
  if (opt$verbose) { print(opt) }
  
  return(opt)
}

#===============================================================================
# LOAD REQUIRED LIBRARIES
#===============================================================================

suppressWarnings(suppressPackageStartupMessages(library(clusterProfiler)))  # Core enrichment analysis
suppressWarnings(suppressPackageStartupMessages(library(org.Hs.eg.db)))     # Human genome annotation
suppressWarnings(suppressPackageStartupMessages(library(enrichplot)))       # Visualization of enrichment results
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))          # Grammar of graphics
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))            # Data manipulation

# Increase timeout for downloading KEGG data
options(timeout = 1000)  # Default is 60 seconds

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
  pro_diff_filter_list <- unique(pro_diff[
    ((pro_diff$logFC >= log2FC | pro_diff$logFC <= -log2FC) & pro_diff$adj.P.Val < p_adj),
  ][[omics1_name]])
  
  # Filter differentially expressed phosphosites
  phos_diff_filter_list <- phos_diff[
    ((phos_diff$logFC >= log2FC | phos_diff$logFC <= -log2FC) & phos_diff$adj.P.Val < p_adj),
  ][[omics2_name]]
  
  # Map phosphosites back to their corresponding proteins
  phos_diff_filter_list <- unique(phos_pro[
    phos_pro[[omics2_name]] %in% phos_diff_filter_list,
  ][[omics1_name]])
  
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
                                      color_gradient_high = "#90362d") {
  
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
  
  #---------------------------------------------------------------------------
  # Create GO enrichment dot plot
  #---------------------------------------------------------------------------
  go_compare_GO <- ggplot(go_compare_clean, aes(x = Cluster, y = Description)) +
    
    # Add points with size = gene count, color = -log10(adjusted p-value)
    geom_point(aes(size = Count, color = log_p_adjust), alpha = 0.8) +
    
    # Facet by GO ontology with free y-axis scales
    facet_grid(ONTOLOGY ~ .,
               scales = "free_y",
               space = "free_y") +
    
    # Color gradient for -log10(p-values)
    scale_color_gradient(low = color_gradient_low, high = color_gradient_high,
                         name = "-log10(Adjusted P-value)") +
    
    # Size scale for point sizes
    scale_size_continuous(range = c(2, 8),
                          name = paste0(omics1_name, "/", omics2_name, " Count"),
                          guide = guide_legend(
                            override.aes = list(color = "grey60")
                          )) +
    
    # Axis labels and title
    labs(x = "Protein Ratio",
         y = "",
         title = paste0("GO Enrichment Comparison: ", omics1_name, " vs ", omics2_name)) +
    
    # Theme customization
    theme_bw() +
    theme(
      # Title settings
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16,
                                margin = margin(b = 20)),
      
      # Facet strip settings
      strip.text = element_text(size = 12, face = "bold",
                                margin = margin(t = 5, r = 5, b = 5, l = 5)),
      strip.background = element_rect(fill = "lightgray"),
      
      # Axis settings
      axis.text.x = element_text(size = 10, color = "grey20", face = "bold"),
      axis.text.y = element_text(size = 10, color = "grey20", face = "bold"),
      axis.title.x = element_text(size = 10, face = "bold",
                                  margin = margin(t = 15)),
      
      # Legend settings
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
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
  ggsave(paste0(outdir, "/", omics1_name, "_", omics2_name, "_GO_enrichment.png"),
         go_compare_GO, width = 12, height = 8, dpi = 300)
}

#===============================================================================
# KEGG ENRICHMENT COMPARISON FUNCTION
#===============================================================================

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
                                        color_gradient_high = "#90362d") {
  
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
    
    # Theme customization
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16,
                                margin = margin(b = 20)),
      axis.text.x = element_text(size = 10, color = "grey20", face = "bold"),
      axis.text.y = element_text(size = 10, color = "grey20", face = "bold"),
      axis.title.x = element_text(size = 10, face = "bold",
                                  margin = margin(t = 15)),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.8, "cm"),
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
#'   - Prepares data from differential expression results
#'   - Performs GO enrichment comparison
#'   - Performs KEGG enrichment comparison
#'
#' @param pro_diff_path Path to proteomics differential expression results
#' @param phos_diff_path Path to phosphoproteomics differential expression results
#' @param phos_pro_path Path to protein-phosphosite mapping file
#' @param outdir Output directory
#' @param log2FC Log2 fold change threshold
#' @param diff_p_adj Adjusted p-value threshold
#' @param omics1_name Name for first omics dataset
#' @param omics2_name Name for second omics dataset
#' @param pvalueCutoff P-value cutoff for enrichment
#' @param GO_showCategory Number of top GO terms to show
#' @param KEGG_showCategory Number of top KEGG pathways to show
#' @param enrich_fromType ID type for enrichment
#' @param color_gradient_low Low end color
#' @param color_gradient_high High end color
omics_enrichment <- function(pro_diff_path, phos_diff_path, phos_pro_path,
                             outdir = "./", log2FC = 1.2, diff_p_adj = 0.05,
                             omics1_name = "Proteomics", omics2_name = "Phosphoproteomics",
                             pvalueCutoff = 0.05, GO_showCategory = 6, KEGG_showCategory = 15,
                             enrich_fromType = 'UNIPROT',
                             color_gradient_low = "#175663",
                             color_gradient_high = "#90362d") {
  
  # Prepare gene lists from differential expression results
  predata <- enrichment_predata(pro_diff_path, phos_diff_path, phos_pro_path,
                                outdir, log2FC, diff_p_adj)
  omics1_list <- predata$protein_list
  omics2_list <- predata$phosphoprotein_list
  
  #---------------------------------------------------------------------------
  # GO Enrichment Comparison
  #---------------------------------------------------------------------------
  two_omicses_GO_enrichment(omics1_list, omics2_list,
                            omics1_name, omics2_name,
                            outdir, pvalueCutoff, GO_showCategory,
                            enrich_fromType = enrich_fromType,
                            color_gradient_low = color_gradient_low,
                            color_gradient_high = color_gradient_high)
  
  #---------------------------------------------------------------------------
  # KEGG Enrichment Comparison
  #---------------------------------------------------------------------------
  two_omicses_KEGG_enrichment(omics1_list, omics2_list,
                              omics1_name, omics2_name,
                              outdir, pvalueCutoff, KEGG_showCategory,
                              enrich_fromType = enrich_fromType,
                              color_gradient_low = color_gradient_low,
                              color_gradient_high = color_gradient_high)
}

#===============================================================================
# PARSE PARAMETERS AND EXECUTE MAIN ANALYSIS
#===============================================================================

# Validate command line arguments
opt <- parameter_validation(opt)

# Extract parameters into variables
pro_diff_path <- opt$pro_diff
phos_diff_path <- opt$phos_diff
phos_pro_path <- opt$phos_pro
outdir <- opt$outdir
log2FC <- opt$log2FC
diff_p_adj <- opt$diff_p_adj
omics1_name <- opt$omics1_name
omics2_name <- opt$omics2_name
GO_showCategory <- opt$GO_showCategory
KEGG_showCategory <- opt$KEGG_showCategory
pvalueCutoff <- opt$pvalueCutoff
color_gradient_low <- opt$color_gradient_low
color_gradient_high <- opt$color_gradient_high
enrich_fromType <- opt$enrich_fromType

# Execute main enrichment analysis
omics_enrichment(pro_diff_path, phos_diff_path, phos_pro_path,
                 outdir, log2FC, diff_p_adj,
                 omics1_name, omics2_name,
                 pvalueCutoff, GO_showCategory, KEGG_showCategory,
                 enrich_fromType = enrich_fromType,
                 color_gradient_low = color_gradient_low,
                 color_gradient_high = color_gradient_high)