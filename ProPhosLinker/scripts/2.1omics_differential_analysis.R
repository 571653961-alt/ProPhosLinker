#!/usr/bin/env Rscript

#===============================================================================
# Multi-omics Differential Expression Analysis
#===============================================================================
# This script performs differential expression analysis for paired proteomics
# and phosphoproteomics data, integrates the results, and visualizes the
# relationship between protein and phosphosite changes using quadrant plots.
#
# Key features:
#   - Differential expression analysis using limma
#   - Integration of protein and phosphosite-level results
#   - Quadrant plots showing coordinated regulation patterns
#   - Support for pre-computed differential expression results
#===============================================================================

#===============================================================================
# COMMAND LINE ARGUMENT DEFINITION
#===============================================================================

suppressWarnings(library(optparse))

option_list <- list(
  #---------------------------------------------------------------------------
  # Input Parameters
  #---------------------------------------------------------------------------
  make_option(c("--profile"), 
              type = "character", 
              default = NULL, 
              metavar = "character", 
              help = "Path to proteomics abundance data file (TSV format)"),
  
  make_option(c("--phosfile"), 
              type = "character", 
              default = NULL, 
              metavar = "character", 
              help = "Path to phosphoproteomics abundance data file (TSV format)"),
  
  make_option(c("--pro_phos_cor"), 
              type = "character", 
              default = NULL, 
              metavar = "character", 
              help = "File path mapping phosphoprotein site IDs to protein IDs"),
  
  make_option(c("--compare_groups"), 
              type = "character", 
              default = NULL, 
              metavar = "character", 
              help = "Path to comparison groups file (TSV format with sample and group columns)"),
  
  make_option(c("--prodiff"), 
              type = "character", 
              default = NULL, 
              metavar = "PATH", 
              help = "Path to pre-computed proteomics differential expression results (optional). 
                      Must contain columns: protein ID, log2 fold change, and adjusted p-value."),
  
  make_option(c("--phosdiff"), 
              type = "character", 
              default = NULL, 
              metavar = "PATH", 
              help = "Path to pre-computed phosphoproteomics differential expression results (optional).
                      Must contain columns: phosphosite ID, log2 fold change, and adjusted p-value."),
  
  #---------------------------------------------------------------------------
  # Output Parameters
  #---------------------------------------------------------------------------
  make_option(c("-o", "--outdir"), 
              type = "character", 
              default = NULL, 
              metavar = "character",
              help = "Output directory for results [Default: current directory]"),
  
  #---------------------------------------------------------------------------
  # Analysis Parameters
  #---------------------------------------------------------------------------
  make_option(c("--log2FC_pro"), 
              type = "numeric", 
              default = 1.2, 
              metavar = "numeric",
              help = "Log2 fold change threshold for proteomics significance [Default: %default]"),
  
  make_option(c("--log2FC_phos"), 
              type = "numeric", 
              default = 1.2, 
              metavar = "numeric",
              help = "Log2 fold change threshold for phosphoproteomics significance [Default: %default]"),
  
  make_option(c("--p_val"), 
              type = "numeric", 
              default = 0.05, 
              metavar = "numeric",
              help = "Adjusted p-value threshold for statistical significance [Default: %default]"),
  
  #---------------------------------------------------------------------------
  # General Parameters
  #---------------------------------------------------------------------------
  make_option(c("--omics1_name"), 
              type = "character", 
              default = "Pro", 
              metavar = "character",
              help = "Name for first omics dataset (proteomics) [Default: %default]"),
  
  make_option(c("--omics2_name"), 
              type = "character", 
              default = "Phos", 
              metavar = "character",
              help = "Name for second omics dataset (phosphoproteomics) [Default: %default]"),
  
  make_option(c("--group_comparing"), 
              type = "character", 
              default = "T:N", 
              metavar = "character",
              help = "Comparison groups in format 'case:control' (e.g., 'T:N') [Default: %default]"),
  
  #---------------------------------------------------------------------------
  # Visualization Parameters
  #---------------------------------------------------------------------------
  make_option(c("--quadrant_plot_up_color"), 
              type = "character", 
              default = "#a03c32", 
              metavar = "character",
              help = "Color for upregulated points in quadrant plot [Default: %default]"),
  
  make_option(c("--quadrant_plot_down_color"), 
              type = "character", 
              default = "#43656C", 
              metavar = "character",
              help = "Color for downregulated points in quadrant plot [Default: %default]"),
  
  #---------------------------------------------------------------------------
  # Verbose Mode
  #---------------------------------------------------------------------------
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = FALSE, 
              metavar = "FLAG",
              help = "Print detailed output messages [Default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(
  option_list = option_list,
  description = "Differential Expression Analysis for Proteomics and Phosphoproteomics Data"
)
opt <- parse_args(opt_parser)

#===============================================================================
# PARAMETER VALIDATION FUNCTION
#===============================================================================

#' Validate and clean input parameters
#'
#' This function performs comprehensive validation of all command-line parameters:
#'   - Checks for required input files
#'   - Cleans and normalizes file paths
#'   - Verifies file existence and non-emptiness
#'   - Creates output directory if needed
#'   - Validates numeric parameters within ranges
#'
#' @param opt List of parsed command-line arguments
#' @return Validated and cleaned parameter list
parameter_validation <- function(opt) {
  
  #---------------------------------------------------------------------------
  # Helper Functions
  #---------------------------------------------------------------------------
  
  #' Check if required file parameter is provided
  #' @param file_path Path to file (may be NULL)
  #' @param file_name Name of file for error message
  check_required_file <- function(file_path, file_name) {
    if (is.null(file_path)) {
      stop(paste0("Error: ", file_name, " file is required."))
    }
  }
  
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
    if (is.null(result_dir) || (length(result_dir) == 0)) {
      result_dir <- getwd()
    }
    if (!dir.exists(result_dir)) {
      # Recursively create directory and all parent directories
      dir.create(result_dir, recursive = TRUE)
      print(paste("✅ Output directory created:", result_dir))
    } else {
      print(paste("✅ Output directory already exists:", result_dir))
    }
    return(result_dir)
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
      stop("❌ Error: ", option_name, " must be within the range [", 
           min_val, ", ", max_val, "].")
    }
  }
  
  #---------------------------------------------------------------------------
  # Required File Validation
  #---------------------------------------------------------------------------
  check_required_file(opt$profile, opt$omics1_name)
  check_required_file(opt$phosfile, opt$omics2_name)
  check_required_file(opt$compare_groups, "compare groups")
  check_required_file(opt$pro_phos_cor, "phos site IDs to protein IDs mapping")
  
  #---------------------------------------------------------------------------
  # Path Cleaning
  #---------------------------------------------------------------------------
  opt$profile <- clean_path(opt$profile)
  opt$phosfile <- clean_path(opt$phosfile)
  opt$pro_phos_cor <- clean_path(opt$pro_phos_cor)
  opt$compare_groups <- clean_path(opt$compare_groups)
  opt$outdir <- clean_path(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Input File Existence Check
  #---------------------------------------------------------------------------
  if (!is.null(opt$profile)) { check_input_file(opt$profile, opt$omics1_name) }
  if (!is.null(opt$phosfile)) { check_input_file(opt$phosfile, opt$omics2_name) }
  if (!is.null(opt$pro_phos_cor)) { check_input_file(opt$pro_phos_cor, "protein-phosphosite mapping") }
  if (!is.null(opt$compare_groups)) { check_input_file(opt$compare_groups, "compare groups") }
  if (!is.null(opt$prodiff)) { check_input_file(opt$prodiff, "pre-computed proteomics DE results") }
  if (!is.null(opt$phosdiff)) { check_input_file(opt$phosdiff, "pre-computed phosphoproteomics DE results") }
  
  #---------------------------------------------------------------------------
  # Output Directory Setup
  #---------------------------------------------------------------------------
  opt$outdir <- check_output_dir(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Numeric Parameter Validation
  #---------------------------------------------------------------------------
  validate_numeric_range(opt$log2FC_pro, 0, 10, "log2FC_pro")
  validate_numeric_range(opt$log2FC_phos, 0, 10, "log2FC_phos")
  validate_numeric_range(opt$p_val, 0, 1, "p_val")
  
  #---------------------------------------------------------------------------
  # Verbose Output
  #---------------------------------------------------------------------------
  if (opt$verbose) { print(opt) }
  
  return(opt)
}

#===============================================================================
# LOAD REQUIRED LIBRARIES
#===============================================================================

suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(limma)))
suppressWarnings(suppressPackageStartupMessages(library(statmod)))
suppressWarnings(suppressPackageStartupMessages(library(ggrepel)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(colorspace)))

#===============================================================================
# PARSE AND VALIDATE PARAMETERS
#===============================================================================

opt <- parameter_validation(opt)

# Extract parameters into variables for easier access
pro_path <- opt$profile
phos_path <- opt$phosfile
prodiff_path <- opt$prodiff
phosdiff_path <- opt$phosdiff
pro_phos_path <- opt$pro_phos_cor
group_path <- opt$compare_groups
log2FC_Pro_cut <- opt$log2FC_pro
log2FC_Phos_cut <- opt$log2FC_phos
PVal_cut <- opt$p_val
outdir <- opt$outdir
omics1_name <- opt$omics1_name
omics2_name <- opt$omics2_name
group_comparing <- strsplit(opt$group_comparing, ':')[[1]]
quadrant_plot_up_color <- opt$quadrant_plot_up_color
quadrant_plot_down_color <- opt$quadrant_plot_down_color

#===============================================================================
# COLOR UTILITY FUNCTIONS
#===============================================================================

#' Convert various color formats to hexadecimal
#'
#' Supports hex colors, color names, and RGB vectors
#'
#' @param color Color specification (hex string, color name, or RGB vector)
#' @return Hexadecimal color string
convert_to_hex <- function(color) {
  if (is.character(color)) {
    if (grepl("^#", color)) {
      # Already in hexadecimal format
      return(color)
    } else {
      # Color name, convert to hex
      tryCatch({
        rgb_val <- col2rgb(color)
        return(rgb(rgb_val[1, ], rgb_val[2, ], rgb_val[3, ], maxColorValue = 255))
      }, error = function(e) {
        stop("Unrecognized color name: ", color)
      })
    }
  } else if (is.numeric(color) && length(color) == 3) {
    # RGB vector, e.g., c(255, 0, 0)
    return(rgb(color[1], color[2], color[3], maxColorValue = 255))
  } else {
    stop("Unsupported color format")
  }
}

#' Adjust color brightness
#'
#' Lightens or darkens a color by a specified percentage
#'
#' @param color Input color (hex, name, or RGB)
#' @param percent Percentage to adjust brightness (positive = lighter, negative = darker)
#' @return Adjusted color in hexadecimal format
adjust_color_brightness <- function(color, percent) {
  # Convert to hexadecimal
  hex_color <- convert_to_hex(color)
  
  # Convert to HSL color space (Hue, Saturation, Lightness)
  hsl_color <- hex2RGB(hex_color)
  hsl_color <- as(hsl_color, "HLS")
  
  # Get current HSL values
  hsl_matrix <- coords(hsl_color)
  
  if (percent > 0) {
    # Lighten: increase lightness
    hsl_matrix[, "L"] <- hsl_matrix[, "L"] + (1 - hsl_matrix[, "L"]) * (percent / 100)
  } else {
    # Darken: decrease lightness
    hsl_matrix[, "L"] <- hsl_matrix[, "L"] * (1 - abs(percent) / 100)
  }
  
  # Ensure lightness stays within valid range [0, 1]
  hsl_matrix[, "L"] <- pmin(pmax(hsl_matrix[, "L"], 0), 1)
  
  # Convert back to hexadecimal
  adjusted_color <- hex(HLS(hsl_matrix))
  
  return(adjusted_color)
}

#===============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS FUNCTION (limma)
#===============================================================================

#' Perform differential expression analysis using limma
#'
#' This function conducts differential expression analysis on omics data
#' using the limma pipeline with empirical Bayes moderation.
#'
#' @param omics_path Path to omics abundance data file
#' @param group_path Path to group information file
#' @param outdir Output directory for results
#' @param PVal Adjusted p-value threshold
#' @param log2FC Log2 fold change threshold
#' @param group_comparing Vector of two group names for comparison
#' @param omics_type Name of omics type for labeling
#' @return Data frame with differential expression results
DEP_analysis <- function(omics_path, group_path, outdir, PVal, log2FC, 
                         group_comparing, omics_type) {
  
  # Load and log-transform omics data
  omics_data <- read.csv(omics_path, sep = '\t', header = TRUE, row.names = 1)
  omics_data <- log2(omics_data + 1)  # Add pseudocount to avoid log(0)
  
  # Load group information
  groups_data <- read.csv(group_path, sep = '\t', header = TRUE)
  
  # Ensure sample order matches between data and group information
  omics_data <- omics_data[, groups_data[, 1]]
  
  # Set factor levels according to comparison order (case vs control)
  groups <- factor(groups_data$group, 
                   levels = unique(c(group_comparing[1], group_comparing[2])))
  
  #---------------------------------------------------------------------------
  # Design Matrix Construction
  #---------------------------------------------------------------------------
  # Create design matrix for linear model (no intercept, one column per group)
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)
  rownames(design) <- colnames(omics_data)
  
  #---------------------------------------------------------------------------
  # Contrast Matrix Construction
  #---------------------------------------------------------------------------
  # Create contrast for specified comparison (e.g., "T - N")
  contrast_matrix <- makeContrasts(
    T_vs_N = paste0(group_comparing[1], " - ", group_comparing[2]),
    levels = design
  )
  
  #---------------------------------------------------------------------------
  # limma Pipeline
  #---------------------------------------------------------------------------
  # Fit linear model for each feature
  fit <- lmFit(omics_data, design)  # limma expects samples in columns
  
  # Apply contrast to compare specified groups
  fit2 <- contrasts.fit(fit, contrast_matrix)
  
  # Empirical Bayes moderation to improve variance estimates
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  
  # Extract all results
  tempOutput <- topTable(fit2, n = Inf, adjust = "fdr")
  nrDEG <- na.omit(tempOutput)
  
  # Add feature identifiers as a column
  nrDEG_with_omics <- nrDEG %>%
    rownames_to_column(omics_type)
  
  # Classify features as UP, DOWN, or Non-significant based on thresholds
  nrDEG_with_omics$class <- ifelse(
    nrDEG_with_omics$logFC >= log2FC & nrDEG_with_omics$adj.P.Val < PVal,
    "UP",
    ifelse(
      nrDEG_with_omics$logFC <= -log2FC & nrDEG_with_omics$adj.P.Val < PVal,
      "DOWN",
      "Non-significant"
    )
  )
  
  # Save results to file
  write.table(
    nrDEG_with_omics, 
    file.path(outdir, paste0(group_comparing[1], "_vs_", group_comparing[2], "_", omics_type, ".tsv")),
    sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  return(nrDEG_with_omics)
}

#===============================================================================
# QUADRANT PLOT FUNCTION
#===============================================================================

#' Create quadrant plot showing relationship between protein and phosphosite changes
#'
#' This function generates a 9-quadrant plot that visualizes the coordinated
#' regulation patterns between proteins and their corresponding phosphosites.
#'
#' @param pro_phos_path Path to protein-phosphosite mapping file
#' @param pro_DEP Proteomics differential expression results
#' @param phos_DEP Phosphoproteomics differential expression results
#' @param outdir Output directory
#' @param log2FC_Pro_cut Log2 FC threshold for proteomics
#' @param log2FC_Phos_cut Log2 FC threshold for phosphoproteomics
#' @param adj_PVal_cut Adjusted p-value threshold
#' @param omics1_type Name of proteomics dataset
#' @param omics2_type Name of phosphoproteomics dataset
#' @param quadrant_plot_up_color Color for up-up quadrant
#' @param quadrant_plot_down_color Color for down-down quadrant
#' @param add_fit_line Whether to add a trend line
#' @param fit_method Method for trend line ("lm", "loess", "gam")
#' @param show_correlation Whether to display correlation coefficient
#' @param line_color Color of trend line
#' @param line_size Size of trend line
#' @param line_type Type of trend line (solid, dashed, etc.)
quadrant_plot <- function(pro_phos_path, pro_DEP, phos_DEP, outdir,
                          log2FC_Pro_cut, log2FC_Phos_cut, adj_PVal_cut,
                          omics1_type, omics2_type,
                          quadrant_plot_up_color, quadrant_plot_down_color,
                          add_fit_line = TRUE,  
                          fit_method = "lm",
                          show_correlation = TRUE,
                          line_color = "black",    
                          line_size = 0.8,        
                          line_type = "dashed") {
  
  #---------------------------------------------------------------------------
  # Data Integration
  #---------------------------------------------------------------------------
  # Load protein-phosphosite mapping file
  pro_phos <- read.csv(pro_phos_path, sep = '\t', header = TRUE)
  colnames(pro_phos) <- c(omics1_type, omics2_type)
  
  # Merge protein mapping with proteomics DE results
  combined_df <- merge(pro_phos, pro_DEP[, c(omics1_type, "logFC", "adj.P.Val")], 
                       by = omics1_type, all.x = TRUE)
  combined_df <- na.omit(combined_df)
  colnames(combined_df)[3] <- paste0(omics1_type, "_log2FC")
  
  # Merge with phosphoproteomics DE results
  combined_df <- merge(combined_df, phos_DEP[, c(omics2_type, "logFC", "adj.P.Val")],
                       by.x = omics2_type, by.y = omics2_type, all.x = TRUE)
  colnames(combined_df)[5] <- paste0(omics2_type, "_log2FC")
  
  combined_df <- na.omit(combined_df)
  
  # Filter by adjusted p-value thresholds
  combined_df <- combined_df[combined_df$adj.P.Val.x < adj_PVal_cut &
                               combined_df$adj.P.Val.y < adj_PVal_cut, ]
  combined_df <- combined_df %>%
    dplyr::select(-adj.P.Val.x, -adj.P.Val.y)
  
  #---------------------------------------------------------------------------
  # Classification into 9 quadrants
  #---------------------------------------------------------------------------
  combined_df <- combined_df %>%
    mutate(class = case_when(
      # UP in both
      get(paste0(omics1_type, "_log2FC")) > log2FC_Pro_cut & 
        get(paste0(omics2_type, "_log2FC")) > log2FC_Phos_cut ~ 'UP|UP',
      
      # UP in proteomics only
      get(paste0(omics1_type, "_log2FC")) > log2FC_Pro_cut & 
        abs(get(paste0(omics2_type, "_log2FC"))) <= log2FC_Phos_cut ~ 'UP|Non-significant',
      
      # UP in proteomics, DOWN in phosphoproteomics
      get(paste0(omics1_type, "_log2FC")) > log2FC_Pro_cut & 
        get(paste0(omics2_type, "_log2FC")) < -log2FC_Phos_cut ~ 'UP|DOWN',
      
      # DOWN in proteomics, UP in phosphoproteomics
      get(paste0(omics1_type, "_log2FC")) < -log2FC_Pro_cut & 
        get(paste0(omics2_type, "_log2FC")) > log2FC_Phos_cut ~ 'DOWN|UP',
      
      # DOWN in both
      get(paste0(omics1_type, "_log2FC")) < -log2FC_Pro_cut & 
        get(paste0(omics2_type, "_log2FC")) < -log2FC_Phos_cut ~ 'DOWN|DOWN',
      
      # DOWN in proteomics only
      get(paste0(omics1_type, "_log2FC")) < -log2FC_Pro_cut & 
        abs(get(paste0(omics2_type, "_log2FC"))) <= log2FC_Phos_cut ~ 'DOWN|Non-significant',
      
      # UP in phosphoproteomics only
      abs(get(paste0(omics1_type, "_log2FC"))) <= log2FC_Pro_cut & 
        get(paste0(omics2_type, "_log2FC")) > log2FC_Phos_cut ~ 'Non-significant|UP',
      
      # DOWN in phosphoproteomics only
      abs(get(paste0(omics1_type, "_log2FC"))) <= log2FC_Pro_cut & 
        get(paste0(omics2_type, "_log2FC")) < -log2FC_Phos_cut ~ 'Non-significant|DOWN',
      
      # Non-significant in both
      TRUE ~ 'Non-significant|Non-significant'
    ))
  
  # Save integrated results
  write.table(
    combined_df, 
    file.path(outdir, paste0(group_comparing[1], "_vs_", group_comparing[2], "_combined.tsv")),
    sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  #---------------------------------------------------------------------------
  # Correlation Analysis
  #---------------------------------------------------------------------------
  # Calculate Spearman correlation between protein and phosphosite log2FC
  cor_test <- cor.test(combined_df[[paste0(omics1_type, "_log2FC")]], 
                       combined_df[[paste0(omics2_type, "_log2FC")]], 
                       method = "spearman",
                       exact = FALSE)
  
  cor_value <- round(cor_test$estimate, 3)
  
  # Save correlation results
  cor_results <- c(cor_value = cor_value)
  write_tsv(as.data.frame(t(cor_results)), 
            file.path(outdir, "cor_value_results.tsv"), 
            col_names = TRUE)
  
  #---------------------------------------------------------------------------
  # Quadrant Plot Construction
  #---------------------------------------------------------------------------
  p <- ggplot(combined_df, aes(x = .data[[paste0(omics1_type, "_log2FC")]], 
                               y = .data[[paste0(omics2_type, "_log2FC")]])) +
    geom_point(aes(color = class, alpha = class), size = 2) +
    
    # Add trend line if requested
    {
      if (add_fit_line) {
        if (fit_method == "lm") {
          geom_smooth(method = "lm", formula = y ~ x, 
                      color = "#664B3A", linewidth = line_size, linetype = line_type,
                      se = TRUE, alpha = 0.2, fill = "#E0552C")
        } else if (fit_method == "loess") {
          geom_smooth(method = "loess", formula = y ~ x,
                      color = "#664B3A", linewidth = line_size, linetype = line_type,
                      se = TRUE, alpha = 0.2, fill = "#E0552C")
        } else if (fit_method == "gam") {
          geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                      color = "#664B3A", linewidth = line_size, linetype = line_type,
                      se = TRUE, alpha = 0.2, fill = "#E0552C")
        }
      }
    } +
    
    # Add threshold lines
    geom_vline(xintercept = c(-log2FC_Pro_cut, log2FC_Pro_cut), 
               color = "darkgrey", linetype = "dashed", linewidth = 0.7) +
    geom_hline(yintercept = c(-log2FC_Phos_cut, log2FC_Phos_cut), 
               color = "darkgrey", linetype = "dashed", linewidth = 0.7) +
    
    # Color mapping with adjusted brightness for different categories
    scale_color_manual(values = c(
      "UP|UP" = quadrant_plot_up_color, 
      "DOWN|DOWN" = quadrant_plot_down_color,
      'UP|Non-significant' = adjust_color_brightness(quadrant_plot_up_color, 40), 
      'Non-significant|UP' = adjust_color_brightness(quadrant_plot_up_color, 40),
      'DOWN|Non-significant' = adjust_color_brightness(quadrant_plot_down_color, 40), 
      'Non-significant|DOWN' = adjust_color_brightness(quadrant_plot_down_color, 40),
      "UP|DOWN" = adjust_color_brightness(quadrant_plot_up_color, -50), 
      "DOWN|UP" = adjust_color_brightness(quadrant_plot_down_color, -50),
      "Non-significant|Non-significant" = "grey"
    )) +
    
    # Alpha transparency mapping
    scale_alpha_manual(values = c(
      "Non-significant|Non-significant" = 0.4,
      'UP|Non-significant' = 0.6, 'DOWN|Non-significant' = 0.6, 
      'Non-significant|UP' = 0.6, 'Non-significant|DOWN' = 0.6,
      "UP|UP" = 1, "DOWN|DOWN" = 0.9, "DOWN|UP" = 0.9, "UP|DOWN" = 0.9
    ), guide = "none") +
    
    # Add threshold labels
    annotate("text",
             x = -1 * log2FC_Pro_cut - 0.3,
             y = max(combined_df[[paste0(omics2_type, "_log2FC")]], na.rm = TRUE) + 0.5,
             label = paste0(omics1_type, " log2FC = -", log2FC_Pro_cut, " "),
             color = "grey40",
             size = 3,
             fontface = 'bold') +
    
    annotate("text",
             x = log2FC_Pro_cut + 0.3,
             y = max(combined_df[[paste0(omics2_type, "_log2FC")]], na.rm = TRUE) + 0.5,
             label = paste0(omics1_type, " log2FC = +", log2FC_Pro_cut, " "),
             color = "grey40",
             size = 3,
             fontface = 'bold') +
    
    annotate("text",
             y = -1 * log2FC_Phos_cut - 0.2,
             x = max(combined_df[[paste0(omics1_type, "_log2FC")]], na.rm = TRUE) + 0.1,
             label = paste0(omics2_type, " log2FC\n = -", log2FC_Phos_cut),
             color = "grey40",
             size = 3,
             fontface = 'bold') +
    
    annotate("text",
             y = log2FC_Phos_cut + 0.2,
             x = max(combined_df[[paste0(omics1_type, "_log2FC")]], na.rm = TRUE) + 0.1,
             label = paste0(omics2_type, " log2FC\n = +", log2FC_Phos_cut),
             color = "grey40",
             size = 3,
             fontface = 'bold') +
    
    # Labels and theme
    labs(title = paste0(omics1_type, " vs. ", omics2_type, " Log2FC Quadrant Plot", 
                        " (", group_comparing[1], " vs. ", group_comparing[2], ")"),
         x = paste0("Log2FC of ", omics1_type),
         y = paste0("Log2FC of ", omics2_type),
         color = "Expression (Pro|Phos)") +
    theme_classic() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  
  # Add correlation annotation if requested
  if (show_correlation) {
    p <- p + annotate("text",
                      x = min(combined_df[[paste0(omics1_type, "_log2FC")]], na.rm = TRUE),
                      y = max(combined_df[[paste0(omics2_type, "_log2FC")]], na.rm = TRUE) - 0.5,
                      label = paste0("r = ", cor_value, "\n"),
                      color = "#E0552C",
                      size = 6,
                      fontface = "bold",
                      hjust = 0,
                      vjust = 1)
  }
  
  # Save plot
  ggsave(file.path(outdir, "omics_diff_quadrant_plot.png"), 
         width = 15, height = 8, dpi = 300)
}

#===============================================================================
# MAIN ANALYSIS FUNCTION
#===============================================================================

#' Main function orchestrating differential expression and quadrant plot analysis
#'
#' @param pro_path Path to proteomics data
#' @param phos_path Path to phosphoproteomics data
#' @param group_path Path to group information
#' @param outdir Output directory
#' @param PVal_cut Adjusted p-value threshold
#' @param log2FC_Pro_cut Log2 FC threshold for proteomics
#' @param log2FC_Phos_cut Log2 FC threshold for phosphoproteomics
#' @param pro_phos_path Path to protein-phosphosite mapping
#' @param omics1_type Name of proteomics dataset
#' @param omics2_type Name of phosphoproteomics dataset
#' @param group_comparing Vector of two group names for comparison
#' @param quadrant_plot_up_color Color for up-up quadrant
#' @param quadrant_plot_down_color Color for down-down quadrant
DE_pro_phos <- function(pro_path, phos_path, group_path, outdir, PVal_cut,
                        log2FC_Pro_cut, log2FC_Phos_cut, pro_phos_path,
                        omics1_type = "Pro",
                        omics2_type = "Phos",
                        group_comparing,
                        quadrant_plot_up_color = "#a03c32",
                        quadrant_plot_down_color = "#1a5f6e") {
  
  #---------------------------------------------------------------------------
  # Proteomics Differential Expression
  #---------------------------------------------------------------------------
  if (is.null(prodiff_path)) {
    # Run limma analysis if no pre-computed results provided
    pro_DEP <- DEP_analysis(pro_path, group_path, outdir, PVal_cut, 
                            log2FC_Pro_cut, group_comparing, omics1_type)
  } else {
    # Load pre-computed results
    pro_DEP <- read.csv(prodiff_path, sep = '\t', header = TRUE)
    write.table(pro_DEP, 
                file.path(outdir, paste0(group_comparing[1], "_vs_", group_comparing[2], "_", omics1_type, ".tsv")),
                sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
  #---------------------------------------------------------------------------
  # Phosphoproteomics Differential Expression
  #---------------------------------------------------------------------------
  if (is.null(phosdiff_path)) {
    # Run limma analysis if no pre-computed results provided
    phos_DEP <- DEP_analysis(phos_path, group_path, outdir, PVal_cut, 
                             log2FC_Phos_cut, group_comparing, omics2_type)
  } else {
    # Load pre-computed results
    phos_DEP <- read.csv(phosdiff_path, sep = '\t', header = TRUE)
    write.table(phos_DEP, 
                file.path(outdir, paste0(group_comparing[1], "_vs_", group_comparing[2], "_", omics2_type, ".tsv")),
                sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
  #---------------------------------------------------------------------------
  # Quadrant Plot Generation
  #---------------------------------------------------------------------------
  quadrant_plot(pro_phos_path, pro_DEP, phos_DEP, outdir,
                log2FC_Pro_cut, log2FC_Phos_cut, adj_PVal_cut = 0.05,
                omics1_type, omics2_type,
                quadrant_plot_up_color, quadrant_plot_down_color,
                add_fit_line = TRUE,
                fit_method = "loess",
                show_correlation = TRUE,
                line_color = "#F6A01D",
                line_size = 1.6,
                line_type = "solid")
}

#===============================================================================
# EXECUTE MAIN ANALYSIS
#===============================================================================

DE_pro_phos(
  pro_path, phos_path, group_path, outdir, PVal_cut,
  log2FC_Pro_cut, log2FC_Phos_cut, pro_phos_path,
  omics1_name, omics2_name, group_comparing,
  quadrant_plot_up_color, quadrant_plot_down_color
)