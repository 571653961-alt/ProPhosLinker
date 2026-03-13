#!/usr/bin/env Rscript

#===============================================================================
# Phosphorylation Rate Quantile Subtyping using Mfuzz Clustering
#===============================================================================
# This script performs clustering analysis on phosphorylation rate data across
# protein intensity quantiles using the Mfuzz soft clustering algorithm.
# It identifies patterns of phosphorylation rate changes and visualizes them
# in a 3x3 facet plot with color-coded trends.
#
# Key features:
#   - Mfuzz soft clustering of phosphorylation rate patterns
#   - Classification of phosphorylation trends (up, down, both, or non-significant)
#   - Multi-panel visualization with customizable colors
#===============================================================================

suppressWarnings(library(optparse))

#===============================================================================
# COMMAND LINE ARGUMENT DEFINITION
#===============================================================================

option_list <- list(
  #---------------------------------------------------------------------------
  # Input Parameters
  #---------------------------------------------------------------------------
  make_option(c("-i", "--input"), 
              type = "character", 
              default = NULL, 
              metavar = "character", 
              help = "Path to denoised phosphorylation rate quantile data file (TSV format)"),
  
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
  make_option(c("-c", "--clusternum"), 
              type = "integer", 
              default = 9, 
              metavar = "integer",
              help = "Number of clusters for Mfuzz clustering [Default: %default]"),
  
  make_option(c("-f", "--phosratefc"), 
              type = "numeric", 
              default = 1, 
              metavar = "numeric",
              help = "Phosphorylation rate fold change threshold for significance [Default: %default]"),
  
  #---------------------------------------------------------------------------
  # Visualization Color Parameters
  #---------------------------------------------------------------------------
  make_option(c("--pluscolor"), 
              type = "character", 
              default = "#a03c32", 
              metavar = "character",
              help = "Color for positive phosphorylation rate threshold line [Default: %default]"),
  
  make_option(c("--minuscolor"), 
              type = "character", 
              default = "#1a5f6e", 
              metavar = "character",
              help = "Color for negative phosphorylation rate threshold line [Default: %default]"),
  
  make_option(c("--allcolor"), 
              type = "character", 
              default = "grey60", 
              metavar = "character",
              help = "Color for all phosphorylation rate lines (deprecated) [Default: %default]"),
  
  make_option(c("--onlyupcolor"), 
              type = "character", 
              default = "#a03c32", 
              metavar = "character",
              help = "Color for lines showing only up-regulation [Default: %default]"),
  
  make_option(c("--onlydowncolor"), 
              type = "character", 
              default = "#1a5f6e", 
              metavar = "character",
              help = "Color for lines showing only down-regulation [Default: %default]"),
  
  make_option(c("--updowncolor"), 
              type = "character", 
              default = "#4D613C", 
              metavar = "character",
              help = "Color for lines showing both up and down regulation [Default: %default]"),
  
  make_option(c("--notsigcolor"), 
              type = "character", 
              default = "grey60", 
              metavar = "character",
              help = "Color for non-significant phosphorylation rate lines [Default: %default]"),
  
  #---------------------------------------------------------------------------
  # General Parameters
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
  description = "Phosphorylation Rate Quantile Subtyping using Mfuzz Clustering"
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
  
  # Check if required input file parameter is provided
  if (is.null(opt$input)) {
    stop("Error: Input file is required. Use -i or --input to specify the input file.")
  }
  
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
  # Path Cleaning and Validation
  #---------------------------------------------------------------------------
  opt$input <- clean_path(opt$input)
  opt$outdir <- clean_path(opt$outdir)
  
  # Create output directory
  opt$outdir <- check_output_dir(opt$outdir)
  
  # Check input file
  check_input_file(opt$input, "the phosRate quantile denoised file")
  
  #---------------------------------------------------------------------------
  # Numeric Parameter Validation
  #---------------------------------------------------------------------------
  validate_numeric_range(opt$clusternum, 2, 20, "clusternum")
  validate_numeric_range(opt$phosratefc, 0, 10, "phosratefc")
  
  #---------------------------------------------------------------------------
  # Verbose Output
  #---------------------------------------------------------------------------
  if (opt$verbose) { print(opt) }
  
  return(opt)
}

# Execute parameter validation
opt <- parameter_validation(opt)

#===============================================================================
# EXTRACT PARAMETERS FOR ANALYSIS
#===============================================================================

phosRate_quantile_denoised_path <- opt$input
outdir <- opt$outdir
setwd(outdir)
cluster_num <- opt$clusternum
phos_rate_FC <- opt$phosratefc

# Color parameters (renamed for clarity)
plus_phos_rate_FC_line_color <- opt$pluscolor
minus_phos_rate_FC_line_color <- opt$minuscolor
all_phosrate_quantile_line_color <- opt$allcolor
only_up_phosrate_quantile_line_color <- opt$onlyupcolor
only_down_phosrate_quantile_line_color <- opt$onlydowncolor
up_down_phosrate_quantile_line_color <- opt$updowncolor
notsig_phosrate_quantile_line_color <- opt$notsigcolor

#===============================================================================
# LOAD REQUIRED LIBRARIES
#===============================================================================

suppressWarnings(suppressPackageStartupMessages(library(Mfuzz)))      # Soft clustering for time-series/expression data
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))    # Grammar of graphics for visualization
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))   # Data reshaping (melt/cast)
suppressWarnings(suppressPackageStartupMessages(library(tibble)))     # Modern data frames
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))      # Data manipulation

#===============================================================================
# MFUZZ CLUSTERING FUNCTION
#===============================================================================

#' Perform Mfuzz soft clustering on phosphorylation rate data
#'
#' This function applies Mfuzz soft clustering algorithm to identify patterns
#' in phosphorylation rate changes across protein intensity quantiles.
#'
#' @param phosRate_quantile_denoised_data Data frame with phosphorylation rates
#'                                        (rows = phosphosites, columns = quantiles)
#' @param outdir Output directory for saving results
#' @param cluster_num Number of clusters to generate
#' @return Data frame with cluster assignments added
mfuzz_cluster <- function(phosRate_quantile_denoised_data, outdir, cluster_num) {
  
  # Convert to matrix format (required for Mfuzz)
  phosRate_quantile_denoised <- as.matrix(phosRate_quantile_denoised_data)
  
  # Create ExpressionSet object (Mfuzz data structure)
  mfuzz_class <- new('ExpressionSet', exprs = phosRate_quantile_denoised)
  
  # Set seed for reproducibility
  set.seed(124)
  
  # Estimate fuzzification parameter (m) and perform clustering
  # mestimate() calculates the optimal fuzzification parameter
  mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
  
  # Add cluster assignments to original data
  phosRate_quantile_denoised_data$cluster <- mfuzz_cluster$cluster[
    match(rownames(phosRate_quantile_denoised_data), names(mfuzz_cluster$cluster))
  ]
  
  # Save clustered data to file
  write.table(
    phosRate_quantile_denoised_data,
    file = file.path(outdir, "fitted_phosRate_quantile_denoised_cluster.tsv"),
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
  
  return(phosRate_quantile_denoised_data)
}

#===============================================================================
# VISUALIZATION FUNCTION (Commented out alternate versions)
#===============================================================================

# Note: The following functions were earlier versions of the plotting code
# and are kept commented for reference:
# - mfuzz_phosrate_subtyping_all_plot: Simple plot with all lines in same color
# - mfuzz_phosrate_subtyping_diff_plot: Plot with color-coded trends (precursor)

#' Create multi-panel visualization of phosphorylation rate patterns
#'
#' This function generates a 3x3 facet plot showing phosphorylation rate trends
#' across protein intensity quantiles for each Mfuzz cluster. Lines are color-coded
#' based on their regulation pattern relative to the threshold.
#'
#' @param phosRate_quantile_denoised_clusters_data Data frame with cluster assignments
#' @param outdir Output directory for saving plots
#' @param cluster_num Number of clusters
#' @param phos_rate_FC Phosphorylation rate fold change threshold
#' @param plus_phos_rate_FC_line_color Color for positive threshold line
#' @param minus_phos_rate_FC_line_color Color for negative threshold line
#' @param all_phosrate_quantile_line_color Color for all lines (deprecated)
#' @param only_up_phosrate_quantile_line_color Color for lines with only up-regulation
#' @param only_down_phosrate_quantile_line_color Color for lines with only down-regulation
#' @param up_down_phosrate_quantile_line_color Color for lines with both up and down regulation
#' @param notsig_phosrate_quantile_line_color Color for non-significant lines
mfuzz_phosrate_subtyping_plot <- function(phosRate_quantile_denoised_clusters_data,
                                          outdir,
                                          cluster_num,
                                          phos_rate_FC,
                                          plus_phos_rate_FC_line_color = "#a03c32",
                                          minus_phos_rate_FC_line_color = "#1a5f6e",
                                          all_phosrate_quantile_line_color = "grey60",
                                          only_up_phosrate_quantile_line_color = "#a03c32",
                                          only_down_phosrate_quantile_line_color = "#1a5f6e",
                                          up_down_phosrate_quantile_line_color = "#4D613C",
                                          notsig_phosrate_quantile_line_color = "grey60") {
  
  #---------------------------------------------------------------------------
  # Data Preparation
  #---------------------------------------------------------------------------
  
  # Initialize empty data frame to store all plot data
  all_plot_data <- data.frame()
  
  # Loop through each cluster and extract data
  for (cluster_id in 1:cluster_num) {
    # Subset data for current cluster
    cluster_value <- phosRate_quantile_denoised_clusters_data[
      phosRate_quantile_denoised_clusters_data$cluster == cluster_id, 
    ]
    
    # Check if cluster has any data
    if (nrow(cluster_value) > 0) {
      # Remove cluster column for plotting
      cluster_value <- cluster_value[, !colnames(cluster_value) %in% "cluster", drop = FALSE]
      
      # Reshape data from wide to long format for ggplot
      plot_data <- cluster_value %>%
        as.data.frame() %>%
        rownames_to_column("phos_site") %>%
        melt(id.vars = "phos_site", variable.name = "quantile", value.name = "phos_rate")
      
      # Convert quantile names to numeric (removing 'X' prefix if present)
      plot_data$quantile <- as.numeric(gsub("X", "", plot_data$quantile))
      plot_data$cluster <- paste("Cluster", cluster_id)
      
      # Append to combined data frame
      all_plot_data <- rbind(all_plot_data, plot_data)
    }
  }
  
  #---------------------------------------------------------------------------
  # Classification of Phosphorylation Patterns
  #---------------------------------------------------------------------------
  
  # Calculate max and min phosphorylation rate for each phosphosite
  line_colors <- all_plot_data %>%
    group_by(phos_site, cluster) %>%
    summarize(
      max_phos_rate = max(phos_rate, na.rm = TRUE),
      min_phos_rate = min(phos_rate, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      # Classify based on threshold crossing
      change_type = case_when(
        max_phos_rate >= phos_rate_FC & min_phos_rate <= -phos_rate_FC ~ 
          paste0("Both contain >= ", phos_rate_FC, " & <= -", phos_rate_FC),
        max_phos_rate >= phos_rate_FC ~ 
          paste0("Only contain >= ", phos_rate_FC),
        min_phos_rate <= -phos_rate_FC ~ 
          paste0("Only contain <= -", phos_rate_FC),
        TRUE ~ "Not Significant"
      ),
      # Assign colors based on classification
      line_color = case_when(
        max_phos_rate >= phos_rate_FC & min_phos_rate <= -phos_rate_FC ~ 
          up_down_phosrate_quantile_line_color,
        max_phos_rate >= phos_rate_FC ~ 
          only_up_phosrate_quantile_line_color,
        min_phos_rate <= -phos_rate_FC ~ 
          only_down_phosrate_quantile_line_color,
        TRUE ~ notsig_phosrate_quantile_line_color
      )
    )
  
  # Merge color and classification information with main data
  all_plot_data <- all_plot_data %>%
    left_join(line_colors %>% select(phos_site, cluster, line_color, change_type), 
              by = c("phos_site", "cluster"))
  
  #---------------------------------------------------------------------------
  # Plot Construction
  #---------------------------------------------------------------------------
  
  p <- ggplot(all_plot_data, aes(x = quantile, y = phos_rate, group = phos_site)) +
    
    # Main lines with conditional aesthetics
    geom_line(
      aes(color = change_type,
          alpha = change_type,
          linewidth = change_type)
    ) +
    
    # Threshold lines
    geom_hline(aes(yintercept = phos_rate_FC, linetype = "Positive Threshold"),
               color = plus_phos_rate_FC_line_color, linewidth = 0.4) +
    geom_hline(aes(yintercept = -phos_rate_FC, linetype = "Negative Threshold"),
               color = minus_phos_rate_FC_line_color, linewidth = 0.4) +
    
    # Color scale for main lines
    scale_color_manual(
      name = "Phosphorylation Pattern",
      values = setNames(
        c(up_down_phosrate_quantile_line_color, 
          only_up_phosrate_quantile_line_color, 
          only_down_phosrate_quantile_line_color, 
          notsig_phosrate_quantile_line_color),
        c(paste0("Both contain >= ", phos_rate_FC, " & <= -", phos_rate_FC),
          paste0("Only contain >= ", phos_rate_FC),
          paste0("Only contain <= -", phos_rate_FC),
          "Not Significant")
      )
    ) +
    
    # Linetype scale for threshold lines
    scale_linetype_manual(
      name = "Phos Rate Threshold Lines",
      values = c("Positive Threshold" = "dashed", "Negative Threshold" = "dashed"),
      labels = c(paste0("Phos Rate = +", phos_rate_FC), 
                 paste0("Phos Rate = -", phos_rate_FC)),
      guide = guide_legend(override.aes = list(
        color = c(plus_phos_rate_FC_line_color, minus_phos_rate_FC_line_color)
      ))
    ) +
    
    # Alpha transparency scale
    scale_alpha_manual(
      name = "Phosphorylation Pattern",
      values = setNames(
        c(0.8, 0.8, 0.8, 0.3),
        c(paste0("Both contain >= ", phos_rate_FC, " & <= -", phos_rate_FC),
          paste0("Only contain >= ", phos_rate_FC),
          paste0("Only contain <= -", phos_rate_FC),
          "Not Significant")
      ),
      guide = "none"
    ) +
    
    # Line width scale
    scale_linewidth_manual(
      name = "Phosphorylation Pattern",
      values = setNames(
        c(0.3, 0.3, 0.3, 0.2),
        c(paste0("Both contain >= ", phos_rate_FC, " & <= -", phos_rate_FC),
          paste0("Only contain >= ", phos_rate_FC),
          paste0("Only contain <= -", phos_rate_FC),
          "Not Significant")
      ),
      guide = "none"
    ) +
    
    # Create 3x3 facet grid
    facet_wrap(~ cluster, nrow = ceiling(cluster_num / 3), ncol = 3) +
    
    # Labels and theme
    labs(x = "Quantile",
         y = "Phosphorylation Rate",
         title = "Phosphorylation Rate by Protein Intensity Quantile") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      strip.background = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  # Save plot
  ggsave(file.path(outdir, "mfuzz_phosrate_subtyping.png"), 
         plot = p, width = 10, height = 6, dpi = 300)
}

#===============================================================================
# MAIN ANALYSIS FUNCTION
#===============================================================================

#' Main function for phosphorylation rate quantile subtyping
#'
#' Orchestrates the entire workflow:
#'   1. Loads denoised phosphorylation rate data
#'   2. Performs Mfuzz clustering
#'   3. Generates visualization plots
#'
#' @param phosRate_quantile_denoised_path Path to input data file
#' @param outdir Output directory for results
#' @param cluster_num Number of clusters for Mfuzz
#' @param phos_rate_FC Phosphorylation rate fold change threshold
#' @param plus_phos_rate_FC_line_color Color for positive threshold line
#' @param minus_phos_rate_FC_line_color Color for negative threshold line
#' @param all_phosrate_quantile_line_color Color for all lines (deprecated)
#' @param only_up_phosrate_quantile_line_color Color for up-regulated only
#' @param only_down_phosrate_quantile_line_color Color for down-regulated only
#' @param up_down_phosrate_quantile_line_color Color for both up and down
#' @param notsig_phosrate_quantile_line_color Color for non-significant
phosRate_quantile_subtyping <- function(phosRate_quantile_denoised_path,
                                        outdir,
                                        cluster_num = 9,
                                        phos_rate_FC = 1,
                                        plus_phos_rate_FC_line_color = "#a03c32",
                                        minus_phos_rate_FC_line_color = "#1a5f6e",
                                        all_phosrate_quantile_line_color = "grey60",
                                        only_up_phosrate_quantile_line_color = "#a03c32",
                                        only_down_phosrate_quantile_line_color = "#1a5f6e",
                                        up_down_phosrate_quantile_line_color = "#4D613C",
                                        notsig_phosrate_quantile_line_color = "grey60") {
  
  # Load input data
  phosRate_quantile_denoised_data <- read.csv(
    phosRate_quantile_denoised_path, 
    sep = '\t', 
    header = TRUE, 
    row.names = 1
  )
  
  # Perform Mfuzz clustering
  phosRate_quantile_denoised_clusters_data <- mfuzz_cluster(
    phosRate_quantile_denoised_data, 
    outdir, 
    cluster_num
  )
  
  # Generate visualization
  mfuzz_phosrate_subtyping_plot(
    phosRate_quantile_denoised_clusters_data,
    outdir,
    cluster_num,
    phos_rate_FC,
    plus_phos_rate_FC_line_color,
    minus_phos_rate_FC_line_color,
    all_phosrate_quantile_line_color,
    only_up_phosrate_quantile_line_color,
    only_down_phosrate_quantile_line_color,
    up_down_phosrate_quantile_line_color,
    notsig_phosrate_quantile_line_color
  )
}

#===============================================================================
# EXECUTE MAIN ANALYSIS
#===============================================================================

# Execute main analysis with parsed parameters
phosRate_quantile_subtyping(
  phosRate_quantile_denoised_path,
  outdir,
  cluster_num,
  phos_rate_FC,
  plus_phos_rate_FC_line_color = plus_phos_rate_FC_line_color,
  minus_phos_rate_FC_line_color = minus_phos_rate_FC_line_color,
  all_phosrate_quantile_line_color = all_phosrate_quantile_line_color,
  only_up_phosrate_quantile_line_color = only_up_phosrate_quantile_line_color,
  only_down_phosrate_quantile_line_color = only_down_phosrate_quantile_line_color,
  up_down_phosrate_quantile_line_color = up_down_phosrate_quantile_line_color,
  notsig_phosrate_quantile_line_color = notsig_phosrate_quantile_line_color
)