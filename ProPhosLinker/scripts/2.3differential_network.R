#!/usr/bin/env Rscript

#===============================================================================
# Stable Differential Network Analysis Pipeline
#===============================================================================
# This script performs a comprehensive differential network analysis between
# two conditions (e.g., tumor vs normal) for multi-omics data (proteomics and
# phosphoproteomics). The pipeline includes:
#   - Data preprocessing and filtering
#   - Differential expression analysis integration
#   - Correlation stability analysis via bootstrap
#   - Conditional network construction
#   - Differential network identification
#   - Functional enrichment analysis
#   - Subnetwork detection and visualization
#   - Mediation analysis
#===============================================================================

options(warn = -1)
suppressWarnings(library(optparse))

#===============================================================================
# COMMAND LINE ARGUMENT DEFINITION
#===============================================================================

option_list <- list(
  #---------------------------------------------------------------------------
  # Script and Input/Output Paths
  #---------------------------------------------------------------------------
  make_option(c("--script_path"),
              type = "character",
              default = NULL,
              help = "Path to the directory containing analysis scripts",
              metavar = "FILE"),
  
  # Input Parameters
  make_option(c("--profile"),
              type = "character",
              default = NULL,
              help = "Path to normalized protein abundance data file (TSV format)",
              metavar = "FILE"),
  
  make_option(c("--phosfile"),
              type = "character",
              default = NULL,
              help = "Path to normalized phosphoprotein abundance data file (TSV format)",
              metavar = "FILE"),
  
  make_option(c("--pro_phos_cor"),
              type = "character",
              default = NULL,
              help = "Path to protein-phosphosite mapping file",
              metavar = "FILE"),
  
  make_option(c("--sample_group"),
              type = "character",
              default = NULL,
              help = "Path to sample group information file",
              metavar = "FILE"),
  
  make_option(c("--diff_pro_path"),
              type = "character",
              default = NULL,
              help = "Path to differential protein expression results file",
              metavar = "FILE"),
  
  make_option(c("--diff_phos_path"),
              type = "character",
              default = NULL,
              help = "Path to differential phosphoprotein expression results file",
              metavar = "FILE"),
  
  # Output Parameters
  make_option(c("-o", "--outdir"),
              type = "character",
              default = getwd(),
              help = "Output directory for analysis results [default: current directory]",
              metavar = "DIR"),
  
  #---------------------------------------------------------------------------
  # Dataset Names
  #---------------------------------------------------------------------------
  make_option(c("--omics1_name"),
              type = "character",
              default = "Pro",
              help = "Name for first omics dataset (proteomics) [default: %default]"),
  
  make_option(c("--omics2_name"),
              type = "character",
              default = "Phos",
              help = "Name for second omics dataset (phosphoproteomics) [default: %default]"),
  
  make_option(c("--group_comparing"),
              type = "character",
              default = "T:N",
              metavar = "character",
              help = "Comparison groups in format 'case:control' (e.g., 'T:N') [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Analysis Parameters
  #---------------------------------------------------------------------------
  make_option(c("--filter_num"),
              type = "integer",
              default = 300,
              help = "Number of top variable features to retain [default: %default]",
              metavar = "INT"),
  
  make_option(c("--FC_threshold"),
              type = "double",
              default = 1.5,
              help = "Fold change threshold for differential expression [default: %default]",
              metavar = "FLOAT"),
  
  make_option(c("--p_threshold"),
              type = "double",
              default = 0.05,
              help = "P-value threshold for statistical significance [default: %default]",
              metavar = "FLOAT"),
  
  make_option(c("--p_value_type"),
              type = "character",
              default = "q_value",
              help = "Type of p-value to use [default: %default], options: p_value, q_value",
              metavar = "STRING"),
  
  make_option(c("--nBoots"),
              type = "integer",
              default = 20,
              help = "Number of bootstrap iterations for stability analysis [default: %default]",
              metavar = "INT"),
  
  make_option(c("--bootnet_R_threshold"),
              type = "double",
              default = 0.3,
              help = "R threshold for bootstrap network stability [default: %default]",
              metavar = "FLOAT"),
  
  make_option(c("--nCores"),
              type = "integer",
              default = 6,
              help = "Number of CPU cores to use for parallel processing [default: %default]",
              metavar = "INT"),
  
  make_option(c("--stability_threshold"),
              type = "double",
              default = 0.4,
              help = "Stability threshold for edge selection [default: %default]",
              metavar = "FLOAT"),
  
  make_option(c("--cor_method"),
              type = "character",
              default = "spearman",
              help = "Correlation method [default: %default], options: pearson, spearman, kendall",
              metavar = "METHOD"),
  
  make_option(c("--edge_FC_threshold"),
              type = "double",
              default = 1.2,
              help = "Fold change threshold for differential edges [default: %default]",
              metavar = "FLOAT"),
  
  make_option(c("--edge_p_threshold"),
              type = "double",
              default = 0.05,
              help = "P-value threshold for differential edges [default: %default]",
              metavar = "FLOAT"),
  
  make_option(c("--max_subnet_num"),
              type = "integer",
              default = 6,
              help = "Maximum number of subnetworks to display [default: %default]",
              metavar = "INT"),
  
  make_option(c("--enrich_fromType"),
              type = "character",
              default = "UNIPROT",
              help = "ID type for enrichment analysis [default: %default], options: UNIPROT, SYMBOL",
              metavar = "STRING"),
  
  make_option(c("--R_threshold"),
              type = "double",
              default = 0.5,
              help = "Correlation threshold for network edges [default: %default]",
              metavar = "FLOAT"),
  
  #---------------------------------------------------------------------------
  # Visualization Color Parameters
  #---------------------------------------------------------------------------
  make_option(c("--edge_color_pos"),
              type = "character",
              default = "#9b6a65",
              metavar = "character",
              help = "Color for positive edges in differential network [Default: %default]"),
  
  make_option(c("--edge_color_neg"),
              type = "character",
              default = "#5d8992",
              metavar = "character",
              help = "Color for negative edges in differential network [Default: %default]"),
  
  make_option(c("--Enhanced_in_N"),
              type = "character",
              default = "#5d8992",
              metavar = "character",
              help = "Color for nodes enhanced in normal group [Default: %default]"),
  
  make_option(c("--Enhanced_in_T"),
              type = "character",
              default = "#9b6a65",
              metavar = "character",
              help = "Color for nodes enhanced in tumor group [Default: %default]"),
  
  make_option(c("--Only_in_N"),
              type = "character",
              default = "#0c2b32",
              metavar = "character",
              help = "Color for nodes only present in normal group [Default: %default]"),
  
  make_option(c("--Only_in_T"),
              type = "character",
              default = "#381512",
              metavar = "character",
              help = "Color for nodes only present in tumor group [Default: %default]"),
  
  make_option(c("--Conflict_relation"),
              type = "character",
              default = "#808080",
              metavar = "character",
              help = "Color for edges showing conflicting correlations between groups (i.e., positive in one group and negative in the other) [Default: %default]"),
  
  make_option(c("--fill_gradientn_color"),
              type = "character",
              default = "#175663;#dce6e9;#90362d",
              metavar = "character",
              help = "Gradient colors for fill (semicolon-separated) [Default: %default]"),
  
  make_option(c("--color_gradient_low"),
              type = "character",
              default = "#175663",
              metavar = "character",
              help = "Color for low end of enrichment gradient [default: %default]"),
  
  make_option(c("--color_gradient_high"),
              type = "character",
              default = "#90362d",
              metavar = "character",
              help = "Color for high end of enrichment gradient [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Verbose Mode
  #---------------------------------------------------------------------------
  make_option(c("-v", "--verbose"),
              action = "store_true",
              default = FALSE,
              help = "Display verbose runtime information [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                           description = "Stable Differential Network Analysis for Multi-omics Data")
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
#'   - Validates categorical parameters against allowed values
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
      stop("❌ Error: Invalid ", option_name, " value '", value, 
           "'. Please use one of the following options: ",
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
    if (!is.numeric(value)) {
      stop("❌ Error: ", option_name, " must be numeric")
    }
    if (value < min_val || value > max_val) {
      stop("❌ Error: ", option_name, " must be within the range [", 
           min_val, ", ", max_val, "].")
    }
  }
  
  #' Parse color vector from semicolon-separated string
  #' @param color_str String like "#175663;#dce6e9;#90362d"
  #' @return Character vector of colors
  parse_color_vector <- function(color_str) {
    # Remove c() wrapper and whitespace if present
    color_str <- gsub("^c\\(|\\)$", "", color_str)
    color_str <- gsub("\\s+", "", color_str)
    
    # Split by comma or semicolon and remove quotes
    colors <- strsplit(color_str, ",|;")[[1]]
    colors <- gsub('\"|\'', '', colors)
    
    return(colors)
  }
  
  #---------------------------------------------------------------------------
  # Path Cleaning
  #---------------------------------------------------------------------------
  opt$script_path <- clean_path(opt$script_path)
  opt$profile <- clean_path(opt$profile)
  opt$phosfile <- clean_path(opt$phosfile)
  opt$sample_group <- clean_path(opt$sample_group)
  opt$pro_phos_cor <- clean_path(opt$pro_phos_cor)
  opt$diff_pro_path <- clean_path(opt$diff_pro_path)
  opt$diff_phos_path <- clean_path(opt$diff_phos_path)
  opt$outdir <- clean_path(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Input File Validation
  #---------------------------------------------------------------------------
  check_input_file(opt$profile, "profile")
  check_input_file(opt$phosfile, "phosfile")
  check_input_file(opt$sample_group, "sample_group")
  check_input_file(opt$pro_phos_cor, "pro_phos_cor")
  check_input_file(opt$diff_pro_path, "diff_pro_path")
  check_input_file(opt$diff_phos_path, "diff_phos_path")
  
  # Create output directory
  opt$outdir <- check_output_dir(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Define Allowed Values for Categorical Parameters
  #---------------------------------------------------------------------------
  allowed_cor_method <- c("pearson", "spearman", "kendall")
  allowed_p_value_type <- c("p_value", "q_value")
  allowed_enrich_fromType <- c('UNIPROT', 'SYMBOL')
  
  #---------------------------------------------------------------------------
  # Validate Categorical Parameters
  #---------------------------------------------------------------------------
  validate_option_choices(opt$cor_method, allowed_cor_method, "cor_method")
  validate_option_choices(opt$p_value_type, allowed_p_value_type, "p_value_type")
  validate_option_choices(opt$enrich_fromType, allowed_enrich_fromType, "enrich_fromType")
  
  #---------------------------------------------------------------------------
  # Validate Numeric Parameters
  #---------------------------------------------------------------------------
  validate_numeric_range(opt$filter_num, 0, 10000, "filter_num")
  validate_numeric_range(opt$FC_threshold, 0, 10, "FC_threshold")
  validate_numeric_range(opt$nBoots, 0, 100, "nBoots")
  validate_numeric_range(opt$bootnet_R_threshold, 0, 1, "bootnet_R_threshold")
  validate_numeric_range(opt$nCores, 0, 100, "nCores")
  validate_numeric_range(opt$stability_threshold, 0, 1, "stability_threshold")
  validate_numeric_range(opt$edge_FC_threshold, 0, 10, "edge_FC_threshold")
  validate_numeric_range(opt$edge_p_threshold, 0, 1, "edge_p_threshold")
  validate_numeric_range(opt$max_subnet_num, 0, 100, "max_subnet_num")
  validate_numeric_range(opt$R_threshold, 0, 1, "R_threshold")
  
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

# Assign parameters to variables with meaningful names
if (!is.null(opt$script_path)) script_path <- opt$script_path
if (!is.null(opt$profile)) pro_path <- opt$profile
if (!is.null(opt$phosfile)) phos_path <- opt$phosfile
if (!is.null(opt$pro_phos_cor)) phos_pro_path <- opt$pro_phos_cor
if (!is.null(opt$sample_group)) samplelist_path <- opt$sample_group
if (!is.null(opt$diff_pro_path)) diff_pro_path <- opt$diff_pro_path
if (!is.null(opt$diff_phos_path)) diff_phos_path <- opt$diff_phos_path
if (!is.null(opt$outdir)) outdir <- opt$outdir
if (!is.null(opt$group_comparing)) group_comparing <- opt$group_comparing
if (!is.null(opt$omics1_name)) omics1_name <- opt$omics1_name
if (!is.null(opt$omics2_name)) omics2_name <- opt$omics2_name
if (!is.null(opt$filter_num)) filter_num <- opt$filter_num
if (!is.null(opt$FC_threshold)) FC_threshold <- opt$FC_threshold
if (!is.null(opt$p_threshold)) p_threshold <- opt$p_threshold
if (!is.null(opt$p_value_type)) p_value_type <- opt$p_value_type
if (!is.null(opt$nBoots)) nBoots <- opt$nBoots
if (!is.null(opt$bootnet_R_threshold)) bootnet_R_threshold <- opt$bootnet_R_threshold
if (!is.null(opt$nCores)) nCores <- opt$nCores
if (!is.null(opt$stability_threshold)) stability_threshold <- opt$stability_threshold
if (!is.null(opt$cor_method)) cor_method <- opt$cor_method
if (!is.null(opt$edge_FC_threshold)) edge_FC_threshold <- opt$edge_FC_threshold
if (!is.null(opt$edge_p_threshold)) edge_p_threshold <- opt$edge_p_threshold
if (!is.null(opt$max_subnet_num)) max_subnet_num <- opt$max_subnet_num
if (!is.null(opt$R_threshold)) R_threshold <- opt$R_threshold
if (!is.null(opt$enrich_fromType)) enrich_fromType <- opt$enrich_fromType
if (!is.null(opt$edge_color_pos)) edge_color_pos <- opt$edge_color_pos
if (!is.null(opt$edge_color_neg)) edge_color_neg <- opt$edge_color_neg
if (!is.null(opt$Enhanced_in_N)) Enhanced_in_N <- opt$Enhanced_in_N
if (!is.null(opt$Enhanced_in_T)) Enhanced_in_T <- opt$Enhanced_in_T
if (!is.null(opt$color_gradient_low)) color_gradient_low <- opt$color_gradient_low
if (!is.null(opt$color_gradient_high)) color_gradient_high <- opt$color_gradient_high
if (!is.null(opt$Only_in_N)) Only_in_N <- opt$Only_in_N
if (!is.null(opt$Only_in_T)) Only_in_T <- opt$Only_in_T
if (!is.null(opt$Conflict_relation)) Conflict_relation <- opt$Conflict_relation
if (!is.null(opt$fill_gradientn_color)) {
  fill_gradientn_color <- unlist(strsplit(opt$fill_gradientn_color, ';')[[1]])
}

#===============================================================================
# LOAD SOURCE SCRIPTS
#===============================================================================

# Load all required function scripts
source(file.path(script_path, "run_color.R"))
source(file.path(script_path, "run_predata.R"))
source(file.path(script_path, "run_diff.R"))
source(file.path(script_path, "run_corStability.R"))
source(file.path(script_path, "run_conditional_network.R"))
source(file.path(script_path, "run_diff_network.R"))
source(file.path(script_path, "run_enrichment.R"))
source(file.path(script_path, "run_cluster.R"))
source(file.path(script_path, "differential_network.R"))
source(file.path(script_path, "run_mediation.R"))
source(file.path(script_path, "network_show.R"))
source(file.path(script_path, "run_prenetwork.R"))
source(file.path(script_path, "run_diff_enrichment.R"))
source(file.path(script_path, "pipline_save.R"))
source(file.path(script_path, "run_diffsubnet_enrichment.R"))
source(file.path(script_path, "differential_subnetwork_plot.R"))
source(file.path(script_path, "diff_net_community_detection_plot.R"))
source(file.path(script_path, "..", "3.2functional_enrichment_function.R"))

#===============================================================================
# LOAD REQUIRED LIBRARIES
#===============================================================================

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(readr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

#===============================================================================
# DATA LOADING AND PREPROCESSING
#===============================================================================

# Load mapping file between phosphosites and proteins
phos_pro <- read.csv(phos_pro_path, sep = '\t')

# Load proteomics and phosphoproteomics data
pro <- read.csv(pro_path, sep = '\t')
phos <- read.csv(phos_path, sep = '\t')

#---------------------------------------------------------------------------
# Filter features based on differential expression results (if provided)
#---------------------------------------------------------------------------

# Filter proteomics data to keep only differentially expressed proteins
if (diff_pro_path != "") {
  diff_pro <- read.csv(diff_pro_path, sep = '\t')
  diff_pro_ids <- diff_pro[diff_pro$class != 'Non-significant', ][[omics1_name]]
  pro <- pro[pro[[omics1_name]] %in% diff_pro_ids, ]
}

# Filter phosphoproteomics data to keep only differentially expressed phosphosites
if (diff_phos_path != "") {
  diff_phos <- read.csv(diff_phos_path, sep = '\t')
  diff_phos_ids <- diff_phos[diff_phos$class != 'Non-significant', ][[omics2_name]]
  phos <- phos[phos[[omics2_name]] %in% diff_phos_ids, ]
}

#---------------------------------------------------------------------------
# Standardize column names
#---------------------------------------------------------------------------
colnames(pro)[1] <- "feature_ID"
colnames(phos)[1] <- "feature_ID"

#---------------------------------------------------------------------------
# Combine proteomics and phosphoproteomics data
#---------------------------------------------------------------------------
count_table <- rbind(pro, phos)
print(paste("Total number of features:", nrow(count_table)))

#---------------------------------------------------------------------------
# Create annotation table with omics class information
#---------------------------------------------------------------------------
pro_annotation <- pro %>%
  dplyr::select(feature_ID) %>%
  dplyr::mutate(Class = omics1_name, KEGG.ID = NA)

phos_annotation <- phos %>%
  dplyr::select(feature_ID) %>%
  dplyr::mutate(Class = omics2_name, KEGG.ID = NA)

annotation_table <- rbind(pro_annotation, phos_annotation)
annotation_table$omics_name <- annotation_table$Class

#---------------------------------------------------------------------------
# Load sample group information
#---------------------------------------------------------------------------
samplelist <- read.csv(samplelist_path, sep = "\t")

#===============================================================================
# MAIN DIFFERENTIAL NETWORK ANALYSIS
#===============================================================================

#' Run the complete differential network analysis pipeline
#'
#' This function orchestrates the entire differential network analysis workflow:
#'   - Network construction for each condition
#'   - Stability analysis via bootstrap
#'   - Differential network identification
#'   - Functional enrichment analysis
#'   - Mediation analysis (optional)
#'
#' @param count_table Combined abundance data for all features
#' @param samplelist Sample grouping information
#' @param compare_group Comparison groups (e.g., "T:N")
#' @param filter_num Number of top variable features to retain
#' @param annotation_table Feature annotations with omics class
#' @param FC_threshold Fold change threshold for differential expression
#' @param p_threshold P-value threshold for statistical significance
#' @param p_value_type Type of p-value to use (p_value or q_value)
#' @param nBoots Number of bootstrap iterations
#' @param bootnet_R_threshold R threshold for bootstrap network stability
#' @param nCores Number of CPU cores for parallel processing
#' @param stability_threshold Stability threshold for edge selection
#' @param cor_method Correlation method
#' @param edge_FC_threshold Fold change threshold for differential edges
#' @param edge_p_threshold P-value threshold for differential edges
#' @param run_enrich Whether to run enrichment analysis
#' @param run_mediation Whether to run mediation analysis
#' @param species Species for enrichment analysis
#' @param run_diffsubnet_enrich Whether to run differential subnetwork enrichment
#'
#' @return Network analysis results object
differential_network1 <- differential_network(
  count_table = count_table,
  samplelist = samplelist,
  compare_group = group_comparing,
  filter_num = filter_num,
  annotation_table = annotation_table,
  FC_threshold = FC_threshold,
  p_threshold = p_threshold,
  p_value_type = p_value_type,
  nBoots = nBoots,
  bootnet_R_threshold = bootnet_R_threshold,
  nCores = nCores,
  stability_threshold = stability_threshold,
  cor_method = "spearman",
  edge_FC_threshold = edge_FC_threshold,
  edge_p_threshold = edge_p_threshold,
  run_enrich = FALSE,
  run_mediation = FALSE,
  species = "hsa",
  run_diffsubnet_enrich = FALSE
)

#===============================================================================
# SAVE PIPELINE RESULTS
#===============================================================================

#' Save all pipeline results to files and generate visualizations
#'
#' This function saves the complete analysis results including:
#'   - Network objects
#'   - Subnetworks
#'   - Enrichment results
#'   - Visualizations
#'
#' @param Network Network analysis results object
#' @param outdir Output directory
#' @param R_threshold Correlation threshold for network edges
#' @param max_subnet_num Maximum number of subnetworks to display
#' @param omics1_name Name of proteomics dataset
#' @param omics2_name Name of phosphoproteomics dataset
#' @param phos_pro Protein-phosphosite mapping
#' @param enrich_fromType ID type for enrichment analysis
#' @param edge_color_pos Color for positive edges
#' @param edge_color_neg Color for negative edges
#' @param Enhanced_in_N Color for nodes enhanced in normal group
#' @param Enhanced_in_T Color for nodes enhanced in tumor group
#' @param Only_in_N Color for nodes only in normal group
#' @param Only_in_T Color for nodes only in tumor group
#' @param Conflict_relation Color for conflicting edges
#' @param color_gradient_low Low end color for enrichment gradient
#' @param color_gradient_high High end color for enrichment gradient
#' @param fill_gradientn_color Gradient colors for fill
pipline_save(
  Network = differential_network1,
  outdir = outdir,
  R_threshold = R_threshold,
  max_subnet_num = max_subnet_num,
  omics1_name = omics1_name,
  omics2_name = omics2_name,
  phos_pro = phos_pro,
  enrich_fromType = enrich_fromType,
  edge_color_pos = edge_color_pos,
  edge_color_neg = edge_color_neg,
  Enhanced_in_N = Enhanced_in_N,
  Enhanced_in_T = Enhanced_in_T,
  Only_in_N = Only_in_N,
  Only_in_T = Only_in_T,
  Conflict_relation = Conflict_relation,
  color_gradient_low = color_gradient_low,
  color_gradient_high = color_gradient_high,
  fill_gradientn_color = fill_gradientn_color
)