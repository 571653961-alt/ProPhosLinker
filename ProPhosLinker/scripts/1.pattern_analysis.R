#!/usr/bin/env Rscript

#===============================================================================
# Pattern Analysis Main Script
#===============================================================================
# This script orchestrates the complete pattern analysis workflow for multi-omics
# data, including:
#   1.1 Procrustes Analysis (PA) - assessing consistency between datasets
#   1.2 Non-negative Matrix Factorization (NMF) - molecular subtyping
#   1.3 Weighted Gene Co-expression Network Analysis (WGCNA) - module detection
#
# The script accepts command-line arguments via optparse and executes all three
# analyses sequentially, generating comprehensive output files and visualizations.
#===============================================================================

suppressWarnings(library(optparse))

#===============================================================================
# COMMAND LINE ARGUMENT DEFINITION
#===============================================================================

# Define all command-line options
option_list <- list(
  #---------------------------------------------------------------------------
  # Script and Input/Output Paths
  #---------------------------------------------------------------------------
  make_option(
    c("--script_path"),
    type = "character",
    default = NULL,
    help = "Path to the directory containing analysis scripts",
    metavar = "FILE"
  ),
  make_option(
    c("--profile"),
    type = "character",
    default = NULL,
    help = "Path to normalized protein abundance data file (TSV format)",
    metavar = "FILE"
  ),
  make_option(
    c("--phosfile"),
    type = "character",
    default = NULL,
    help = "Path to normalized phosphoprotein abundance data file (TSV format)",
    metavar = "FILE"
  ),
  make_option(
    c("--sample_group"),
    type = "character",
    default = NULL,
    help = "Path to sample group information file (TSV format)",
    metavar = "FILE"
  ),
  make_option(
    c("--metadatafile"),
    type = "character",
    default = NULL,
    help = "Path to metadata/clinical data file (TSV format)",
    metavar = "FILE"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = getwd(),
    help = "Output directory for analysis results [default: current working directory]",
    metavar = "DIR"
  ),
  
  #---------------------------------------------------------------------------
  # Dataset Names
  #---------------------------------------------------------------------------
  make_option(
    c("--omics1_name"),
    type = "character",
    default = "Pro",
    help = "Name for first omics dataset (e.g., Proteomics) [default: %default]"
  ),
  make_option(
    c("--omics2_name"),
    type = "character",
    default = "Phos",
    help = "Name for second omics dataset (e.g., Phosphoproteomics) [default: %default]"
  ),
  
  #---------------------------------------------------------------------------
  # 1.1 Procrustes Analysis Parameters
  #---------------------------------------------------------------------------
  make_option(
    c("--dim_rd_mtd"),
    type = "character",
    default = "PCA",
    help = "Dimensionality reduction method for Procrustes [default: %default], options: PCA, PCOA",
    metavar = "METHOD"
  ),
  make_option(
    c("--PA_group1_color"),
    type = "character",
    default = "#a03c32",
    help = "Color for first omics dataset in Procrustes plots [default: %default]"
  ),
  make_option(
    c("--PA_group2_color"),
    type = "character",
    default = "#43656C",
    help = "Color for second omics dataset in Procrustes plots [default: %default]"
  ),
  
  #---------------------------------------------------------------------------
  # 1.2 NMF (Molecular Subtyping) Parameters
  #---------------------------------------------------------------------------
  make_option(
    c("--NMF_pro_filter_num"),
    type = "integer",
    default = 3000,
    help = "Number of top variable features to retain for proteomics NMF [default: %default]",
    metavar = "INT"
  ),
  make_option(
    c("--NMF_phos_filter_num"),
    type = "integer",
    default = 3000,
    help = "Number of top variable features to retain for phosphoproteomics NMF [default: %default]",
    metavar = "INT"
  ),
  
  #---------------------------------------------------------------------------
  # 1.3 WGCNA Parameters - Feature Filtering
  #---------------------------------------------------------------------------
  make_option(
    c("--WGCNA_pro_filter_num"),
    type = "integer",
    default = 5000,
    help = "Number of top ANOVA features to retain for proteomics WGCNA [default: %default]",
    metavar = "INT"
  ),
  make_option(
    c("--WGCNA_phos_filter_num"),
    type = "integer",
    default = 5000,
    help = "Number of top ANOVA features to retain for phosphoproteomics WGCNA [default: %default]",
    metavar = "INT"
  ),
  
  #---------------------------------------------------------------------------
  # 1.3 WGCNA Parameters - Proteomics
  #---------------------------------------------------------------------------
  make_option(
    c("--protein_SoftPower"),
    type = "integer",
    default = NULL,
    help = "Soft threshold power for proteomics (NULL = auto-select) [default: %default]",
    metavar = "INT"
  ),
  make_option(
    c("--protein_RsquareCut"),
    type = "double",
    default = 0.890,
    help = "R² cutoff for scale-free topology fit (proteomics) [default: %default]",
    metavar = "FLOAT"
  ),
  make_option(
    c("--protein_cor_method"),
    type = "character",
    default = "spearman",
    help = "Correlation method for proteomics [default: %default], options: pearson, spearman, kendall, bicor",
    metavar = "METHOD"
  ),
  make_option(
    c("--protein_corFun_tmp"),
    type = "character",
    default = "bicor",
    help = "Correlation function for proteomics [default: %default], options: bicor, cor",
    metavar = "STRING"
  ),
  make_option(
    c("--protein_cluster_method"),
    type = "character",
    default = "average",
    help = "Hierarchical clustering method for proteomics [default: %default], options: average, complete, ward.D, ward.D2, single, mcquitty, median, centroid",
    metavar = "STRING"
  ),
  make_option(
    c("--protein_corOptions"),
    type = "character",
    default = "pairwise.complete.obs",
    help = "Correlation options for proteomics [default: '%default'], options: pairwise.complete.obs, everything, all.obs, complete.obs, na.or.complete",
    metavar = "STRING"
  ),
  make_option(
    c("--protein_networkType"),
    type = "character",
    default = "signed",
    help = "Network type for proteomics [default: %default], options: signed, unsigned, signed hybrid",
    metavar = "STRING"
  ),
  make_option(
    c("--protein_mergingThresh"),
    type = "double",
    default = 0.20,
    help = "Module merging threshold (1-correlation) for proteomics [default: %default]",
    metavar = "FLOAT"
  ),
  make_option(
    c("--protein_minModuleSize"),
    type = "integer",
    default = 60,
    help = "Minimum module size for proteomics [default: %default]",
    metavar = "INT"
  ),
  
  #---------------------------------------------------------------------------
  # 1.3 WGCNA Parameters - Phosphoproteomics
  #---------------------------------------------------------------------------
  make_option(
    c("--phosphoprotein_SoftPower"),
    type = "integer",
    default = NULL,
    help = "Soft threshold power for phosphoproteomics (NULL = auto-select) [default: %default]",
    metavar = "INT"
  ),
  make_option(
    c("--phosphoprotein_RsquareCut"),
    type = "double",
    default = 0.890,
    help = "R² cutoff for scale-free topology fit (phosphoproteomics) [default: %default]",
    metavar = "FLOAT"
  ),
  make_option(
    c("--phosphoprotein_networkType"),
    type = "character",
    default = "signed",
    help = "Network type for phosphoproteomics [default: %default], options: signed, unsigned, signed hybrid",
    metavar = "METHOD"
  ),
  make_option(
    c("--phosphoprotein_cor_method"),
    type = "character",
    default = "spearman",
    help = "Correlation method for phosphoproteomics [default: %default], options: pearson, spearman, kendall",
    metavar = "METHOD"
  ),
  make_option(
    c("--phosphoprotein_corFun_tmp"),
    type = "character",
    default = "bicor",
    help = "Correlation function for phosphoproteomics [default: %default]",
    metavar = "STRING"
  ),
  make_option(
    c("--phosphoprotein_cluster_method"),
    type = "character",
    default = "average",
    help = "Hierarchical clustering method for phosphoproteomics [default: %default], options: average, complete, ward.D, ward.D2, single, mcquitty, median, centroid",
    metavar = "STRING"
  ),
  make_option(
    c("--phosphoprotein_corOptions"),
    type = "character",
    default = "pairwise.complete.obs",
    help = "Correlation options for phosphoproteomics [default: '%default'], options: pairwise.complete.obs, everything, all.obs, complete.obs, na.or.complete",
    metavar = "STRING"
  ),
  make_option(
    c("--phosphoprotein_mergingThresh"),
    type = "double",
    default = 0.20,
    help = "Module merging threshold (1-correlation) for phosphoproteomics [default: %default]",
    metavar = "FLOAT"
  ),
  make_option(
    c("--phosphoprotein_minModuleSize"),
    type = "integer",
    default = 60,
    help = "Minimum module size for phosphoproteomics [default: %default]",
    metavar = "INT"
  ),
  
  #---------------------------------------------------------------------------
  # Module Correlation Parameters
  #---------------------------------------------------------------------------
  make_option(
    c("--module_cor_threshhold"),
    type = "double",
    default = 0.8,
    help = "Correlation threshold for significant module pairs [default: %default]",
    metavar = "NUMERIC"
  ),
  make_option(
    c("--module_cor_p_adj"),
    type = "double",
    default = 0.05,
    help = "Adjusted p-value threshold for module correlations [default: %default]",
    metavar = "NUMERIC"
  ),
  make_option(
    c("--module_cor_method"),
    type = "character",
    default = "spearman",
    help = "Correlation method for module analysis [default: %default], options: pearson, spearman, kendall",
    metavar = "METHOD"
  ),
  
  #---------------------------------------------------------------------------
  # Visualization Color Parameters
  #---------------------------------------------------------------------------
  make_option(
    c("--pro_ME_color"),
    type = "character",
    default = "#D12128",
    help = "Color for proteomics module eigengenes [default: %default]"
  ),
  make_option(
    c("--phos_ME_color"),
    type = "character",
    default = "#01344F",
    help = "Color for phosphoproteomics module eigengenes [default: %default]"
  ),
  make_option(
    c("--pro_color"),
    type = "character",
    default = "#a03c32",
    help = "Color for proteomics data in plots [default: %default]"
  ),
  make_option(
    c("--phos_color"),
    type = "character",
    default = "#43656C",
    help = "Color for phosphoproteomics data in plots [default: %default]"
  ),
  make_option(
    c("--pheatmap_color"),
    type = "character",
    default = 'c("#43656C", "white", "#a03c32")',
    help = "Color palette for pheatmap visualizations [default: %default]"
  ),
  
  #---------------------------------------------------------------------------
  # Verbose Mode
  #---------------------------------------------------------------------------
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Display verbose runtime information [default: %default]"
  )
)

# Create option parser with description
opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nPattern Analysis Pipeline\nIncludes: Procrustes Analysis, Molecular Subtyping (NMF), and WGCNA"
)

# Parse command-line arguments
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
  
  #' Create output directory if it doesn't exist
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
  
  #' Parse color vector from string representation
  #' @param color_str String like 'c("#43656C", "white", "#a03c32")'
  #' @return Character vector of colors
  parse_color_vector <- function(color_str) {
    # Remove c() wrapper and whitespace
    color_str <- gsub("^c\\(|\\)$", "", color_str)
    color_str <- gsub("\\s+", "", color_str)
    
    # Split by comma and remove quotes
    colors <- strsplit(color_str, ",")[[1]]
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
  
  # Handle metadata file (can be 'None')
  if (opt$metadatafile != 'None') {
    opt$metadatafile <- clean_path(opt$metadatafile)
  }
  opt$outdir <- clean_path(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Input File Validation
  #---------------------------------------------------------------------------
  check_input_file(opt$profile, "profile")
  check_input_file(opt$phosfile, "phosfile")
  check_input_file(opt$sample_group, "sample_group")
  
  if (opt$metadatafile != 'None') {
    check_input_file(opt$metadatafile, "metadatafile")
  }
  
  # Create output directory
  opt$outdir <- check_output_dir(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Parse Color Parameters
  #---------------------------------------------------------------------------
  opt$pheatmap_color <- parse_color_vector(opt$pheatmap_color)
  
  #---------------------------------------------------------------------------
  # Define Allowed Values for Categorical Parameters
  #---------------------------------------------------------------------------
  
  # Dimensionality reduction methods
  allowed_dim_rd_mtd <- c("PCA", "PCOA")
  validate_option_choices(opt$dim_rd_mtd, allowed_dim_rd_mtd, "dim_rd_mtd")
  
  # Correlation methods
  allowed_protein_cor_method <- c("pearson", "spearman", "kendall")
  allowed_phosphoprotein_cor_method <- c("pearson", "spearman", "kendall")
  allowed_protein_corFun_tmp <- c("bicor", "cor")
  allowed_phosphoprotein_corFun_tmp <- c("bicor", "cor")
  
  # Clustering methods
  allowed_protein_cluster_method <- c("average", "complete", "ward.D", "ward.D2", 
                                      "single", "mcquitty", "median", "centroid")
  allowed_phosphoprotein_cluster_method <- c("average", "complete", "ward.D", "ward.D2", 
                                             "single", "mcquitty", "median", "centroid")
  
  # Correlation options
  allowed_protein_corOptions <- c("pairwise.complete.obs", "everything", "all.obs", 
                                  "complete.obs", "na.or.complete")
  allowed_phosphoprotein_corOptions <- c("pairwise.complete.obs", "everything", "all.obs", 
                                         "complete.obs", "na.or.complete")
  
  # Network types
  allowed_protein_networkType <- c("signed", "unsigned", "signed hybrid")
  allowed_phosphoprotein_networkType <- c("signed", "unsigned", "signed hybrid")
  
  # Module correlation methods
  allowed_module_cor_method <- c("pearson", "spearman", "kendall")
  
  #---------------------------------------------------------------------------
  # Validate Categorical Parameters
  #---------------------------------------------------------------------------
  validate_option_choices(opt$protein_cor_method, allowed_protein_cor_method, 
                          "protein_cor_method")
  validate_option_choices(opt$phosphoprotein_cor_method, allowed_phosphoprotein_cor_method, 
                          "phosphoprotein_cor_method")
  validate_option_choices(opt$protein_corFun_tmp, allowed_protein_corFun_tmp, 
                          "protein_corFun_tmp")
  validate_option_choices(opt$phosphoprotein_corFun_tmp, allowed_phosphoprotein_corFun_tmp, 
                          "phosphoprotein_corFun_tmp")
  validate_option_choices(opt$protein_cluster_method, allowed_protein_cluster_method, 
                          "protein_cluster_method")
  validate_option_choices(opt$phosphoprotein_cluster_method, allowed_phosphoprotein_cluster_method, 
                          "phosphoprotein_cluster_method")
  validate_option_choices(opt$protein_corOptions, allowed_protein_corOptions, 
                          "protein_corOptions")
  validate_option_choices(opt$phosphoprotein_corOptions, allowed_phosphoprotein_corOptions, 
                          "phosphoprotein_corOptions")
  validate_option_choices(opt$protein_networkType, allowed_protein_networkType, 
                          "protein_networkType")
  validate_option_choices(opt$phosphoprotein_networkType, allowed_phosphoprotein_networkType, 
                          "phosphoprotein_networkType")
  validate_option_choices(opt$module_cor_method, allowed_module_cor_method, 
                          "module_cor_method")
  
  #---------------------------------------------------------------------------
  # Validate Numeric Parameters
  #---------------------------------------------------------------------------
  
  # Feature count parameters
  validate_numeric_range(opt$NMF_pro_filter_num, 0, 100000, "NMF_pro_filter_num")
  validate_numeric_range(opt$NMF_phos_filter_num, 0, 100000, "NMF_phos_filter_num")
  validate_numeric_range(opt$WGCNA_pro_filter_num, 0, 100000, "WGCNA_pro_filter_num")
  validate_numeric_range(opt$WGCNA_phos_filter_num, 0, 100000, "WGCNA_phos_filter_num")
  
  # R-squared cutoffs
  validate_numeric_range(opt$protein_RsquareCut, 0, 1, "protein_RsquareCut")
  validate_numeric_range(opt$phosphoprotein_RsquareCut, 0, 1, "phosphoprotein_RsquareCut")
  
  # Module merging thresholds
  validate_numeric_range(opt$protein_mergingThresh, 0, 1, "protein_mergingThresh")
  validate_numeric_range(opt$phosphoprotein_mergingThresh, 0, 1, "phosphoprotein_mergingThresh")
  
  # Minimum module sizes
  validate_numeric_range(opt$protein_minModuleSize, 0, 1000, "protein_minModuleSize")
  validate_numeric_range(opt$phosphoprotein_minModuleSize, 0, 1000, "phosphoprotein_minModuleSize")
  
  # Correlation threshold
  validate_numeric_range(opt$module_cor_threshhold, 0, 1, "module_cor_threshhold")
  
  #---------------------------------------------------------------------------
  # Verbose Output
  #---------------------------------------------------------------------------
  if (opt$verbose) {
    print(opt)
  }
  
  return(opt)
}

# Execute parameter validation
opt <- parameter_validation(opt)

#===============================================================================
# EXTRACT PARAMETERS FOR ANALYSIS
#===============================================================================

#---------------------------------------------------------------------------
# Path Parameters
#---------------------------------------------------------------------------
script_path <- opt$script_path
pro_path <- opt$profile
phos_path <- opt$phosfile
group_path <- opt$sample_group

# Handle metadata path (can be NULL)
if (opt$metadatafile == 'None') {
  metadata_path <- NULL
} else {
  metadata_path <- opt$metadatafile
}

outdir <- opt$outdir

#---------------------------------------------------------------------------
# Dataset Names
#---------------------------------------------------------------------------
omics1_name <- opt$omics1_name
omics2_name <- opt$omics2_name

#---------------------------------------------------------------------------
# 1.1 Procrustes Analysis Parameters
#---------------------------------------------------------------------------
dim_rd_mtd <- opt$dim_rd_mtd
PA_group1_color <- opt$PA_group1_color
PA_group2_color <- opt$PA_group2_color

#---------------------------------------------------------------------------
# 1.2 NMF Parameters
#---------------------------------------------------------------------------
NMF_pro_filter_num <- opt$NMF_pro_filter_num
NMF_phos_filter_num <- opt$NMF_phos_filter_num

#---------------------------------------------------------------------------
# 1.3 WGCNA Parameters
#---------------------------------------------------------------------------
WGCNA_pro_filter_num <- opt$WGCNA_pro_filter_num
WGCNA_phos_filter_num <- opt$WGCNA_phos_filter_num

# Proteomics parameters
protein_SoftPower <- opt$protein_SoftPower
protein_RsquareCut <- opt$protein_RsquareCut
protein_cor_method <- opt$protein_cor_method
protein_corFun_tmp <- opt$protein_corFun_tmp
protein_cluster_method <- opt$protein_cluster_method
protein_corOptions_str <- paste0("use = '", opt$protein_corOptions, "'")
protein_networkType <- opt$protein_networkType
protein_mergingThresh <- opt$protein_mergingThresh
protein_minModuleSize <- opt$protein_minModuleSize

# Phosphoproteomics parameters
phosphoprotein_SoftPower <- opt$phosphoprotein_SoftPower
phosphoprotein_RsquareCut <- opt$phosphoprotein_RsquareCut
phosphoprotein_networkType <- opt$phosphoprotein_networkType
phosphoprotein_cor_method <- opt$phosphoprotein_cor_method
phosphoprotein_corFun_tmp <- opt$phosphoprotein_corFun_tmp
phosphoprotein_cluster_method <- opt$phosphoprotein_cluster_method
phosphoprotein_corOptions_str <- paste0("use = '", opt$phosphoprotein_corOptions, "'")
phosphoprotein_mergingThresh <- opt$phosphoprotein_mergingThresh
phosphoprotein_minModuleSize <- opt$phosphoprotein_minModuleSize

# Module correlation parameters
module_cor_threshhold <- opt$module_cor_threshhold
module_cor_p_adj <- opt$module_cor_p_adj
module_cor_method <- opt$module_cor_method

# Visualization colors
pro_ME_color <- opt$pro_ME_color
phos_ME_color <- opt$phos_ME_color
pro_color <- opt$pro_color
phos_color <- opt$phos_color
pheatmap_color <- strsplit(opt$pheatmap_color, ';')[[1]]

#===============================================================================
# LOAD ANALYSIS SCRIPTS
#===============================================================================

source(file.path(script_path, "1.1procrustes.R"))
source(file.path(script_path, "1.2clusterbyNMF.R"))
source(file.path(script_path, "1.3WGCNA.R"))

#===============================================================================
# 1.1 PROCRUSTES ANALYSIS
#===============================================================================

cat("\n========================================================================\n")
cat("Running 1.1 Procrustes Analysis...\n")
cat("========================================================================\n")

# Create output directory for Procrustes results
PA_outdir <- file.path(outdir, '1.1Procrustes')
if (dir.exists(PA_outdir)) {
  unlink(PA_outdir, recursive = TRUE)
  dir.create(PA_outdir)
} else {
  dir.create(PA_outdir)
}

# Execute Procrustes Analysis
Procrustes_Analysis(
  pro_path,
  phos_path,
  group_path,
  PA_outdir,
  dim_rd_mtd,
  omics1_name,
  omics2_name,
  PA_group1_color = PA_group1_color,
  PA_group2_color = PA_group2_color
)

cat("✅ Procrustes Analysis completed. Results saved to:", PA_outdir, "\n")

#===============================================================================
# 1.2 NMF MOLECULAR SUBTYPING
#===============================================================================

cat("\n========================================================================\n")
cat("Running 1.2 NMF Molecular Subtyping...\n")
cat("========================================================================\n")

# Create output directory for NMF results
NMF_outdir <- file.path(outdir, '1.2Molecular_Subtyping')
if (dir.exists(NMF_outdir)) {
  unlink(NMF_outdir, recursive = TRUE)
  dir.create(NMF_outdir)
} else {
  dir.create(NMF_outdir)
}

# Load data for NMF analysis
pro <- read_tsv(pro_path)
phos <- read_tsv(phos_path)
group_info <- read_tsv(group_path)

# Extract sample names (first column is feature names)
name <- colnames(pro)
name <- name[2:length(name)]

# Set filter parameters
pro_filter_top_var_num <- NMF_pro_filter_num
phos_filter_top_var_num <- NMF_phos_filter_num

# Perform NMF clustering for both datasets
cat("  - Performing NMF on proteomics data...\n")
res_pro <- nmf_cluster(pro, pro_filter_top_var_num)

cat("  - Performing NMF on phosphoproteomics data...\n")
res_phos <- nmf_cluster(phos, phos_filter_top_var_num)

# Create data frame with cluster assignments
clusters_df <- tibble(
  sample = name,  # Assume same sample order for both datasets
  Pro_cluster = paste0(omics1_name, "_", res_pro$sample_cluster),
  Phos_cluster = paste0(omics2_name, "_", res_phos$sample_cluster)
)

# Save cluster assignments
write_tsv(clusters_df, file.path(NMF_outdir, "molecular_subtyping_clusters.tsv"))

# Generate and save Sankey diagram
plt <- plotSankey(clusters_df, omics1_name = omics1_name, omics2_name = omics2_name)
ggsave(file.path(NMF_outdir, "molecular_subtyping_sankey.png"), plt, width = 8, height = 6)

#---------------------------------------------------------------------------
# Group-specific NMF analysis
#---------------------------------------------------------------------------
cat("  - Generating group-specific cluster assignments...\n")

# Create data frame with sample, proteomics cluster, and phosphoproteomics cluster
clusters_df_by_group <- tibble(
  sample = name,
  Pro_cluster = paste0(omics1_name, "_", res_pro$sample_cluster),
  Phos_cluster = paste0(omics2_name, "_", res_phos$sample_cluster)
)

# Add group information
clusters_df_by_group <- clusters_df_by_group %>%
  left_join(group_info, by = "sample")

# Create combined group-cluster identifiers
clusters_df_by_group2 <- tibble(
  sample = clusters_df_by_group$sample,
  group_Pro_cluster = paste0(clusters_df_by_group$group, "_", clusters_df_by_group$Pro_cluster),
  group_Phos_cluster = paste0(clusters_df_by_group$group, "_", clusters_df_by_group$Phos_cluster)
)

# Save group-specific cluster assignments
write_tsv(clusters_df_by_group2, file.path(NMF_outdir, "molecular_subtyping_clusters_by_group.tsv"))

# Generate and save group-specific Sankey diagram
plt_by_group <- plotSankey_by_group(clusters_df_by_group2, 
                                    omics1_name = omics1_name, 
                                    omics2_name = omics2_name)
ggsave(file.path(NMF_outdir, "molecular_subtyping_sankey_by_group.png"), 
       plt_by_group, width = 8, height = 6)

cat("✅ NMF Molecular Subtyping completed. Results saved to:", NMF_outdir, "\n")

#===============================================================================
# 1.3 WGCNA ANALYSIS
#===============================================================================

cat("\n========================================================================\n")
cat("Running 1.3 WGCNA Analysis...\n")
cat("========================================================================\n")

# Create output directory for WGCNA results
WGCNA_outdir <- file.path(outdir, '1.3WGCNA')
if (dir.exists(WGCNA_outdir)) {
  unlink(WGCNA_outdir, recursive = TRUE)
  dir.create(WGCNA_outdir)
} else {
  dir.create(WGCNA_outdir)
}

# Execute WGCNA Analysis
WGCNA_Analysis(
  pro_path = pro_path,
  phos_path = phos_path,
  sample_group = group_path,
  metadata_path = metadata_path,
  out_dir = WGCNA_outdir,
  omics1_name = omics1_name,
  omics2_name = omics2_name,
  module_cor_threshhold = module_cor_threshhold,
  module_cor_p_adj = module_cor_p_adj,
  module_cor_method = module_cor_method,
  WGCNA_pro_filter_num = WGCNA_pro_filter_num,
  WGCNA_phos_filter_num = WGCNA_phos_filter_num,
  pro_ME_color = pro_ME_color,
  phos_ME_color = phos_ME_color,
  pro_color = pro_color,
  phos_color = phos_color,
  pheatmap_color = pheatmap_color
)

cat("✅ WGCNA Analysis completed. Results saved to:", WGCNA_outdir, "\n")
cat("\n========================================================================\n")
cat("🎉 All pattern analyses completed successfully!\n")
cat("========================================================================\n")