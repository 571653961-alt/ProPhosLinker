#===============================================================================
# NMF-based Clustering and Visualization Functions
#===============================================================================
# This script contains functions for performing Non-negative Matrix Factorization
# (NMF) clustering on multi-omics data and visualizing the results using Sankey
# diagrams to show cluster concordance between different omics datasets.
#===============================================================================

# Load necessary libraries with suppressed startup messages
suppressWarnings(suppressPackageStartupMessages(library(NMF)))        # For NMF clustering
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))  # For data manipulation
suppressWarnings(suppressPackageStartupMessages(library(ggsankeyfier)))# For Sankey diagram visualization
suppressWarnings(suppressPackageStartupMessages(library(doParallel)))  # For parallel computing

#===============================================================================
# Core NMF Clustering Function
#===============================================================================

#' Perform NMF clustering on omics data
#'
#' This function performs Non-negative Matrix Factorization (NMF) clustering on
#' omics data, including feature filtering, rank estimation, and extraction of
#' sample clusters and feature programs.
#'
#' @param df Data frame with first column as feature identifiers (e.g., gene/protein names)
#'           and remaining columns as samples with expression values
#' @param filter_top_var_num Integer, number of top variable features to retain (default: 3000)
#' @param cutoff Numeric, proportion of top features to retain based on SD (optional)
#' @param parallel Logical, whether to use parallel computing (default: FALSE)
#' @param reserved_cores Integer, number of cores to reserve when using parallel computing (default: 1)
#' @param top_features_per_cluster Integer, number of top features to extract per cluster (default: 50)
#'
#' @return List containing:
#'   \item{nmf}{NMF object from the nmf() function}
#'   \item{W}{Basis matrix (features × clusters)}
#'   \item{H}{Coefficient matrix (clusters × samples)}
#'   \item{sample_cluster}{Factor vector assigning each sample to a cluster}
#'   \item{feature_programs}{List of top features for each cluster}
#' @export
#'
#' @examples
#' \dontrun{
#' result <- nmf_cluster(
#'   df = expression_data,
#'   filter_top_var_num = 2000,
#'   parallel = TRUE
#' )
#' }
nmf_cluster <- function(
    df,
    filter_top_var_num = 3000,
    cutoff = NULL,
    parallel = FALSE,
    reserved_cores = 1,
    top_features_per_cluster = 50
) {
  #---------------------------------------------------------------------------
  # Parallel computing setup (if enabled)
  #---------------------------------------------------------------------------
  if (isTRUE(parallel)) {
    # Create cluster using all available cores except reserved ones
    cl <- makeCluster(max(1, parallel::detectCores() - max(1, reserved_cores)))
    # Ensure cluster is stopped on function exit
    on.exit({
      if (!is.null(cl) && inherits(cl, "cluster")) {
        parallel::stopCluster(cl)
      }
    }, add = TRUE)
    doParallel::registerDoParallel(cl)
  }
  
  #---------------------------------------------------------------------------
  # Data preprocessing
  #---------------------------------------------------------------------------
  # Filter features by standard deviation if cutoff is provided
  data <- if (!is.null(cutoff)) filter_top_sd(df, cutoff) else df
  
  # Convert data frame to matrix format
  omics_mat <- matrix_from_df(data)
  
  # Filter top variable features
  omics_mat <- filter_top_var(omics_mat, filter_top_var_num)
  
  # Add small constant to avoid zero values (NMF requires non-negative data)
  omics_mat <- omics_mat + 1e-7
  
  #---------------------------------------------------------------------------
  # NMF analysis
  #---------------------------------------------------------------------------
  # Estimate optimal number of clusters (rank)
  k_omics <- estimate_rank(omics_mat)
  
  # Perform NMF with estimated rank
  nmf_res <- run_nmf(omics_mat, k_omics)
  
  #---------------------------------------------------------------------------
  # Extract results
  #---------------------------------------------------------------------------
  ## Extract basis and coefficient matrices
  # W: basis matrix (features × clusters) - represents feature programs
  W <- NMF::basis(nmf_res)
  
  # H: coefficient matrix (clusters × samples) - represents sample membership strengths
  H <- NMF::coef(nmf_res)
  
  ## Assign samples to clusters based on maximum coefficient
  sample_cluster <- apply(H, 2, which.max)
  sample_cluster <- factor(
    sample_cluster,
    levels = seq_len(k_omics),
    labels = paste0("Cluster_", seq_len(k_omics))
  )
  
  ## Extract top features for each cluster
  feature_programs <- lapply(seq_len(k_omics), function(i) {
    # Get feature weights for cluster i
    w_i <- W[, i]
    # Sort in descending order
    w_i <- sort(w_i, decreasing = TRUE)
    # Keep top features
    head(w_i, top_features_per_cluster)
  })
  names(feature_programs) <- paste0("Cluster_", seq_len(k_omics))
  
  # Return comprehensive results
  return(list(
    nmf = nmf_res,
    W = W,
    H = H,
    sample_cluster = sample_cluster,
    feature_programs = feature_programs
  ))
}

#===============================================================================
# Helper Functions
#===============================================================================

#' Filter top features by standard deviation
#'
#' @param data Data frame with features in rows and samples in columns
#' @param ratio Numeric, proportion of top features to retain based on SD
#' @return Filtered data frame with top variable features
filter_top_sd <- function(data, ratio) {
  data |>
    # Calculate standard deviation for each feature (excluding first column)
    mutate(sd = apply(across(-1), 1, sd)) |>
    # Keep top proportion of features with highest SD
    slice_max(sd, prop = ratio) |>
    # Remove SD column
    select(-sd)
}

#' Run NMF with specified parameters
#'
#' @param mat Numeric matrix with non-negative values (features × samples)
#' @param k Integer, rank (number of clusters) for NMF
#' @param nrun Integer, number of runs for NMF (default: 30)
#' @param seed Integer, random seed for reproducibility (default: 123)
#' @param method Character, NMF algorithm to use (default: "brunet")
#' @return NMF object from the nmf() function
run_nmf <- function(mat, k, nrun = 30, seed = 123, method = "brunet") {
  # Ensure non-negative input matrix (NMF requirement)
  mat[mat < 0] <- 0
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Run NMF and return the NMF object
  # .options = "v" enables verbose output to monitor progress
  nmf_res <- NMF::nmf(
    mat,
    rank = k,
    nrun = nrun,
    method = method,
    .options = "v"  # verbose mode
  )
  
  return(nmf_res)
}

#' Estimate optimal NMF rank using cophenetic correlation coefficient
#'
#' The cophenetic correlation coefficient measures how well the clustering
#' structure is preserved. The optimal rank is typically where this coefficient
#' starts to decrease.
#'
#' @param mat Numeric matrix (features × samples)
#' @param max_rank Integer, maximum rank to consider (default: 8)
#' @return Integer, estimated optimal rank
estimate_rank <- function(mat, max_rank = 8) {
  # Test ranks from 2 to max_rank
  ranks <- 2:max_rank
  
  # Run NMF for each rank with 10 runs
  est <- nmf(mat, ranks, nrun = 10, method = "brunet", seed = 123)
  
  # Return rank with minimum cophenetic correlation
  # +1 compensates for starting at rank 2 (index 1 corresponds to rank 2)
  which.min(est$measures$cophenetic) + 1
}

#' Filter top variance features from a matrix
#'
#' @param mat Numeric matrix (features × samples)
#' @param n Integer, number of top features to retain
#' @return Matrix with only top n features by variance
filter_top_var <- function(mat, n = 3000) {
  # If matrix already has fewer features than requested, return as is
  if (nrow(mat) <= n) {
    return(mat)
  }
  
  # Calculate variance for each feature
  vars <- apply(mat, 1, var)
  
  # Return matrix with top n features by variance
  mat[order(vars, decreasing = TRUE)[1:n], ]
}

#' Convert data frame to matrix with row names
#'
#' @param df Data frame with first column as row identifiers
#' @return Matrix with row names set from first column
matrix_from_df <- function(df) {
  # Convert all but first column to matrix
  mat <- as.matrix(df[, -1])
  # Set row names from first column
  rownames(mat) <- df[[1]]
  mat
}

#===============================================================================
# Visualization Functions
#===============================================================================

#' Create Sankey diagram showing cluster concordance between two omics datasets
#'
#' This function creates a Sankey diagram to visualize how samples are clustered
#' in two different omics datasets and the flow of samples between clusters.
#'
#' @param data Data frame with sample clustering results for two omics datasets
#' @param v_space Character or numeric, vertical space between nodes in Sankey plot (default: "auto")
#' @param order Character, ordering of nodes (default: "ascending")
#' @param align Character, alignment of stages (default: "justify")
#' @param omics1_name Character, name of first omics dataset (default: 'Pro')
#' @param omics2_name Character, name of second omics dataset (default: 'Phos')
#' @return ggplot2 object
plotSankey = function(
    data = data,
    v_space = "auto",
    order = "ascending",
    align = "justify",
    omics1_name = 'Pro',
    omics2_name = 'Phos'
) {
  # Add weight column for Sankey plotting (each sample has weight 1)
  wide_df <- data |> mutate(weight = 1)
  
  # Convert wide format to long format required by ggsankeyfier
  sankey_long <- pivot_stages_longer(
    data = wide_df,
    stages_from = colnames(data)[-1],
    values_from = "weight"
  )
  
  # Set Sankey position parameters
  pos <- position_sankey(v_space = v_space, order = order, align = align)
  
  # Create Sankey plot
  plt <- ggplot(
    sankey_long,
    aes(
      x = stage,
      y = weight,
      group = node,
      connector = connector,
      edge_id = edge_id,
      fill = node
    )
  ) +
    # Add Sankey nodes
    geom_sankeynode(position = pos, colour = NA, alpha = 0.9) +
    # Add Sankey edges (flows)
    geom_sankeyedge(position = pos, alpha = 0.7, colour = NA) +
    # Add node labels
    geom_text(
      data = sankey_long,
      aes(label = node), 
      stat = "sankeynode",
      position = pos,
      hjust = 0.5, 
      size = 3.5, 
      fontface = "bold", 
      color = "black"
    ) +
    # Apply classic theme
    theme_classic() +
    # Add labels and title
    labs(
      title = paste0(omics1_name, " vs. ", omics2_name, " NMF cluster concordance"),
      x = NULL,
      y = "Sample count",
      fill = "Clusters"
    ) +
    # Use viridis color scale
    scale_colour_viridis_d() +
    # Customize theme elements
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold")
    )
  
  return(plt)
}

#' Create Sankey diagram with custom node ordering by group
#'
#' This function creates a Sankey diagram similar to plotSankey but with
#' custom ordering of nodes (by group prefix and numeric suffix).
#'
#' @param data Data frame with sample clustering results for two omics datasets
#' @param v_space Character or numeric, vertical space between nodes (default: "auto")
#' @param order Character, ordering of nodes (default: "as_is")
#' @param align Character, alignment of stages (default: "justify")
#' @param omics1_name Character, name of first omics dataset (default: 'Pro')
#' @param omics2_name Character, name of second omics dataset (default: 'Phos')
#' @return ggplot2 object
plotSankey_by_group = function(
    data = data,
    v_space = "auto",
    order = "as_is",
    align = "justify",
    omics1_name = 'Pro',
    omics2_name = 'Phos'
) {
  # Add weight column for Sankey plotting
  wide_df <- data |> mutate(weight = 1)
  
  # Convert to long format
  sankey_long <- pivot_stages_longer(
    data = wide_df,
    stages_from = colnames(data)[-1],
    values_from = "weight"
  )
  
  #---------------------------------------------------------------------------
  # Custom node sorting function
  #---------------------------------------------------------------------------
  #' Sort nodes by prefix then numeric value
  #' @param nodes Character vector of node names (format: "Prefix_Number")
  #' @return Sorted character vector
  custom_sort <- function(nodes) {
    # Convert to character
    nodes <- as.character(nodes)
    
    # Extract prefix (e.g., "N_Pro") and numeric part
    prefixes <- gsub("_\\d+$", "", nodes)  # Get prefix (e.g., "N_Pro")
    numbers <- as.numeric(gsub("^.*_", "", nodes))  # Get numeric part
    
    # Sort by prefix first, then by number (ascending)
    sorted_indices <- order(prefixes, numbers)
    return(nodes[sorted_indices])
  }
  
  # Apply custom sorting to nodes
  unique_nodes <- unique(sankey_long$node)
  sorted_nodes <- custom_sort(unique_nodes)
  
  # Set node factor levels for ordered display
  sankey_long$node <- factor(sankey_long$node, levels = sorted_nodes)
  
  # Set Sankey position parameters
  pos <- position_sankey(v_space = v_space, order = order, align = align)
  
  # Create Sankey plot
  plt <- ggplot(
    sankey_long,
    aes(
      x = stage,
      y = weight,
      group = node,
      connector = connector,
      edge_id = edge_id,
      fill = node
    )
  ) +
    # Add nodes and edges
    geom_sankeynode(position = pos, colour = NA, alpha = 0.9) +
    geom_sankeyedge(position = pos, alpha = 0.7, colour = NA) +
    theme_classic() +
    # Add labels and title
    labs(
      title = paste0(omics1_name, " vs. ", omics2_name, " NMF cluster concordance by groups"),
      x = NULL,
      y = "Sample count",
      fill = "Clusters"
    ) +
    # Use viridis color scale
    scale_colour_viridis_d() +
    # Add node labels
    geom_text(
      data = sankey_long,
      aes(label = node), 
      stat = "sankeynode",
      position = pos,
      hjust = 0.5, 
      size = 3.5, 
      fontface = "bold", 
      color = "black"
    ) +
    # Customize theme
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold")
    )
  
  plt
  return(plt)
}