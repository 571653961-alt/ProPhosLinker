#!/usr/bin/env Rscript

#' @title Visualize a phosphoprotein network
#' @description This function visualizes a phosphoprotein network using the igraph package.

#===============================================================================
# Functional Interaction Network Visualization
#===============================================================================
# This script performs comprehensive network visualization for protein-protein
# interaction networks with phosphoproteomics data integration, focusing on
# module-specific subnetworks identified through community detection. Key features:
#   - Module-specific subnetwork extraction
#   - Multiple filtering modes (ALL, KS, LR, TF, DEP, DEPUP, DEPDOWN, SEEDNODE)
#   - Centrality-based node filtering
#   - Phosphosite display with pie charts
#   - Disease-related node highlighting
#   - Multiple layout options
#   - Edge type color coding
#   - Functional enrichment analysis for subnetworks
#===============================================================================

suppressWarnings(library(optparse))

#===============================================================================
# COMMAND LINE ARGUMENT DEFINITION
#===============================================================================

option_list <- list(
  #---------------------------------------------------------------------------
  # Script Path
  #---------------------------------------------------------------------------
  make_option(c("--function_enrich_rscript_path"), type = "character", default = NULL,
              help = "Path to enrichment function script", metavar = "FILE"),
  
  #---------------------------------------------------------------------------
  # Input Parameters
  #---------------------------------------------------------------------------
  make_option(c("-n", "--nodes"), type = "character", default = NULL, metavar = "FILE",
              help = "Input nodes file (TSV format)"),
  make_option(c("-e", "--edges"), type = "character", default = NULL, metavar = "FILE",
              help = "Input edges file (TSV format)"),
  make_option(c("-m", "--module_info"), type = "character", default = NULL, metavar = "FILE",
              help = "Module information file (TSV format)"),
  make_option(c("-i", "--module_id"), type = "integer", default = 0, metavar = "INT",
              help = "Module ID to visualize [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Output Parameters
  #---------------------------------------------------------------------------
  make_option(c("-o", "--outdir"), type = "character", default = getwd(), metavar = "DIR",
              help = "Output directory [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Enrichment Parameters
  #---------------------------------------------------------------------------
  make_option(c("--enrich_fromType"), type = "character", default = "UNIPROT",
              help = "ID type for enrichment [default: %default], options: UNIPROT, SYMBOL", metavar = "STRING"),
  
  #---------------------------------------------------------------------------
  # Network Visualization Parameters
  #---------------------------------------------------------------------------
  make_option(c("-d", "--max_phosphoSite_displayed"), type = "integer", default = 5, metavar = "INT",
              help = "Maximum number of phosphosites to display per protein [default: %default]"),
  make_option(c("-t", "--network_mode"), type = "character", default = "ALL", metavar = "METHOD",
              help = "Subnetwork filter mode: ALL, KS, LR, TF, DEP, DEPUP, DEPDOWN, SEEDNODE [default: %default]"),
  make_option(c("-s", "--SEEDNODEID"), type = "character", default = NULL, metavar = "STRING",
              help = "Seed node ID for network expansion [default: %default]"),
  make_option(c("-f", "--node_filtering"), type = "character", default = "betweenness", metavar = "STRING",
              help = "Node filter method: degree, betweenness, closeness, harmonic, eigenvector, pagerank, alpha, hub, authority [default: %default]"),
  make_option(c("-l", "--network_layout"), type = "character", default = "kk", metavar = "STRING",
              help = "Layout method: fr, kk, dh, stress, tree, gem, graphopt, lgl, circle, grid [default: %default]"),
  make_option(c("-p", "--top_nodes_visualization_num"), type = "integer", default = 20, metavar = "INT",
              help = "Number of top nodes to visualize based on centrality"),
  
  #---------------------------------------------------------------------------
  # Omics Names
  #---------------------------------------------------------------------------
  make_option(c("--omics1_name"), type = "character", default = "Proteomics", metavar = "STR",
              help = "Name for first omics dataset [default: %default]"),
  make_option(c("--omics2_name"), type = "character", default = "Phosphoproteomics", metavar = "STR",
              help = "Name for second omics dataset [default: %default]"),
  
  #---------------------------------------------------------------------------
  # Color Parameters
  #---------------------------------------------------------------------------
  make_option(c("--node_up"), type = "character", default = "#B83A2D", metavar = "character",
              help = "Color for upregulated nodes [default: %default]"),
  make_option(c("--node_down"), type = "character", default = "#657f68", metavar = "character",
              help = "Color for downregulated nodes [default: %default]"),
  make_option(c("--node_nonsig"), type = "character", default = "#E3C79F", metavar = "character",
              help = "Color for non-significant nodes [default: %default]"),
  make_option(c("--node_notdet"), type = "character", default = "grey90", metavar = "character",
              help = "Color for not detected nodes [default: %default]"),
  make_option(c("--node_disease_border_color"), type = "character", default = "#535c54", metavar = "character",
              help = "Border color for disease-related nodes [default: %default]"),
  make_option(c("--node_notdisease_border_color"), type = "character", default = "#acabab", metavar = "character",
              help = "Border color for non-disease nodes [default: %default]"),
  make_option(c("--function_enrich_color_gradient_low"), type = "character", default = "#175663", metavar = "character",
              help = "Low end color for functional enrichment gradient [default: %default]"),
  make_option(c("--function_enrich_color_gradient_high"), type = "character", default = "#90362d", metavar = "character",
              help = "High end color for functional enrichment gradient [default: %default]"),
  make_option(c("--edge_type_colors"), type = "character",
              default = "association:#d7d6d6;physical association:#838181;binding:#838181;direct interaction:#d7d6d6;activation:#d62c0c;catalysis:#7f00ff;proximity:#FF9301;reaction:#114335;phosphorylation:#2FBE95;dephosphorylation:#6B8E23;ptmod:#8C97D6;inhibition:#0cb6d6;expression:#FCF402;regulation:#e4dadd;colocalization:#4c95cd;covalent binding:#716F74;ubiquitination:#FF4500;multiRel:#9a6728",
              metavar = "character",
              help = "Edge type colors as semicolon-separated key:value pairs [default: %default]"),
  
  #---------------------------------------------------------------------------
  # General Parameters
  #---------------------------------------------------------------------------
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE, metavar = "FLAG",
              help = "Print detailed output messages [Default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                           description = "Visualize a protein-phosphoprotein network using igraph package.")
opt <- parse_args(opt_parser)

#===============================================================================
# PARAMETER VALIDATION FUNCTION
#===============================================================================

#' Validate and clean input parameters
#'
#' @param opt List of parsed command-line arguments
#' @return Validated and cleaned parameter list
parameter_validation <- function(opt) {
  
  #---------------------------------------------------------------------------
  # Helper Functions
  #---------------------------------------------------------------------------
  
  #' Clean and normalize file paths
  clean_path <- function(path) {
    path <- gsub("^['\" ]+|['\" ]+$", "", path)
    path <- gsub("\\\\", "/", path)
    return(path)
  }
  
  #' Check if input file exists and is not empty
  check_input_file <- function(file_path, file_name) {
    if (!file.exists(file_path)) {
      stop("❌ Error: The file '", file_name, "' does not exist at the specified path: ", file_path)
    }
    if (file.size(file_path) == 0) {
      stop("❌ Error: The file '", file_name, "' is empty at the specified path: ", file_path)
    }
  }
  
  #' Create or validate output directory
  check_output_dir <- function(result_dir) {
    # Create subdirectory named after module file and module ID
    result_dir <- paste0(result_dir, "/", tools::file_path_sans_ext(basename(opt$module_info)), '_', opt$module_id, '_visualization')
    if (!dir.exists(result_dir)) {
      dir.create(result_dir, recursive = TRUE)
      print(paste("✅ Output directory created:", result_dir))
    } else {
      print(paste("✅ Output directory already exists:", result_dir))
    }
    return(result_dir)
  }
  
  #' Validate categorical parameter against allowed values
  validate_option_choices <- function(value, allowed_values, option_name) {
    if (!value %in% allowed_values) {
      stop("❌ Error: Invalid", option_name, "value '", value, "'. Please use one of the following options: ",
           paste(allowed_values, collapse = ", "))
    }
  }
  
  #' Validate numeric parameter within specified range
  validate_numeric_range <- function(value, min_val, max_val, option_name) {
    if (is.null(value)) {
      return()
    }
    if (value < min_val || value > max_val) {
      stop("❌ Error: ", option_name, " must be within the range [", min_val, ", ", max_val, "].")
    }
  }
  
  #' Parse edge colors from string format
  parse_edge_colors <- function(color_str) {
    color_str <- gsub('"', '', color_str)
    pairs <- strsplit(color_str, ";")[[1]]
    color_list <- list()
    for (pair in pairs) {
      key_val <- strsplit(pair, ":")[[1]]
      if (length(key_val) == 2) {
        color_list[[key_val[1]]] <- key_val[2]
      }
    }
    return(color_list)
  }
  
  #---------------------------------------------------------------------------
  # Path Cleaning
  #---------------------------------------------------------------------------
  opt$nodes <- clean_path(opt$nodes)
  opt$edges <- clean_path(opt$edges)
  opt$module_info <- clean_path(opt$module_info)
  opt$function_enrich_rscript_path <- clean_path(opt$function_enrich_rscript_path)
  opt$outdir <- clean_path(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Input File Validation
  #---------------------------------------------------------------------------
  check_input_file(opt$nodes, "nodes")
  check_input_file(opt$edges, "edges")
  check_input_file(opt$module_info, "module info")
  opt$outdir <- check_output_dir(opt$outdir)
  
  #---------------------------------------------------------------------------
  # Parameter Validation
  #---------------------------------------------------------------------------
  allowed_network_mode <- c("ALL", "KS", "LR", "TF", "DEP", "DEPUP", "DEPDOWN", "SEEDNODE")
  allowed_node_filtering <- c("betweenness", "degree", "closeness", "harmonic", "eigenvector", "pagerank", "alpha", "hub", "authority")
  allowed_network_layout_type <- c("fr", "kk", "dh", "stress", "tree", "gem", "graphopt", "lgl", "circle", "grid")
  allowed_enrich_fromType <- c('UNIPROT', 'SYMBOL')
  
  validate_option_choices(opt$network_mode, allowed_network_mode, "network mode")
  validate_option_choices(opt$node_filtering, allowed_node_filtering, "node filtering")
  validate_option_choices(opt$network_layout, allowed_network_layout_type, "network layout")
  validate_option_choices(opt$enrich_fromType, allowed_enrich_fromType, "enrich_fromType")
  validate_numeric_range(opt$max_phosphoSite_displayed, 1, 5, "max_phosphoSite_displayed")
  validate_numeric_range(opt$top_nodes_visualization_num, 3, 100, "top_nodes_visualization_num")
  
  opt$SEEDNODEID <- clean_path(opt$SEEDNODEID)
  opt$edge_type_colors <- parse_edge_colors(opt$edge_type_colors)
  
  if (opt$verbose) { print(opt) }
  
  return(opt)
}

# Set warning level to ignore
options(warn = -1)

# Load required libraries
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(graphlayouts)))
suppressWarnings(suppressPackageStartupMessages(library(ggforce)))
suppressWarnings(suppressPackageStartupMessages(library(scatterpie)))
suppressWarnings(suppressPackageStartupMessages(library(ggsci)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(tidygraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ggnewscale)))  # For handling multiple color scales
suppressWarnings(suppressPackageStartupMessages(library(svglite)))

# Source enrichment functions
source(opt$function_enrich_rscript_path)

# Number of top nodes to visualize (from parameters)
top_nodes_visualization_num <- opt$top_nodes_visualization_num

#===============================================================================
# SUBNETWORK CONSTRUCTION FUNCTION
#===============================================================================

#' Construct subnetwork based on module ID and filtering criteria
#'
#' @param opt List of parameters
#' @param top_nodes_visualization_num Number of top nodes to retain
#' @return igraph object of the constructed subnetwork, or error message string
subnetwork_construction <- function(opt, top_nodes_visualization_num) {
  setwd(opt$outdir)
  
  #=============================================================================
  # 1. DATA LOADING
  #=============================================================================
  f_Node_path <- opt$nodes
  f_Edge_path <- opt$edges
  f_Module_info_path <- opt$module_info
  f_Module_id <- opt$module_id
  network_mode <- opt$network_mode
  SEEDNODEID <- opt$SEEDNODEID
  node_filtering <- opt$node_filtering
  
  # Load nodes, edges, and module information
  f_Node_rawdata <- read.csv(f_Node_path, sep = '\t', stringsAsFactors = FALSE)
  f_Edge_rawdata <- read.csv(f_Edge_path, sep = '\t', stringsAsFactors = FALSE)
  f_Module_info <- read.csv(f_Module_info_path, sep = '\t', stringsAsFactors = FALSE)
  
  # Extract nodes belonging to the specified module
  module_nodesID_list <- strsplit(f_Module_info[f_Module_info$ModuleID == f_Module_id, 'ModuleItems'], ';')[[1]]
  f_Node <- f_Node_rawdata[f_Node_rawdata$NodeName %in% module_nodesID_list, ]
  
  # Validate seed node if in SEEDNODE mode
  if (network_mode == 'SEEDNODE') {
    if (length(SEEDNODEID) == 0) {
      return("❌ Error: SEEDNODEID value is empty!")
    } else if (!any(SEEDNODEID %in% module_nodesID_list)) {
      return(paste0("❌ Error: SEEDNODEID ", SEEDNODEID, " not in this module."))
    }
  }
  
  #=============================================================================
  # 2. NODE FILTERING BASED ON MODE
  #=============================================================================
  max_phosphoSite_num_displayed <- opt$max_phosphoSite_displayed
  
  # Filter nodes based on selected mode
  if (network_mode == 'DEP') {
    # Keep nodes with both protein and phosphosite differentially expressed
    filter_node_data <- f_Node[
      (f_Node$Phosphopro_class == 'UP' | f_Node$Phosphopro_class == 'DOWN') &
        (f_Node$Pro_class == 'UP' | f_Node$Pro_class == 'DOWN'), ] %>%
      group_by(NodeName) %>%
      arrange(desc(abs(Phosphopro_FC)), .by_group = TRUE) %>%
      mutate(phosphoSite_filter_index = row_number()) %>%
      filter(phosphoSite_filter_index <= max_phosphoSite_num_displayed) %>%
      ungroup()
    
  } else if (network_mode == 'DEPUP') {
    # Keep nodes with both protein and phosphosite upregulated
    filter_node_data <- f_Node[f_Node$Phosphopro_class == 'UP' & f_Node$Pro_class == 'UP', ] %>%
      group_by(NodeName) %>%
      arrange(desc(abs(Phosphopro_FC)), .by_group = TRUE) %>%
      mutate(phosphoSite_filter_index = row_number()) %>%
      filter(phosphoSite_filter_index <= max_phosphoSite_num_displayed) %>%
      ungroup()
    
  } else if (network_mode == 'DEPDOWN') {
    # Keep nodes with both protein and phosphosite downregulated
    filter_node_data <- f_Node[f_Node$Phosphopro_class == 'DOWN' & f_Node$Pro_class == 'DOWN', ] %>%
      group_by(NodeName) %>%
      arrange(desc(abs(Phosphopro_FC)), .by_group = TRUE) %>%
      mutate(phosphoSite_filter_index = row_number()) %>%
      filter(phosphoSite_filter_index <= max_phosphoSite_num_displayed) %>%
      ungroup()
    
  } else {
    # Default: keep all nodes, limit phosphosites per protein
    filter_node_data <- f_Node %>%
      group_by(NodeName) %>%
      arrange(desc(abs(Phosphopro_FC)), .by_group = TRUE) %>%
      mutate(phosphoSite_filter_index = row_number()) %>%
      filter(phosphoSite_filter_index <= max_phosphoSite_num_displayed) %>%
      ungroup()
  }
  
  # Check if any nodes remain after filtering
  if (all(is.na(filter_node_data$NodeName))) {
    return("❌ Error: This module does not contain any nodes with the specified filtering criteria.")
  }
  
  #=============================================================================
  # 3. EDGE FILTERING
  #=============================================================================
  nodes <- unique(filter_node_data$NodeName)
  edges <- f_Edge_rawdata %>%
    filter(Source %in% nodes, Target %in% nodes)
  
  #=============================================================================
  # 4. INTERACTION TYPE CLEANING
  #=============================================================================
  
  #' Clean interaction types by removing redundant information
  clean_interaction_type <- function(types) {
    result <- types
    
    # Helper function: Safe replacement
    safe_replace <- function(pattern, replacement, x) {
      new_x <- gsub(pattern, replacement, x)
      ifelse(new_x == "", x, new_x)
    }
    
    result <- safe_replace("(^|;)association(;|$)", "\\1\\2", result)
    result <- safe_replace("/association", "", result)
    
    # Handle 'physical association'
    result <- ifelse(
      grepl("physical association", result) & grepl(";", result),
      safe_replace(";physical association|physical association;", "", result),
      result
    )
    
    # Handle 'binding'
    result <- ifelse(
      grepl("binding", result) & grepl(";", result),
      safe_replace(";binding|binding;", "", result),
      result
    )
    
    # Handle 'direct interaction'
    result <- ifelse(
      grepl("direct interaction", result) & grepl(";", result),
      safe_replace(";direct interaction|direct interaction;", "", result),
      result
    )
    
    # Handle 'expression' coexisting with 'activation' or 'inhibition'
    has_expression <- grepl("expression", result)
    has_activation <- grepl("activation", result)
    has_inhibition <- grepl("inhibition", result)
    
    result <- ifelse(
      has_expression & (has_activation | has_inhibition),
      safe_replace(";expression|expression;", "", result),
      result
    )
    
    # Handle 'reaction' with 'catalysis'
    has_reaction <- grepl("(^reaction;)|(;reaction;)|(;reaction$)", result)
    has_catalysis <- grepl("catalysis", result)
    
    result <- ifelse(
      has_reaction & has_catalysis,
      safe_replace(";reaction|reaction;", "", result),
      result
    )
    
    # Clean up semicolons
    result <- safe_replace("^;|;$", "", result)
    result <- safe_replace(";;+", ";", result)
    
    # Ensure not empty
    result <- ifelse(result == "", types, result)
    
    return(result)
  }
  
  # Apply interaction type cleaning
  edges$Interaction_Type_cleaned <- clean_interaction_type(edges$Interaction_Type)
  edges$Interaction_Type_simplified <- ifelse(
    grepl(';', edges$Interaction_Type_cleaned),
    "multiRel",
    edges$Interaction_Type_cleaned
  )
  
  #=============================================================================
  # 5. INITIAL GRAPH CONSTRUCTION
  #=============================================================================
  g <- graph_from_data_frame(
    d = edges,
    directed = TRUE,
    vertices = nodes
  )
  
  # Attach edge attributes
  E(g)$Interaction_Type <- edges$Interaction_Type
  E(g)$Interaction_Type_cleaned <- edges$Interaction_Type_cleaned
  E(g)$Interaction_Type_simplified <- edges$Interaction_Type_simplified
  
  #=============================================================================
  # 6. EDGE FILTERING BY NETWORK MODE
  #=============================================================================
  
  if (network_mode == 'KS') {
    # Keep only kinase-substrate interactions
    edges_phospho <- E(g)[grepl("phospho", E(g)$Interaction_Type)]
    g <- subgraph_from_edges(
      graph = g,
      eids = edges_phospho,
      delete.vertices = TRUE
    )
    
  } else if (network_mode == 'LR') {
    # Keep only ligand-receptor interactions (binding)
    edges_LR <- E(g)[grepl("binding", E(g)$Interaction_Type)]
    E(g)$Interaction_Type_simplified <- ifelse(E(g)$Interaction_Type_simplified == 'binding', 'binding', 'multiRel')
    g <- subgraph_from_edges(
      graph = g,
      eids = edges_LR,
      delete.vertices = TRUE
    )
    
  } else if (network_mode == 'TF') {
    # Keep only transcription factor related interactions
    edges_TF <- E(g)[
      grepl(
        pattern = "activation|regulation|repression|inhibition",
        x = E(g)$Interaction_Type,
        ignore.case = TRUE
      )
    ]
    g <- subgraph_from_edges(
      graph = g,
      eids = edges_TF,
      delete.vertices = TRUE
    )
    
  } else if (network_mode == 'SEEDNODE') {
    # Expand from seed node to include neighbors
    focus_node_edges <- E(g)[.inc(V(g)[name == SEEDNODEID])]
    connected_vertices <- ends(g, focus_node_edges, names = TRUE)
    connected_vertices <- unique(as.vector((connected_vertices)))
    
    # Include second-degree neighbors if needed
    if (length(connected_vertices) >= 10) {
      final_node <- connected_vertices
    } else {
      final_node <- NULL
      for (i in connected_vertices) {
        node_edges <- E(g)[.inc(V(g)[name == i])]
        vertices <- ends(g, node_edges, names = TRUE)
      }
      final_node <- union(final_node, vertices)
    }
    
    g <- induced_subgraph(
      graph = g,
      vids = final_node
    )
  }
  
  # Check if graph is too small
  if (length(V(g)) == 1) {
    return("❌ Error: This module only contains one node with the specified filtering criteria.")
  }
  
  #=============================================================================
  # 7. CENTRALITY-BASED NODE FILTERING
  #=============================================================================
  
  # Calculate selected centrality measure
  if (node_filtering == 'betweenness') {
    values <- betweenness(g, directed = FALSE, normalized = FALSE)
  } else if (node_filtering == 'degree') {
    values <- degree(g, mode = "all", normalized = FALSE)
  } else if (node_filtering == 'closeness') {
    values <- closeness(g, mode = "all", normalized = FALSE)
  } else if (node_filtering == 'harmonic') {
    values <- harmonic_centrality(g, mode = "all", normalized = FALSE)
  } else if (node_filtering == 'eigenvector') {
    values <- evcent(g, directed = FALSE)$vector
  } else if (node_filtering == 'pagerank') {
    values <- page_rank(g, directed = FALSE)$vector
  } else if (node_filtering == 'alpha') {
    values <- alpha_centrality(g, alpha = 0.85)
  } else if (node_filtering == 'hub') {
    values <- hub_score(g)$vector
  } else if (node_filtering == 'authority') {
    values <- authority_score(g)$vector
  }
  
  # Select top nodes based on centrality
  top_nodes_index <- order(values, decreasing = TRUE)[1:min(top_nodes_visualization_num, vcount(g))]
  top_nodes_ids <- V(g)$name[top_nodes_index]
  
  if (all(is.na(top_nodes_ids))) {
    return("❌ Error: This module does not contain any nodes with the specified filtering criteria.")
  }
  
  # Create subgraph with top nodes
  topnodes_subgraph <- induced_subgraph(g, vids = top_nodes_ids)
  
  # Remove isolated nodes for non-DEP modes
  if (!(network_mode %in% c("DEP", "DEPUP", "DEPDOWN"))) {
    topnodes_subgraph <- delete_vertices(
      topnodes_subgraph,
      v = which(degree(topnodes_subgraph, mode = "all") == 0)
    )
  }
  
  #=============================================================================
  # 8. ADD NODE ATTRIBUTES FOR PIE CHARTS
  #=============================================================================
  
  # Initialize vertex attributes
  vertex_count <- vcount(topnodes_subgraph)
  
  # Initialize disease related attributes
  V(topnodes_subgraph)$DiseaseRelated <- rep("no", vertex_count)
  
  # Initialize all attribute combinations for phosphosites
  prefixes <- c("", paste0("PhosphoPro_site", LETTERS[1:5], '_'))
  suffixes <- c("UP", "DOWN", "Non_significant", "Not_Detected")
  
  for (prefix in prefixes) {
    for (suffix in suffixes) {
      attr_name <- paste0(prefix, suffix)
      V(topnodes_subgraph)$attr_name <- rep(0, vertex_count)
      names(vertex_attr(topnodes_subgraph))[length(vertex_attr(topnodes_subgraph))] <- attr_name
    }
  }
  
  # Assign values to attributes based on filtered node data
  for (nodeName_index in 1:length(V(topnodes_subgraph)$name)) {
    nodeName <- V(topnodes_subgraph)$name[nodeName_index]
    node_phosphoprot_info <- filter_node_data[filter_node_data$NodeName == nodeName, ]
    node_protein_class <- filter_node_data[filter_node_data$NodeName == nodeName, ]$Pro_class[1]
    Node_DiseaseRelated <- filter_node_data[filter_node_data$NodeName == nodeName, ]$DiseaseRelated[1]
    V(topnodes_subgraph)$DiseaseRelated[nodeName_index] <- Node_DiseaseRelated
    
    # Set phosphosite attributes for each filtered site
    for (phosphoSite_filter_index in node_phosphoprot_info$phosphoSite_filter_index) {
      Phosphopro_class <- node_phosphoprot_info[node_phosphoprot_info$phosphoSite_filter_index == phosphoSite_filter_index, ]$Phosphopro_class
      
      if (is.na(Phosphopro_class) || Phosphopro_class == "") {
        Phosphopro_class <- "None"
      }
      
      # Set attributes for sites A through E
      if (phosphoSite_filter_index == 1) {
        if (Phosphopro_class == 'UP') {
          V(topnodes_subgraph)$PhosphoPro_siteA_UP[nodeName_index] <- 1
        } else if (Phosphopro_class == 'DOWN') {
          V(topnodes_subgraph)$PhosphoPro_siteA_DOWN[nodeName_index] <- 1
        } else if (Phosphopro_class == 'Non-significant') {
          V(topnodes_subgraph)$PhosphoPro_siteA_Non_significant[nodeName_index] <- 1
        } else if (Phosphopro_class == 'None') {
          V(topnodes_subgraph)$PhosphoPro_siteA_Not_Detected[nodeName_index] <- 1
        }
      }
      if (phosphoSite_filter_index == 2) {
        if (Phosphopro_class == 'UP') {
          V(topnodes_subgraph)$PhosphoPro_siteB_UP[nodeName_index] <- 1
        } else if (Phosphopro_class == 'DOWN') {
          V(topnodes_subgraph)$PhosphoPro_siteB_DOWN[nodeName_index] <- 1
        } else if (Phosphopro_class == 'Non-significant') {
          V(topnodes_subgraph)$PhosphoPro_siteB_Non_significant[nodeName_index] <- 1
        } else if (Phosphopro_class == 'None') {
          V(topnodes_subgraph)$PhosphoPro_siteB_Not_Detected[nodeName_index] <- 1
        }
      }
      if (phosphoSite_filter_index == 3) {
        if (Phosphopro_class == 'UP') {
          V(topnodes_subgraph)$PhosphoPro_siteC_UP[nodeName_index] <- 1
        } else if (Phosphopro_class == 'DOWN') {
          V(topnodes_subgraph)$PhosphoPro_siteC_DOWN[nodeName_index] <- 1
        } else if (Phosphopro_class == 'Non-significant') {
          V(topnodes_subgraph)$PhosphoPro_siteC_Non_significant[nodeName_index] <- 1
        } else if (Phosphopro_class == 'None') {
          V(topnodes_subgraph)$PhosphoPro_siteC_Not_Detected[nodeName_index] <- 1
        }
      }
      if (phosphoSite_filter_index == 4) {
        if (Phosphopro_class == 'UP') {
          V(topnodes_subgraph)$PhosphoPro_siteD_UP[nodeName_index] <- 1
        } else if (Phosphopro_class == 'DOWN') {
          V(topnodes_subgraph)$PhosphoPro_siteD_DOWN[nodeName_index] <- 1
        } else if (Phosphopro_class == 'Non-significant') {
          V(topnodes_subgraph)$PhosphoPro_siteD_Non_significant[nodeName_index] <- 1
        } else if (Phosphopro_class == 'None') {
          V(topnodes_subgraph)$PhosphoPro_siteD_Not_Detected[nodeName_index] <- 1
        }
      }
      if (phosphoSite_filter_index == 5) {
        if (Phosphopro_class == 'UP') {
          V(topnodes_subgraph)$PhosphoPro_siteE_UP[nodeName_index] <- 1
        } else if (Phosphopro_class == 'DOWN') {
          V(topnodes_subgraph)$PhosphoPro_siteE_DOWN[nodeName_index] <- 1
        } else if (Phosphopro_class == 'Non-significant') {
          V(topnodes_subgraph)$PhosphoPro_siteE_Non_significant[nodeName_index] <- 1
        } else if (Phosphopro_class == 'None') {
          V(topnodes_subgraph)$PhosphoPro_siteE_Not_Detected[nodeName_index] <- 1
        }
      }
    }
    
    # Calculate summary counts based on protein class
    total_count <- V(topnodes_subgraph)$PhosphoPro_siteA_UP[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteB_UP[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteC_UP[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteD_UP[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteE_UP[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteA_DOWN[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteB_DOWN[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteC_DOWN[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteD_DOWN[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteE_DOWN[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteA_Non_significant[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteB_Non_significant[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteC_Non_significant[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteD_Non_significant[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteE_Non_significant[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteA_Not_Detected[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteB_Not_Detected[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteC_Not_Detected[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteD_Not_Detected[nodeName_index] +
      V(topnodes_subgraph)$PhosphoPro_siteE_Not_Detected[nodeName_index]
    
    # Assign to appropriate summary column based on protein class
    if (node_protein_class == 'UP') {
      V(topnodes_subgraph)$UP[nodeName_index] <- total_count
    } else if (node_protein_class == 'DOWN') {
      V(topnodes_subgraph)$DOWN[nodeName_index] <- total_count
    } else if (node_protein_class == 'Non-significant') {
      V(topnodes_subgraph)$Non_significant[nodeName_index] <- total_count
    } else if (node_protein_class == 'None') {
      V(topnodes_subgraph)$Not_Detected[nodeName_index] <- total_count
    }
  }
  
  return(topnodes_subgraph)
}

#===============================================================================
# NETWORK VISUALIZATION FUNCTION
#===============================================================================

#' Visualize the constructed network
#'
#' @param g igraph object to visualize
#' @param opt List of parameters
#' @return Success message
network_visualiztion <- function(g, opt) {
  layout_type <- opt$network_layout
  network_mode <- opt$network_mode
  SEEDNODEID <- opt$SEEDNODEID
  
  #=============================================================================
  # 1. COMPUTE LAYOUT COORDINATES
  #=============================================================================
  
  if (layout_type == 'fr') {
    xy <- layout_with_fr(g)
  } else if (layout_type == 'kk') {
    xy <- layout_with_kk(g)
  } else if (layout_type == 'dh') {
    xy <- layout_with_dh(g)
  } else if (layout_type == 'stress') {
    xy <- layout_with_stress(g)
  } else if (layout_type == 'gem') {
    xy <- layout_with_gem(g)
  } else if (layout_type == 'graphopt') {
    xy <- layout_with_graphopt(g)
  } else if (layout_type == 'tree') {
    xy <- layout_as_tree(g)
  } else if (layout_type == 'lgl') {
    xy <- layout_with_lgl(g)
  } else if (layout_type == 'circle') {
    xy <- layout_in_circle(g)
  } else if (layout_type == 'grid') {
    xy <- layout_on_grid(g)
  }
  
  # Assign layout coordinates
  V(g)$x <- xy[, 1]
  V(g)$y <- xy[, 2]
  
  #=============================================================================
  # 2. EXTRACT NODE AND EDGE DATA
  #=============================================================================
  nodes <- as_data_frame(g, what = "vertices")
  edges <- as_data_frame(g, what = "edges")
  
  # Add border styling based on disease relevance
  nodes <- nodes %>%
    mutate(
      node_border_color = ifelse(DiseaseRelated == "yes", opt$node_disease_border_color, opt$node_notdisease_border_color),
      node_border_linetype = ifelse(DiseaseRelated == "yes", "longdash", "solid")
    )
  
  # Add coordinate information to edges
  edges <- edges %>%
    dplyr::left_join(nodes %>% dplyr::select(name, x, y), by = c("from" = "name")) %>%
    dplyr::rename(x_start = x, y_start = y) %>%
    dplyr::left_join(nodes %>% dplyr::select(name, x, y), by = c("to" = "name")) %>%
    dplyr::rename(x_end = x, y_end = y)
  
  #=============================================================================
  # 3. DEFINE COLOR MAPPINGS
  #=============================================================================
  
  # Pie chart column names
  pie_cols <- c(
    "PhosphoPro_siteA_UP", "PhosphoPro_siteA_DOWN", "PhosphoPro_siteA_Non_significant", "PhosphoPro_siteA_Not_Detected",
    "PhosphoPro_siteB_UP", "PhosphoPro_siteB_DOWN", "PhosphoPro_siteB_Non_significant", "PhosphoPro_siteB_Not_Detected",
    "PhosphoPro_siteC_UP", "PhosphoPro_siteC_DOWN", "PhosphoPro_siteC_Non_significant", "PhosphoPro_siteC_Not_Detected",
    "PhosphoPro_siteD_UP", "PhosphoPro_siteD_DOWN", "PhosphoPro_siteD_Non_significant", "PhosphoPro_siteD_Not_Detected",
    "PhosphoPro_siteE_UP", "PhosphoPro_siteE_DOWN", "PhosphoPro_siteE_Non_significant", "PhosphoPro_siteE_Not_Detected",
    "UP", "DOWN", "Non_significant", "Not_Detected"
  )
  
  # Node fill colors
  fill_values <- c(
    "UP" = opt$node_up,
    "DOWN" = opt$node_down,
    "Non_significant" = opt$node_nonsig,
    "Not_Detected" = opt$node_notdet
  )
  
  # Add colors for individual phosphosites
  for (site in c("E", "D", "C", "B", "A")) {
    fill_values[paste0("PhosphoPro_site", site, "_UP")] <- opt$node_up
    fill_values[paste0("PhosphoPro_site", site, "_DOWN")] <- opt$node_down
    fill_values[paste0("PhosphoPro_site", site, "_Non_significant")] <- opt$node_nonsig
    fill_values[paste0("PhosphoPro_site", site, "_Not_Detected")] <- opt$node_notdet
  }
  
  # Edge type colors from parameters
  edge_type_colors <- opt$edge_type_colors
  
  # Edge line types
  edge_type_linetypes <- c(
    "association" = "dotted",
    "physical association" = "dotted",
    "binding" = "dotted",
    "direct interaction" = "dotted",
    "activation" = "solid",
    "catalysis" = "solid",
    "proximity" = "solid",
    "reaction" = "solid",
    "phosphorylation" = "solid",
    "ptmod" = "solid",
    "inhibition" = "solid",
    "expression" = "solid",
    "colocalization" = "solid",
    "covalent binding" = "solid",
    "dephosphorylation" = "solid",
    "ubiquitination" = "solid",
    "multiRel" = "solid",
    "regulation" = "solid"
  )
  
  #=============================================================================
  # 4. CALCULATE EDGE ADJUSTMENTS FOR VISUALIZATION
  #=============================================================================
  
  # Calculate shortening ratio based on layout type
  shorten_ratio <- 0.1 * length(rownames(nodes)) / 10
  loop_offset <- 0.1 * length(rownames(nodes)) / 12
  
  # Adjust parameters for different layouts
  if (layout_type == "gem") {
    shorten_ratio <- 0.8 * length(rownames(nodes)) * 2
    loop_offset <- 0.1 * length(rownames(nodes)) / (top_nodes_visualization_num + 5) * 300
  } else if (layout_type == "graphopt") {
    shorten_ratio <- 0.2 * length(rownames(nodes)) / (top_nodes_visualization_num) * 40
    loop_offset <- 0.1 * length(rownames(nodes)) / (top_nodes_visualization_num + 5) * 120
  } else if (layout_type == "tree") {
    shorten_ratio <- 0.1 * length(rownames(nodes)) / 5
  } else if (layout_type == "circle") {
    shorten_ratio <- 0.1 * length(rownames(nodes)) / (top_nodes_visualization_num + 1)
    loop_offset <- 0.1 * length(rownames(nodes)) / (top_nodes_visualization_num + 5)
  } else if (layout_type == "dh") {
    shorten_ratio <- 0.1 * length(rownames(nodes)) / (top_nodes_visualization_num) * 10
  } else if (layout_type == "lgl") {
    shorten_ratio <- 0.1 * length(rownames(nodes)) / (top_nodes_visualization_num) * 8
    loop_offset <- 0.1 * length(rownames(nodes)) / (top_nodes_visualization_num) * 6
  }
  
  # Identify self-loop edges
  edges <- edges %>%
    mutate(is_self_loop = (from == to))
  
  # Calculate adjusted coordinates
  edges <- edges %>%
    mutate(
      dx = x_end - x_start,
      dy = y_end - y_start,
      length = sqrt(dx^2 + dy^2),
      x_start_adj = ifelse(!is_self_loop, x_start + dx * shorten_ratio / length, x_start),
      y_start_adj = ifelse(!is_self_loop, y_start + dy * shorten_ratio / length, y_start),
      x_end_adj = ifelse(!is_self_loop, x_end - dx * shorten_ratio / length, x_end),
      y_end_adj = ifelse(!is_self_loop, y_end - dy * shorten_ratio / length, y_end)
    ) %>%
    mutate(
      x_start_adj = ifelse(is_self_loop, x_start - loop_offset, x_start_adj),
      y_start_adj = ifelse(is_self_loop, y_start + loop_offset * 0.5, y_start_adj),
      x_end_adj = ifelse(is_self_loop, x_start + loop_offset, x_end_adj),
      y_end_adj = ifelse(is_self_loop, y_start + loop_offset * 0.5, y_end_adj)
    )
  
  # Separate self-loops from regular edges
  self_loops_all <- edges %>% filter(is_self_loop)
  non_self_edges_all <- edges %>% filter(!is_self_loop)
  
  # Add small epsilon to ensure all pie slices are visible
  epsilon <- 0.00001
  nodes$UP <- nodes$UP + epsilon
  nodes$DOWN <- nodes$DOWN + epsilon
  nodes$Non_significant <- nodes$Non_significant + epsilon
  nodes$Not_Detected <- nodes$Not_Detected + epsilon
  
  #=============================================================================
  # 5. GENERATE PLOTS FOR DIFFERENT MODES
  #=============================================================================
  
  if (network_mode == 'ALL') {
    # Generate plots for all mode types
    for (network_mode in c('ALL', 'LR', 'TF', 'KS')) {
      
      # Extract network type name from output directory
      nettype <- strsplit(opt$outdir, "/")[[1]]
      nettype_name <- strsplit(nettype[length(nettype) - 1], '_')[[1]][1]
      if (nettype_name == 'ALL') {
        nettype_name <- ''
      }
      
      # Set title based on mode
      if (network_mode == 'KS') {
        title_str <- paste0("Kinase-Substrate Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
      } else if (network_mode == 'LR') {
        title_str <- paste0("Ligand-Receptor Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
      } else if (network_mode == 'TF') {
        title_str <- paste0("Transcriptional Regulation Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
      } else if (network_mode == 'DEP') {
        title_str <- paste0("Differential Protein and Phosphosite Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
      } else if (network_mode == 'DEPUP') {
        title_str <- paste0("Upregulated Protein and Phosphosite Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
      } else if (network_mode == 'DEPDOWN') {
        title_str <- paste0("Downregulated Protein and Phosphosite Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
      } else {
        title_str <- paste0("Comprehensive Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
      }
      
      non_self_edges <- non_self_edges_all
      self_loops <- self_loops_all
      
      # Filter edges based on mode
      if (network_mode == 'KS') {
        mode_edges <- E(g)[grepl("phospho", E(g)$Interaction_Type)]
        sub_g <- subgraph_from_edges(
          graph = g,
          eids = mode_edges,
          delete.vertices = TRUE
        )
        non_self_edges$alpha_mask <- ifelse(grepl(pattern = "phospho", x = non_self_edges_all$Interaction_Type, ignore.case = TRUE), 0.5, 0)
        self_loops$alpha_mask <- ifelse(grepl(pattern = "phospho", x = self_loops_all$Interaction_Type, ignore.case = TRUE), 0.5, 0)
        
        drop_col <- c("alpha_mask", "color_mask")
        filter_edges <- rbind(non_self_edges[non_self_edges$alpha_mask != 0, ], self_loops[self_loops$alpha_mask != 0, ])
        
        # Save filtered data
        write_tsv(nodes[nodes$name %in% V(sub_g)$name, ][, !(colnames(nodes) %in% drop_col)],
                  file = file.path(opt$outdir, paste0("nodes_", network_mode, ".tsv")))
        write_tsv(filter_edges[, !(colnames(filter_edges) %in% drop_col)],
                  file = file.path(opt$outdir, paste0("edges_", network_mode, ".tsv")))
        
      } else if (network_mode == 'LR') {
        mode_edges <- E(g)[grepl("binding", E(g)$Interaction_Type)]
        sub_g <- subgraph_from_edges(
          graph = g,
          eids = mode_edges,
          delete.vertices = TRUE
        )
        non_self_edges$alpha_mask <- ifelse(grepl(pattern = "binding", x = non_self_edges_all$Interaction_Type, ignore.case = TRUE), 0.5, 0)
        self_loops$alpha_mask <- ifelse(grepl(pattern = "binding", x = self_loops_all$Interaction_Type, ignore.case = TRUE), 0.5, 0)
        
        non_self_edges$Interaction_Type_simplified <- ifelse(grepl(pattern = "binding", x = non_self_edges_all$Interaction_Type),
                                                             ifelse(non_self_edges$Interaction_Type_simplified == 'binding', 'binding', 'multiRel'),
                                                             non_self_edges$Interaction_Type_simplified)
        self_loops$Interaction_Type_simplified <- ifelse(grepl(pattern = "binding", x = self_loops$Interaction_Type),
                                                         ifelse(self_loops$Interaction_Type_simplified == 'binding', 'binding', 'multiRel'),
                                                         self_loops$Interaction_Type_simplified)
        
        drop_col <- c("alpha_mask", "color_mask")
        filter_edges <- rbind(non_self_edges[non_self_edges$alpha_mask != 0, ], self_loops[self_loops$alpha_mask != 0, ])
        
        write_tsv(nodes[nodes$name %in% V(sub_g)$name, ][, !(colnames(nodes) %in% drop_col)],
                  file = file.path(opt$outdir, paste0("nodes_", network_mode, ".tsv")))
        write_tsv(filter_edges[, !(colnames(filter_edges) %in% drop_col)],
                  file = file.path(opt$outdir, paste0("edges_", network_mode, ".tsv")))
        
      } else if (network_mode == 'TF') {
        mode_edges <- E(g)[
          grepl(
            pattern = "activation|regulation|repression|inhibition",
            x = E(g)$Interaction_Type,
            ignore.case = TRUE
          )
        ]
        sub_g <- subgraph_from_edges(
          graph = g,
          eids = mode_edges,
          delete.vertices = TRUE
        )
        non_self_edges$alpha_mask <- ifelse(grepl(pattern = "activation|regulation|repression|inhibition", x = non_self_edges_all$Interaction_Type, ignore.case = TRUE), 0.5, 0)
        self_loops$alpha_mask <- ifelse(grepl(pattern = "activation|regulation|repression|inhibition", x = self_loops_all$Interaction_Type, ignore.case = TRUE), 0.5, 0)
        
        drop_col <- c("alpha_mask", "color_mask")
        filter_edges <- rbind(non_self_edges[non_self_edges$alpha_mask != 0, ], self_loops[self_loops$alpha_mask != 0, ])
        
        write_tsv(nodes[nodes$name %in% V(sub_g)$name, ][, !(colnames(nodes) %in% drop_col)],
                  file = file.path(opt$outdir, paste0("nodes_", network_mode, ".tsv")))
        write_tsv(filter_edges[, !(colnames(filter_edges) %in% drop_col)],
                  file = file.path(opt$outdir, paste0("edges_", network_mode, ".tsv")))
        
      } else {
        non_self_edges <- non_self_edges_all
        self_loops <- self_loops_all
        sub_g <- g
        write_tsv(nodes, file = file.path(opt$outdir, paste0("nodes_", network_mode, ".tsv")))
        write_tsv(edges, file = file.path(opt$outdir, paste0("edges_", network_mode, ".tsv")))
      }
      
      mode_nodes <- V(sub_g)$name
      
      if (length(mode_nodes) > 0) {
        # Set node visibility masks
        nodes$alpha_mask <- rep(0, nrow(nodes))
        nodes$alpha_mask <- ifelse(nodes$name %in% mode_nodes, 1, 0)
        nodes$color_mask <- rep(NA, nrow(nodes))
        nodes$color_mask <- ifelse(nodes$name %in% mode_nodes, nodes$node_border_color, NA)
        
        # Create plot based on mode
        if (network_mode == 'ALL') {
          # Full plot with all elements visible
          sub_p <- ggplot() +
            geom_segment(
              data = non_self_edges,
              aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
                  color = Interaction_Type_simplified, linetype = Interaction_Type_simplified),
              alpha = 0.5, linewidth = 0.9,
              arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
            ) +
            geom_curve(
              data = self_loops,
              aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
                  color = Interaction_Type_simplified, linetype = Interaction_Type_simplified),
              curvature = -1, angle = 90, alpha = 0.5, linewidth = 0.9,
              arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
            ) +
            scale_color_manual(name = "Edge Type", values = edge_type_colors, na.value = "#d7d6d6") +
            scale_linetype_manual(name = "Edge Type", values = edge_type_linetypes) +
            ggnewscale::new_scale_color() +
            ggnewscale::new_scale("linetype") +
            scatterpie::geom_scatterpie(
              data = nodes,
              aes(x = x, y = y, group = name, colour = node_border_color, linetype = node_border_linetype),
              cols = pie_cols, pie_scale = 1.1, lwd = 0.85
            ) +
            ggrepel::geom_text_repel(
              data = nodes,
              aes(x = x, y = y, label = name),
              size = 4, color = "grey30", family = "sans", fontface = "bold",
              box.padding = 1, max.overlaps = Inf, min.segment.length = 0,
              segment.color = "grey50", segment.size = 0.5, force = 2
            ) +
            scale_size_identity() +
            scale_fill_manual(
              name = "Protein/PhosphoProtein Expression",
              values = fill_values,
              breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
              drop = FALSE
            ) +
            scale_linetype_identity(
              name = "Node Border", guide = "legend",
              breaks = c("longdash", "solid"),
              labels = c("Disease-related", "Disease-Notrelated")
            ) +
            scale_colour_identity(
              name = "Node Border", guide = "legend",
              breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
              labels = c("Disease-related", "Disease-Notrelated")
            ) +
            coord_fixed() +
            theme_void() +
            labs(title = title_str) +
            theme(
              plot.title = element_text(
                size = 20, face = "bold", family = "sans",
                hjust = 0, vjust = 1,
                margin = margin(b = 15, t = 10), lineheight = 0.8
              ),
              plot.title.position = "plot",
              plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
              legend.position = "right",
              text = element_text(family = "sans"),
              panel.background = element_blank(),
              plot.background = element_rect(fill = "white", color = NA),
              legend.background = element_rect(fill = "white", color = NA),
              legend.key = element_blank()
            )
        } else if (network_mode == 'SEEDNODE') {
          # SEEDNODE mode plot with seed node highlighting
          sub_p <- ggplot() +
            geom_segment(
              data = non_self_edges,
              aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
                  color = Interaction_Type_simplified, linetype = Interaction_Type_simplified),
              alpha = 0.5, linewidth = 0.9,
              arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
            ) +
            geom_curve(
              data = self_loops,
              aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
                  color = Interaction_Type_simplified, linetype = Interaction_Type_simplified),
              curvature = -1, angle = 90, alpha = 0.5, linewidth = 0.9,
              arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
            ) +
            scale_color_manual(name = "Edge Type", values = edge_type_colors, na.value = "grey80") +
            scale_linetype_manual(name = "Edge Type", values = edge_type_linetypes) +
            ggnewscale::new_scale_color() +
            ggnewscale::new_scale("linetype") +
            scatterpie::geom_scatterpie(
              data = nodes,
              aes(x = x, y = y, group = name, colour = node_border_color, linetype = node_border_linetype),
              cols = pie_cols, pie_scale = 1.1, lwd = 0.85
            ) +
            ggrepel::geom_text_repel(
              data = nodes,
              aes(x = x, y = y, label = name, size = ifelse(name == SEEDNODEID, 6, 4)),
              color = "grey30", family = "sans", fontface = "bold",
              box.padding = 1, max.overlaps = Inf, min.segment.length = 0,
              segment.color = "grey50", segment.size = 0.5, force = 2
            ) +
            scale_size_identity() +
            scale_fill_manual(
              name = "Protein/PhosphoProtein Expression",
              values = fill_values,
              breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
              drop = FALSE
            ) +
            scale_linetype_identity(
              name = "Node Border", guide = "legend",
              breaks = c("longdash", "solid"),
              labels = c("Disease-related", "Disease-Notrelated")
            ) +
            scale_colour_identity(
              name = "Node Border", guide = "legend",
              breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
              labels = c("Disease-related", "Disease-Notrelated")
            ) +
            coord_fixed() +
            theme_void() +
            labs(title = title_str) +
            theme(
              plot.title = element_text(
                size = 20, face = "bold", family = "sans",
                hjust = 0, vjust = 1,
                margin = margin(b = 15, t = 10), lineheight = 0.8
              ),
              plot.title.position = "plot",
              plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
              legend.position = "right",
              text = element_text(family = "sans"),
              panel.background = element_blank(),
              plot.background = element_rect(fill = "white", color = NA),
              legend.background = element_rect(fill = "white", color = NA),
              legend.key = element_blank()
            )
        } else {
          # Filtered mode plot with alpha masking
          sub_p <- ggplot() +
            geom_segment(
              data = non_self_edges,
              aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
                  color = Interaction_Type_simplified, alpha = alpha_mask),
              linewidth = 0.9,
              arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
            ) +
            geom_curve(
              data = self_loops,
              aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
                  color = Interaction_Type_simplified, alpha = alpha_mask),
              curvature = -1, angle = 90, linewidth = 0.9,
              arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
            ) +
            scale_color_manual(name = "Edge Type", values = edge_type_colors, na.value = "grey80") +
            ggnewscale::new_scale_color() +
            ggnewscale::new_scale("linetype") +
            scatterpie::geom_scatterpie(
              data = nodes,
              aes(x = x, y = y, group = name, colour = color_mask,
                  linetype = node_border_linetype, alpha = alpha_mask),
              cols = pie_cols, pie_scale = 1.1, lwd = 0.85
            ) +
            ggrepel::geom_text_repel(
              data = nodes,
              aes(x = x, y = y, label = name, alpha = alpha_mask),
              size = 4, color = "grey30", family = "sans", fontface = "bold",
              box.padding = 1, max.overlaps = Inf, min.segment.length = 0,
              segment.color = "grey50", segment.size = 0.5, force = 2
            ) +
            scale_fill_manual(
              name = "Protein/PhosphoProtein Expression",
              values = fill_values,
              breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
              drop = FALSE
            ) +
            scale_linetype_identity(
              name = "Node Border", guide = "legend",
              breaks = c("longdash", "solid"),
              labels = c("Disease-related", "Disease-Notrelated")
            ) +
            scale_colour_identity(
              name = "Node Border", guide = "legend",
              breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
              labels = c("Disease-related", "Disease-Notrelated")
            ) +
            coord_fixed() +
            theme_void() +
            labs(title = title_str) +
            theme(
              plot.title = element_text(
                size = 20, face = "bold", family = "sans",
                hjust = 0, vjust = 1,
                margin = margin(b = 15, t = 10), lineheight = 0.8
              ),
              plot.title.position = "plot",
              plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
              legend.position = "right",
              text = element_text(family = "sans"),
              panel.background = element_blank(),
              plot.background = element_rect(fill = "white", color = NA),
              legend.background = element_rect(fill = "white", color = NA),
              legend.key = element_blank()
            ) +
            guides(alpha = "none")
        }
        
        # Save plot with appropriate dimensions
        if (layout_type == 'tree') {
          ggsave(filename = paste0("pro_phosphopro_network_", network_mode, ".png"),
                 plot = sub_p, width = 16, height = 12, units = "in", dpi = 300)
        } else {
          ggsave(filename = paste0("pro_phosphopro_network_", network_mode, ".png"),
                 plot = sub_p, width = 15, height = 12, units = "in", dpi = 300)
        }
      }
    }
  } else {
    # Single mode plot (non-ALL)
    non_self_edges <- non_self_edges_all
    self_loops <- self_loops_all
    
    # Save node and edge data
    write_tsv(nodes, file = file.path(opt$outdir, paste0("nodes_", network_mode, ".tsv")))
    write_tsv(edges, file = file.path(opt$outdir, paste0("edges_", network_mode, ".tsv")))
    
    # Extract network type name for title
    nettype <- strsplit(opt$outdir, "/")[[1]]
    nettype_name <- strsplit(nettype[length(nettype) - 1], '_')[[1]][1]
    if (nettype_name == 'ALL') {
      nettype_name <- ''
    }
    
    # Set title based on mode
    if (network_mode == 'KS') {
      title_str <- paste0("Kinase-Substrate Visualization Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
    } else if (network_mode == 'LR') {
      title_str <- paste0("Ligand-Receptor Visualization Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
    } else if (network_mode == 'TF') {
      title_str <- paste0("Transcriptional Regulation Visualization Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
    } else if (network_mode == 'DEP') {
      title_str <- paste0("Differential Protein and Phosphosite Visualization Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
    } else if (network_mode == 'DEPUP') {
      title_str <- paste0("Upregulated Protein and Phosphosite Visualization Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
    } else if (network_mode == 'DEPDOWN') {
      title_str <- paste0("Downregulated Protein and Phosphosite Visualization Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
    } else {
      title_str <- paste0("Comprehensive Visualization Subnetwork of Proteomics and Phosphoproteomics ", nettype_name, " Networks")
    }
    
    # Create plot based on mode
    if (network_mode == 'SEEDNODE') {
      sub_p <- ggplot() +
        geom_segment(
          data = non_self_edges,
          aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
              color = Interaction_Type_simplified, linetype = Interaction_Type_simplified),
          alpha = 0.5, linewidth = 0.9,
          arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
        ) +
        geom_curve(
          data = self_loops,
          aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
              color = Interaction_Type_simplified, linetype = Interaction_Type_simplified),
          curvature = -1, angle = 90, alpha = 0.5, linewidth = 0.9,
          arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
        ) +
        scale_color_manual(name = "Edge Type", values = edge_type_colors, na.value = "grey80") +
        scale_linetype_manual(name = "Edge Type", values = edge_type_linetypes) +
        ggnewscale::new_scale_color() +
        ggnewscale::new_scale("linetype") +
        scatterpie::geom_scatterpie(
          data = nodes,
          aes(x = x, y = y, group = name, colour = node_border_color, linetype = node_border_linetype),
          cols = pie_cols, pie_scale = 1.1, lwd = 0.85
        ) +
        ggrepel::geom_text_repel(
          data = nodes,
          aes(x = x, y = y, label = name, size = ifelse(name == SEEDNODEID, 6, 4)),
          color = "grey30", family = "sans", fontface = "bold",
          box.padding = 1, max.overlaps = Inf, min.segment.length = 0,
          segment.color = "grey50", segment.size = 0.5, force = 2
        ) +
        scale_size_identity() +
        scale_fill_manual(
          name = "Protein/PhosphoProtein Expression",
          values = fill_values,
          breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
          drop = FALSE
        ) +
        scale_linetype_identity(
          name = "Node Border", guide = "legend",
          breaks = c("longdash", "solid"),
          labels = c("Disease-related", "Disease-Notrelated")
        ) +
        scale_colour_identity(
          name = "Node Border", guide = "legend",
          breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
          labels = c("Disease-related", "Disease-Notrelated")
        ) +
        coord_fixed() +
        theme_void() +
        labs(title = title_str) +
        theme(
          plot.title = element_text(
            size = 20, face = "bold", family = "sans",
            hjust = 0, vjust = 1,
            margin = margin(b = 15, t = 10), lineheight = 0.8
          ),
          plot.title.position = "plot",
          plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
          legend.position = "right",
          text = element_text(family = "sans"),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_blank()
        )
    } else {
      sub_p <- ggplot() +
        geom_segment(
          data = non_self_edges,
          aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
              color = Interaction_Type_simplified),
          alpha = 0.5, linewidth = 0.9,
          arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
        ) +
        geom_curve(
          data = self_loops,
          aes(x = x_start_adj, y = y_start_adj, xend = x_end_adj, yend = y_end_adj,
              color = Interaction_Type_simplified),
          curvature = -1, angle = 90, alpha = 0.5, linewidth = 0.9,
          arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 30)
        ) +
        scale_color_manual(name = "Edge Type", values = edge_type_colors, na.value = "grey80") +
        ggnewscale::new_scale_color() +
        ggnewscale::new_scale("linetype") +
        scatterpie::geom_scatterpie(
          data = nodes,
          aes(x = x, y = y, group = name, colour = node_border_color, linetype = node_border_linetype),
          cols = pie_cols, pie_scale = 1.1, lwd = 0.85
        ) +
        ggrepel::geom_text_repel(
          data = nodes,
          aes(x = x, y = y, label = name),
          size = 4, color = "grey30", family = "sans", fontface = "bold",
          box.padding = 1, max.overlaps = Inf, min.segment.length = 0,
          segment.color = "grey50", segment.size = 0.5, force = 2
        ) +
        scale_fill_manual(
          name = "Protein/PhosphoProtein Expression",
          values = fill_values,
          breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
          drop = FALSE
        ) +
        scale_linetype_identity(
          name = "Node Border", guide = "legend",
          breaks = c("longdash", "solid"),
          labels = c("Disease-related", "Disease-Notrelated")
        ) +
        scale_colour_identity(
          name = "Node Border", guide = "legend",
          breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
          labels = c("Disease-related", "Disease-Notrelated")
        ) +
        coord_fixed() +
        theme_void() +
        labs(title = title_str) +
        theme(
          plot.title = element_text(
            size = 20, face = "bold", family = "sans",
            hjust = 0, vjust = 1,
            margin = margin(b = 15, t = 10), lineheight = 0.8
          ),
          plot.title.position = "plot",
          plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
          legend.position = "right",
          text = element_text(family = "sans"),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_blank()
        )
    }
    
    # Save plot
    if (layout_type == 'tree') {
      ggsave(filename = paste0("pro_phosphopro_network_", network_mode, ".png"),
             plot = sub_p, width = 16, height = 12, units = "in", dpi = 300)
    } else {
      ggsave(filename = paste0("pro_phosphopro_network_", network_mode, ".png"),
             plot = sub_p, width = 15, height = 12, units = "in", dpi = 300)
    }
  }
  
  return("✅ Success!")
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

# Validate parameters
opt <- parameter_validation(opt)

# Construct subnetwork
topnodes_subgraph <- subnetwork_construction(opt, top_nodes_visualization_num)

# Visualize if network has enough nodes
if (is.character(topnodes_subgraph)) {
  print(paste0("❌ No subnetworks can be filtered under ", network_mode_new, " mode."))
} else if (length(V(topnodes_subgraph)$name) > 2) {
  network_visualiztion(topnodes_subgraph, opt)
} else {
  print(paste0("❌ No subnetworks can be filtered under ", network_mode_new, " mode."))
}

#===============================================================================
# FUNCTIONAL ENRICHMENT FOR THE SUBNETWORK
#===============================================================================

# Load node data and extract protein/phosphosite lists
f_Node_rawdata <- read.csv(opt$nodes, sep = '\t', stringsAsFactors = FALSE)
pro_data_raw <- f_Node_rawdata[f_Node_rawdata$NodeName %in% V(topnodes_subgraph)$name, ]

pro_list <- unique(pro_data_raw[pro_data_raw$Pro_class != "None", ]$NodeName)
phos_list <- unique(pro_data_raw[pro_data_raw$Phosphopro_class != "None", ]$NodeName)

print(file.path(opt$outdir, "Functional_Enrichment"))

# Run enrichment analysis
omics_enrichment_list(pro_list, phos_list, file.path(opt$outdir, "Functional_Enrichment"),
                      omics1_name = opt$omics1_name, omics2_name = opt$omics2_name,
                      enrich_fromType = opt$enrich_fromType,
                      pvalueCutoff = 0.05, GO_showCategory = 6, KEGG_showCategory = 15,
                      color_gradient_low = opt$function_enrich_color_gradient_low,
                      color_gradient_high = opt$function_enrich_color_gradient_high)