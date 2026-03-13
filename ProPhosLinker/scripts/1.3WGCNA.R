#' WGCNA: Weighted Gene Co-expression Network Analysis for Protein and Phosphoprotein Data
#' 
#' This function performs Weighted Gene Co-expression Network Analysis (WGCNA) for
#' proteomics and phosphoproteomics data. It includes soft-threshold selection,
#' network construction, module identification, module correlation analysis,
#' and trait association analysis.
#'
#' @param pro_path Character vector, path to proteomics data file (preprocessed).
#'                 TSV format, tab-delimited, first column contains sample names,
#'                 first row contains protein names.
#' @param phos_path Character vector, path to phosphoproteomics data file (preprocessed).
#'                  TSV format, tab-delimited, first column contains sample names,
#'                  first row contains phosphosite names.
#' @param metadata_path Character vector, path to metadata file containing sample
#'                      clinical or phenotypic information. TSV format, tab-delimited.
#'                      Set to "NULL" if no metadata available.
#' @param out_dir Character vector, output directory path. Default is current directory.
#'                Analysis results will be saved here.
#' 
#' @param omics1_name Character vector, name for proteomics dataset, default 'Pro'
#' @param omics2_name Character vector, name for phosphoproteomics dataset, default 'Phos'
#' 
#' @param protein_cor_method Character vector, correlation method for protein data.
#'                           Options: "pearson", "spearman", "kendall".
#' @param protein_corFun_tmp Character vector, WGCNA correlation function for protein data.
#'                           Usually "bicor" or "cor".
#' @param protein_cluster_method Character vector, clustering method for protein data.
#'                              e.g., "average", "complete", "ward.D", etc.
#' @param protein_corOptions_list List, options list for protein correlation calculation.
#' @param protein_corOptions_str Character vector, options string for protein correlation calculation.
#' @param protein_networkType Character vector, network type for protein data.
#'                            Options: "signed", "unsigned", "signed hybrid".
#' @param protein_RsquareCut_val Numeric, R-squared cutoff for protein soft-threshold selection.
#'                               Usually ranges 0.80-0.95.
#' @param protein_mergingThresh Numeric, module merging threshold for protein data.
#'                              Maximum dissimilarity of module eigengenes (1-correlation).
#' @param protein_minModuleSize Integer, minimum module size for protein data.
#'                              Minimum number of proteins constituting a module.
#' @param protein_SoftPower Numeric or NULL, soft-threshold power for protein data.
#'                          If NULL, automatically calculated.
#' 
#' @param phosphoprotein_cor_method Character vector, correlation method for phosphoprotein data.
#' @param phosphoprotein_corFun_tmp Character vector, WGCNA correlation function for phosphoprotein data.
#' @param phosphoprotein_cluster_method Character vector, clustering method for phosphoprotein data.
#' @param phosphoprotein_corOptions_list List, options list for phosphoprotein correlation calculation.
#' @param phosphoprotein_corOptions_str Character vector, options string for phosphoprotein correlation calculation.
#' @param phosphoprotein_networkType Character vector, network type for phosphoprotein data.
#' @param phosphoprotein_RsquareCut_val Numeric, R-squared cutoff for phosphoprotein soft-threshold selection.
#' @param phosphoprotein_mergingThresh Numeric, module merging threshold for phosphoprotein data.
#' @param phosphoprotein_minModuleSize Integer, minimum module size for phosphoprotein data.
#' @param phosphoprotein_SoftPower Numeric or NULL, soft-threshold power for phosphoprotein data.
#' 
#' @param module_cor_threshhold Numeric, threshold for module correlation significance.
#' @param module_cor_p_adj Numeric, adjusted p-value threshold for module correlation tests.
#' @param module_cor_method Character vector, correlation method for module analysis.
#' 
#' @param WGCNA_pro_filter_num Integer, number of top proteins to retain based on ANOVA F-statistic.
#' @param WGCNA_phos_filter_num Integer, number of top phosphosites to retain based on ANOVA F-statistic.
#' 
#' @param pro_ME_color Character string, color for proteomics module eigengenes in plots.
#' @param phos_ME_color Character string, color for phosphoproteomics module eigengenes in plots.
#' @param pro_color Character string, color for proteomics data in plots.
#' @param phos_color Character string, color for phosphoproteomics data in plots.
#' @param pheatmap_color Character vector, color palette for pheatmap visualizations.
#'
#' @return No direct return value. The function generates the following output files:
#'   - Soft-threshold selection plots for protein and phosphoprotein networks
#'   - Module correlation heatmap showing relationships between protein and phosphoprotein modules
#'   - Module-trait association heatmap (if metadata provided)
#'   - Sample dendrograms with phenotype heatmaps for both datasets (if metadata provided)
#'   - Module member gene lists for significantly correlated modules
#'   - Table of significantly correlated module pairs
#'   
#' @export
#' @import WGCNA
#' @import ggplot2
#' @import patchwork
#' @import pheatmap
#' @import dplyr
#' @import plyr
#' @import readr
#' @import viridis
#' @import grid
#' 
#===============================================================================
# Load Required Libraries
#===============================================================================
suppressWarnings(suppressPackageStartupMessages(library(WGCNA)))   # Clustering software for weighted gene co-expression network analysis
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))  # Grammar of graphics for data visualization
suppressWarnings(suppressPackageStartupMessages(library(patchwork))) # Combine multiple ggplot2 plots
suppressWarnings(suppressPackageStartupMessages(library(pheatmap)))  # Pretty heatmaps
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))     # Data manipulation with pipe-friendly syntax
suppressWarnings(suppressPackageStartupMessages(library(plyr)))      # Data transformations and splitting
suppressWarnings(suppressPackageStartupMessages(library(readr)))     # Fast reading of rectangular data
suppressWarnings(suppressPackageStartupMessages(library(viridis)))   # Colorblind-friendly color palettes
suppressWarnings(suppressPackageStartupMessages(library(grid)))      # Grid graphics for custom annotations

#===============================================================================
# Main WGCNA Analysis Function
#===============================================================================

#' Perform comprehensive WGCNA analysis on paired proteomics and phosphoproteomics datasets
#'
#' This function orchestrates the entire WGCNA workflow for two omics datasets:
#' 1. Data loading and preprocessing
#' 2. Feature selection based on ANOVA F-statistic
#' 3. Soft-threshold power selection for network construction
#' 4. Network construction and module detection
#' 5. Module-trait association analysis (if metadata provided)
#' 6. Cross-omics module correlation analysis
#' 7. Visualization and result export
#'
#' @param ... See parameter descriptions above
#' @return None, generates files in output directory
WGCNA_Analysis <- function(pro_path, phos_path, sample_group, metadata_path, out_dir = './',
                           omics1_name = 'Pro', omics2_name = 'Phos',
                           protein_cor_method = 'spearman', protein_corFun_tmp = 'bicor',
                           protein_cluster_method = 'average', protein_corOptions_list = list(use = 'pairwise.complete.obs'),
                           protein_corOptions_str = "use = 'pairwise.complete.obs'",
                           protein_networkType = 'signed',
                           protein_RsquareCut_val = 0.890,
                           protein_mergingThresh = 0.20, protein_minModuleSize = 60, protein_SoftPower = NULL,
                           phosphoprotein_cor_method = 'spearman', phosphoprotein_corFun_tmp = 'bicor',
                           phosphoprotein_cluster_method = 'average',
                           phosphoprotein_corOptions_list = list(use = 'pairwise.complete.obs'),
                           phosphoprotein_corOptions_str = "use = 'pairwise.complete.obs'",
                           phosphoprotein_networkType = 'signed',
                           phosphoprotein_RsquareCut_val = 0.890,
                           phosphoprotein_mergingThresh = 0.20, phosphoprotein_minModuleSize = 60,
                           phosphoprotein_SoftPower = NULL,
                           module_cor_threshhold = 0.8, module_cor_p_adj = 0.05, module_cor_method = 'spearman',
                           WGCNA_pro_filter_num = 5000,
                           WGCNA_phos_filter_num = 5000,
                           pro_ME_color = "#01344F",
                           phos_ME_color = "#D12128",
                           pro_color = "#a03c32",
                           phos_color = "#1a5f6e",
                           pheatmap_color = c("#1a5f6e", "white", "#a03c32")) {
  
  #=============================================================================
  # 1. INITIALIZATION AND SETUP
  #=============================================================================
  
  # Set global options
  options(stringsAsFactors = FALSE)
  
  # Close any existing connections and enable multi-threading for WGCNA
  closeAllConnections()  # Close all existing connections before parallel code
  enableWGCNAThreads()   # Enable multi-threading for faster computation
  
  # Initialize results vector
  wgcna_results <- c()
  
  #=============================================================================
  # 2. DATA LOADING AND PREPROCESSING
  #=============================================================================
  
  # Load sample group information
  sample_group <- read.table(sample_group, header = TRUE, row.names = 1, sep = '\t', 
                             check.names = FALSE)
  
  # Load proteomics abundance data
  input_pro_abundance <- read.table(pro_path, header = TRUE, row.names = 1, sep = '\t')
  
  #-----------------------------------------------------------------------------
  # 2.1 Feature Selection for Proteomics Data (if too many features)
  #-----------------------------------------------------------------------------
  if (nrow(input_pro_abundance) > WGCNA_pro_filter_num) {
    # Log2 transform the data (adding pseudocount to avoid log(0))
    log_input <- log2(input_pro_abundance + 1)
    
    # Find common samples between expression data and group information
    common_samples <- intersect(colnames(log_input), rownames(sample_group))
    if (length(common_samples) == 0) {
      stop("No common samples between expression matrix and sample_group")
    }
    
    # Reorder both expression matrix columns and sample_group rows to match
    common_samples <- common_samples[order(match(common_samples, colnames(log_input)))]
    log_input <- log_input[, common_samples, drop = FALSE]
    sample_group_sub <- sample_group[common_samples, , drop = FALSE]
    
    # Verify sample alignment
    if (!all(colnames(log_input) == rownames(sample_group_sub))) {
      stop("Sample names are not aligned after subsetting: check column/rownames of expression and sample_group")
    }
    
    # Extract group factor for ANOVA
    group <- factor(sample_group_sub[, "group"])
    
    # Calculate one-way ANOVA F-statistic for each protein
    # This identifies proteins that vary significantly between groups
    f_values <- apply(log_input, 1, function(x) {
      # Check sample length consistency
      if (length(x) != length(group)) {
        stop("Sample length mismatch between expression row and group vector")
      }
      anova_model <- aov(x ~ group)
      anova_summary <- summary(anova_model)
      f_stat <- anova_summary[[1]]$"F value"[1]
      return(f_stat)
    })
    
    # Sort by F-value (descending) and select top proteins
    top_proteins <- names(sort(f_values, decreasing = TRUE)[1:WGCNA_pro_filter_num])
    
    # Filter original abundance matrix to keep only top proteins
    input_pro_abundance <- input_pro_abundance[top_proteins, ]
  }
  
  # Transpose data: convert to format with samples in rows, proteins in columns (required for WGCNA)
  proteinData <- t(input_pro_abundance)
  
  #-----------------------------------------------------------------------------
  # 2.2 Load and Process Phosphoproteomics Data (similar to proteomics)
  #-----------------------------------------------------------------------------
  input_phosphopro_abundance <- read.table(phos_path, header = TRUE, row.names = 1, sep = '\t')
  
  if (nrow(input_phosphopro_abundance) > WGCNA_phos_filter_num) {
    log_input <- log2(input_phosphopro_abundance + 1)
    common_samples <- intersect(colnames(log_input), rownames(sample_group))
    if (length(common_samples) == 0) {
      stop("No common samples between expression matrix and sample_group")
    }
    
    common_samples <- common_samples[order(match(common_samples, colnames(log_input)))]
    log_input <- log_input[, common_samples, drop = FALSE]
    sample_group_sub <- sample_group[common_samples, , drop = FALSE]
    
    if (!all(colnames(log_input) == rownames(sample_group_sub))) {
      stop("Sample names are not aligned after subsetting: check column/rownames of expression and sample_group")
    }
    
    group <- factor(sample_group_sub[, "group"])
    
    # Calculate ANOVA F-statistic for phosphosites
    f_values <- apply(log_input, 1, function(x) {
      if (length(x) != length(group)) {
        stop("Sample length mismatch between expression row and group vector")
      }
      anova_model <- aov(x ~ group)
      anova_summary <- summary(anova_model)
      f_stat <- anova_summary[[1]]$"F value"[1]
      return(f_stat)
    })
    
    # Select top phosphosites
    top_phosproteins <- names(sort(f_values, decreasing = TRUE)[1:WGCNA_phos_filter_num])
    
    input_phosphopro_abundance <- input_phosphopro_abundance[top_phosproteins, ]
  }
  
  # Transpose phosphoprotein data
  phosphoproteinData <- t(input_phosphopro_abundance)
  
  #=============================================================================
  # 3. WGCNA ANALYSIS
  #=============================================================================
  
  #-----------------------------------------------------------------------------
  # 3.1 Adjust correlation method for small sample sizes
  #-----------------------------------------------------------------------------
  if (length(rownames(proteinData)) < 30) {
    if (module_cor_method != 'kendall') {
      module_cor_method = 'kendall'  # Kendall's tau is more robust for small samples
    }
  }
  
  #-----------------------------------------------------------------------------
  # 3.2 Soft-Threshold Power Selection Function
  #-----------------------------------------------------------------------------
  
  #' Determine optimal soft-threshold power for network construction
  #'
  #' @param omicsData Expression data matrix (samples × features)
  #' @param RsquareCut_val R-squared cutoff for scale-free topology fit
  #' @param networkType Type of network ("signed", "unsigned", etc.)
  #' @param datatype Name of data type for plot labeling
  #' @return Optimal soft-threshold power value
  get_SoftPower_val <- function(omicsData, RsquareCut_val, networkType, datatype) {
    # Set color based on data type
    if (datatype == omics1_name) {
      datatype_color = pro_color
    } else {
      datatype_color = phos_color
    }
    
    # Define powers to test
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    nSamples = length(rownames(omicsData))
    
    # Calculate soft-thresholding criteria
    sft <- pickSoftThreshold(omicsData, powerVector = powers, verbose = 1)
    
    # Prepare data for plotting
    df <- data.frame(
      Power = sft$fitIndices[, 1],
      R2 = sft$fitIndices[, 2],
      Connectivity = sft$fitIndices[, 5]
    )
    
    # Create scale independence plot
    p1 <- ggplot(df, aes(Power, R2)) +
      geom_text(aes(label = Power), color = datatype_color, size = 4, fontface = "bold") +
      geom_hline(yintercept = RsquareCut_val, linetype = "dashed", color = datatype_color) +
      annotate(
        "text",
        x = max(df$Power) * 0.95,
        y = RsquareCut_val * 0.9,
        label = paste("R² =", RsquareCut_val),
        color = datatype_color,
        size = 5,
        hjust = 1,
        fontface = "bold"
      ) +
      labs(title = paste0("Scale independence (", datatype, ")"),
           x = "Soft Threshold (power)",
           y = expression(paste("Scale Free Topology Model Fit, ", R^2))) +
      theme_minimal(base_size = 12)
    
    # Create mean connectivity plot
    p2 <- ggplot(df, aes(Power, Connectivity)) +
      geom_text(aes(label = Power), color = datatype_color, size = 4, fontface = "bold") +
      labs(title = paste0("Mean Connectivity (", datatype, ")"),
           x = "Soft Threshold (power)",
           y = "Mean Connectivity") +
      theme_minimal(base_size = 12)
    
    # Combine plots
    combined <- p1 + p2 + plot_layout(ncol = 2)
    
    # Save plot
    ggsave(file.path(out_dir, paste0(datatype, "_soft_threshold.png")), 
           combined, width = 9, height = 5, dpi = 300)
    
    # Extract estimated power
    softPower = sft$powerEstimate
    
    # If no power meets criteria, use heuristic based on sample size
    if (is.na(softPower)) {
      power = ifelse(nSamples < 20, ifelse(networkType == "unsigned", 9, 18),
                     ifelse(nSamples < 30, ifelse(networkType == "unsigned", 8, 16),
                            ifelse(nSamples < 40, ifelse(networkType == "unsigned", 7, 14),
                                   ifelse(networkType == "unsigned", 20, 20))))
    }
    
    return(softPower)
  }
  
  #-----------------------------------------------------------------------------
  # 3.3 Determine Soft Thresholds for Both Datasets
  #-----------------------------------------------------------------------------
  if (is.null(protein_SoftPower)) {
    protein_SoftPower <- get_SoftPower_val(proteinData, protein_RsquareCut_val, 
                                           protein_networkType, omics1_name)
    wgcna_results['protein_SoftPower'] <- protein_SoftPower
  }
  if (is.null(phosphoprotein_SoftPower)) {
    phosphoprotein_SoftPower <- get_SoftPower_val(phosphoproteinData, 
                                                  phosphoprotein_RsquareCut_val, 
                                                  phosphoprotein_networkType, omics2_name)
    wgcna_results['phosphoprotein_SoftPower'] <- phosphoprotein_SoftPower
  }
  
  # Save soft-power values to file
  write_tsv(as.data.frame(t(wgcna_results)), 
            file.path(out_dir, "WGCNA_softpower.tsv"), col_names = TRUE)
  
  #-----------------------------------------------------------------------------
  # 3.4 Network Construction Function
  #-----------------------------------------------------------------------------
  
  #' Construct co-expression network and identify modules
  #'
  #' @param omicsData Expression data matrix (samples × features)
  #' @param SoftPower Soft-threshold power
  #' @param networkType Network type ("signed", "unsigned", etc.)
  #' @param corFun_tmp Correlation function to use
  #' @param corOptions_str Correlation options as string
  #' @param corOptions_list Correlation options as list
  #' @param cluster_method Hierarchical clustering method
  #' @param minClusterSize Minimum module size
  #' @param mergingThresh Module merging threshold
  #' @param datatype Data type for labeling
  #' @return List containing module labels, colors, and eigengenes
  network_construction <- function(omicsData, SoftPower, networkType, corFun_tmp, 
                                   corOptions_str, corOptions_list, cluster_method, 
                                   minClusterSize, mergingThresh, datatype) {
    # Calculate adjacency matrix (co-expression similarity)
    A <- adjacency(omicsData, power = SoftPower, type = networkType, 
                   corFnc = corFun_tmp, corOptions = corOptions_str)
    colnames(A) = rownames(A) = colnames(omicsData)
    
    # Calculate topological overlap matrix (TOM) based dissimilarity
    disTOM <- TOMdist(A, TOMType = networkType)
    colnames(disTOM) = rownames(disTOM) = colnames(omicsData)
    
    # Hierarchical clustering on TOM dissimilarity
    metaTree <- flashClust::flashClust(as.dist(disTOM), method = cluster_method)
    
    # Dynamic tree cutting to identify modules
    moduleLabels1 <- cutreeDynamic(dendro = metaTree,
                                   distM = disTOM,
                                   method = "hybrid",
                                   deepSplit = 4,      # Maximum splitting depth (0-4, 4 most sensitive)
                                   pamRespectsDendro = T,  # Preserve dendrogram structure
                                   minClusterSize = minClusterSize  # Minimum module size
    )
    
    # Merge similar modules based on eigengene correlation
    merge <- mergeCloseModules(omicsData,
                               moduleLabels1,
                               corFnc = corFun_tmp,
                               corOptions = corOptions_list,
                               cutHeight = mergingThresh)
    
    moduleColor <- merge$colors
    MEs <- merge$newMEs
    names(moduleColor) <- colnames(omicsData)
    rownames(MEs) = rownames(omicsData)
    
    return(list(moduleLabels1 = moduleLabels1, moduleColor = moduleColor, MEs = MEs))
  }
  
  #-----------------------------------------------------------------------------
  # 3.5 Construct Networks for Both Datasets
  #-----------------------------------------------------------------------------
  protein_network <- network_construction(proteinData, protein_SoftPower, 
                                          protein_networkType, protein_corFun_tmp, 
                                          protein_corOptions_str, protein_corOptions_list, 
                                          protein_cluster_method, protein_minModuleSize, 
                                          protein_mergingThresh, "protein")
  
  phosphoprotein_network <- network_construction(phosphoproteinData, phosphoprotein_SoftPower,
                                                 phosphoprotein_networkType, phosphoprotein_corFun_tmp,
                                                 phosphoprotein_corOptions_str, phosphoprotein_corOptions_list,
                                                 phosphoprotein_cluster_method, phosphoprotein_minModuleSize,
                                                 phosphoprotein_mergingThresh, "phosphoprotein")
  
  # Extract results from protein network
  protein_moduleLabels1 <- protein_network$moduleLabels1
  protein_moduleColor <- protein_network$moduleColor
  protein_MEs <- protein_network$MEs
  
  # Extract results from phosphoprotein network
  phosphoprotein_moduleLabels1 <- phosphoprotein_network$moduleLabels1
  phosphoprotein_moduleColor <- phosphoprotein_network$moduleColor
  phosphoprotein_MEs <- phosphoprotein_network$MEs
  
  #=============================================================================
  # 4. CROSS-OMICS MODULE CORRELATION ANALYSIS
  #=============================================================================
  
  # Rename module eigengenes to include omics prefix
  colnames(protein_MEs) <- paste0(omics1_name, "_", colnames(protein_MEs))
  colnames(phosphoprotein_MEs) <- paste0(omics2_name, "_", colnames(phosphoprotein_MEs))
  
  # Combine module eigengenes from both datasets
  all_MEs <- cbind(protein_MEs, phosphoprotein_MEs)
  
  # Calculate correlations between all modules
  module_cor <- cor(all_MEs, use = "pairwise.complete.obs")
  
  # Calculate p-values and adjust for multiple testing
  module_p <- corPvalueStudent(module_cor, nrow(all_MEs))
  module_p_adj <- p.adjust(module_p, method = "BH")
  
  # Reshape adjusted p-values to matrix format
  module_p_adj_matrix <- matrix(module_p_adj,
                                nrow = nrow(module_p),
                                ncol = ncol(module_p),
                                dimnames = dimnames(module_p))
  
  #-----------------------------------------------------------------------------
  # 4.1 Create Module Correlation Heatmap
  #-----------------------------------------------------------------------------
  plot_matrix <- module_cor
  
  # Calculate module sizes for annotation
  pro_counts <- table(protein_moduleColor)
  pro_module_names <- paste0(omics1_name, "_ME", names(pro_counts))
  pro_counts <- as.numeric(pro_counts)
  names(pro_counts) <- pro_module_names
  
  # Create color gradient for proteomics modules (based on module size)
  pro_color_fun <- colorRampPalette(c(pro_ME_color, "white"))
  pro_colors <- pro_color_fun(length(pro_counts))
  pro_colors <- pro_colors[rank(-pro_counts)]  # Darker colors for larger modules
  names(pro_colors) <- names(pro_counts)
  
  # Calculate phospho module sizes
  phospho_counts <- table(phosphoprotein_moduleColor)
  phospho_module_names <- paste0(omics2_name, "_ME", names(phospho_counts))
  phospho_counts <- as.numeric(phospho_counts)
  names(phospho_counts) <- phospho_module_names
  
  # Create color gradient for phosphoproteomics modules
  phospho_color_fun <- colorRampPalette(c(phos_ME_color, "white"))
  phospho_colors <- phospho_color_fun(length(phospho_counts))
  phospho_colors <- phospho_colors[rank(-phospho_counts)]
  names(phospho_colors) <- names(phospho_counts)
  
  # Create annotation data frame for heatmap
  annotation_col <- data.frame(
    ModuleNumber = c(names(pro_counts), names(phospho_counts)),
    ModuleType = c(rep("Proteome", length(pro_counts)), rep("Phosphoproteome", length(phospho_counts))),
    Counts = c(pro_counts, phospho_counts),
    Color = c(pro_colors, phospho_colors),
    row.names = c(names(pro_counts), names(phospho_counts))
  )
  
  # Ensure annotation matches matrix column order
  annotation_col <- annotation_col[colnames(plot_matrix), ]
  
  # Define annotation colors for heatmap
  annotation_colors <- list(
    ModuleType = c(Proteome = pro_ME_color, Phosphoproteome = phos_ME_color),
    ModuleNumber = c(pro_colors, phospho_colors)
  )
  
  annotation_col <- annotation_col[, c('ModuleType', 'ModuleNumber')]
  
  # Calculate dynamic plot height based on number of modules
  cell_height <- 12  # mm per cell
  total_height_mm <- nrow(plot_matrix) * cell_height + 2 + 1 + 0
  height_inch <- total_height_mm / (25.4 * 2) + 2
  
  # Save heatmap to PNG
  png(file.path(out_dir, "module_correlation_pheatmap.png"),
      width = height_inch + 3.5,
      height = height_inch,
      units = "in",
      res = 300)
  
  pheatmap(
    plot_matrix,
    color = colorRampPalette(pheatmap_color)(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_row = annotation_col,
    show_colnames = TRUE,
    angle_col = 45,
    angle_row = 45,
    main = paste0(omics1_name, " vs. ", omics2_name, " Modules Relationships\n"),
    fontsize = 9,
    fontsize_row = 8,
    fontsize_col = 8,
    cellwidth = 12,
    cellheight = 12,
    zlim = c(-1, 1),
    border_color = NA,
    annotation_legend = FALSE
  )
  
  # Add custom legend annotations using grid graphics
  grid.text("Module Type", x = 0.95, y = 0.95, just = "right", gp = gpar(fontface = "bold"))
  grid.rect(x = 0.95, y = 0.90, width = 0.04, height = 0.04, just = "right", 
            gp = gpar(fill = pro_ME_color, col = NA))
  grid.text(omics1_name, x = 0.90, y = 0.90, just = "right", gp = gpar(fontsize = 8))
  grid.rect(x = 0.95, y = 0.85, width = 0.04, height = 0.04, just = "right", 
            gp = gpar(fill = phos_ME_color, col = NA))
  grid.text(omics2_name, x = 0.90, y = 0.85, just = "right", gp = gpar(fontsize = 8))
  
  # Add gradient legend for module sizes
  grid.text("Module Number", x = 0.95, y = 0.79, just = "right", gp = gpar(fontface = "bold"))
  grid.text(paste0(omics1_name, "_ME"), x = 0.86, y = 0.75, just = "right", gp = gpar(fontsize = 8))
  grid.text(paste0(omics2_name, "_ME"), x = 0.86, y = 0.69, just = "right", gp = gpar(fontsize = 8))
  
  # Helper function to draw gradient legend
  gradient_legend <- function(x, y, width, height, colors, num_element) {
    n <- 100
    colors <- colorRampPalette(colors)(n)
    
    for(i in 1:n) {
      grid.rect(
        x = x + (i-1)/n * width,
        y = y,
        width = width/n,
        height = height,
        just = c("left", "bottom"),
        gp = gpar(fill = colors[i], col = NA))
    }
    
    # Add labels at 0%, 50%, and 100%
    label_positions <- c(0, 0.5, 1)
    for(pos in label_positions) {
      grid.text(
        as.character(pos * num_element),
        x = x + pos * width,
        y = y - unit(0.5, "lines"),
        just = c("center", "top"),
        gp = gpar(fontsize = 7)
      )
    }
  }
  
  # Add proteomics module size gradient
  gradient_legend(
    x = unit(0.87, "npc"),
    y = unit(0.74, "npc"),
    width = unit(0.08, "npc"),
    height = unit(0.02, "npc"),
    colors = c("white", pro_ME_color),
    num_element = ncol(proteinData)
  )
  
  # Add phosphoproteomics module size gradient
  gradient_legend(
    x = unit(0.87, "npc"),
    y = unit(0.68, "npc"),
    width = unit(0.08, "npc"),
    height = unit(0.02, "npc"),
    colors = c("white", phos_ME_color),
    num_element = ncol(phosphoproteinData)
  )
  
  dev.off()
  
  #-----------------------------------------------------------------------------
  # 4.2 Filter Significantly Correlated Module Pairs
  #-----------------------------------------------------------------------------
  
  #' Filter module pairs based on correlation threshold and adjusted p-value
  #'
  #' @param threshold Correlation threshold
  #' @param p_adj Adjusted p-value threshold
  #' @return Data frame of significant module pairs
  module_cor_filtering <- function(threshold, p_adj) {
    # Find indices of significant correlations
    sig_pairs <- which(abs(module_cor) > threshold & module_p_adj < p_adj, arr.ind = TRUE)
    
    # Create data frame with module pair information
    core_modules <- data.frame(
      Module1 = colnames(all_MEs)[sig_pairs[, 1]],
      Module2 = colnames(all_MEs)[sig_pairs[, 2]],
      Correlation = module_cor[sig_pairs],
      p_adj = module_p_adj[sig_pairs]
    ) %>%
      filter(Module1 != Module2) %>%  # Remove self-correlations
      arrange(desc(abs(Correlation)))
    
    # Filter to keep only cross-omics pairs (different prefixes)
    core_modules_filtered <- core_modules %>%
      mutate(
        Prefix1 = sapply(strsplit(Module1, "_"), `[`, 1),
        Prefix2 = sapply(strsplit(Module2, "_"), `[`, 1)
      ) %>%
      filter(Prefix1 != Prefix2) %>%
      select(-Prefix1, -Prefix2)
    
    # Exclude modules labeled as "ME0" (grey/unassigned)
    core_modules_filtered_NOTEM0 <- core_modules_filtered %>%
      filter(!grepl("_ME0", .[[1]]),
             !grepl("_ME0", .[[2]]))
    
    # Remove duplicate module pairs (keep the one with smallest p-value)
    core_modules_filtered_NOTEM0_noduplicate <- core_modules_filtered_NOTEM0 %>%
      mutate(
        pair_key = purrr::map2_chr(Module1, Module2, ~ paste(sort(c(.x, .y)), collapse = "|"))
      ) %>%
      group_by(pair_key) %>%
      arrange(p_adj) %>%
      slice(1) %>%
      ungroup() %>%
      select(-pair_key)
    
    return(core_modules_filtered_NOTEM0_noduplicate)
  }
  
  # Apply filtering to get significant module pairs
  module_cor_filtering_info <- module_cor_filtering(module_cor_threshhold, module_cor_p_adj)
  
  # Save filtered results
  write_tsv(
    module_cor_filtering_info,
    file.path(out_dir, "module_correlation_filtering.tsv")
  )
  
  # Extract unique module IDs from significant pairs
  module_cor_filtering_modulids <- sort(unique(c(
    as.vector(module_cor_filtering_info[, "Module1"])$Module1,
    as.vector(module_cor_filtering_info[, "Module2"])$Module2
  )), decreasing = TRUE)
  
  # Save feature lists for each significant module
  for (module_cor_filtering_modulid in module_cor_filtering_modulids) {
    con <- file(file.path(out_dir, paste0(module_cor_filtering_modulid, "_list.tsv")), "w")
    module_type <- strsplit(module_cor_filtering_modulid, '_ME')[[1]][1]
    module_num <- strsplit(module_cor_filtering_modulid, '_ME')[[1]][2]
    
    if (module_type == omics1_name) {
      module_cor_filtering_modulid_list <- unique(names(protein_moduleColor)[protein_moduleLabels1 == module_num])
    } else {
      module_cor_filtering_modulid_list <- unique(names(phosphoprotein_moduleColor)[phosphoprotein_moduleLabels1 == module_num])
    }
    
    writeLines(module_cor_filtering_modulid_list, con)
    close(con)
  }
  
  #=============================================================================
  # 5. METADATA ANALYSIS (if provided)
  #=============================================================================
  
  if (length(metadata_path) != 0) {
    
    #-----------------------------------------------------------------------------
    # 5.1 Load and Prepare Metadata
    #-----------------------------------------------------------------------------
    metadata <- read.csv(metadata_path, sep = '\t')
    common_id <- intersect(metadata$sample, rownames(proteinData))
    
    # Subset all datasets to common samples
    metadata <- metadata[metadata$sample %in% common_id, ]
    proteinData <- proteinData[rownames(proteinData) %in% common_id, ]
    phosphoproteinData <- phosphoproteinData[rownames(phosphoproteinData) %in% common_id, ]
    
    # Set row names from sample column
    metadata_rownames <- metadata$sample
    metadata <- metadata[, -1]
    rownames(metadata) <- metadata_rownames
    
    # Keep only numeric columns for correlation analysis
    metadata <- metadata[, sapply(metadata, is.numeric)]
    
    #-----------------------------------------------------------------------------
    # 5.2 Sample Dendrogram with Metadata Heatmap - Proteomics
    #-----------------------------------------------------------------------------
    sampleTree2 = hclust(dist(proteinData), method = "average")
    
    # Create color gradient for metadata visualization
    customColorRamp <- colorRampPalette(c(pro_color, "white"))
    
    # Generate color matrix for metadata values
    traitColorMatrix <- apply(metadata, 2, function(col) {
      col_numeric <- as.numeric(col)
      normalized <- (col_numeric - min(col_numeric, na.rm = TRUE)) /
        (max(col_numeric, na.rm = TRUE) - min(col_numeric, na.rm = TRUE))
      customColorRamp(100)[as.numeric(cut(normalized, breaks = 100, include.lowest = TRUE))]
    })
    
    traitColors <- as.data.frame(traitColorMatrix)
    rownames(traitColors) <- rownames(metadata)
    
    # Truncate long column names for better visualization
    old_rownames <- names(metadata)
    grouplabels <- ifelse(
      nchar(old_rownames) > 15,
      paste0(substr(old_rownames, 1, 15), "..."),
      old_rownames
    )
    
    # Clear any open graphics devices
    while (!is.null(dev.list())) dev.off()
    
    # Save dendrogram with trait heatmap
    png(file.path(out_dir, paste0(omics1_name, "_sample_dendrogram_and_metadata_heatmap.png")),
        width = 12, height = 6, units = "in", res = 300)
    plotDendroAndColors(sampleTree2, traitColors,
                        groupLabels = grouplabels,
                        main = paste0(omics1_name, " Sample Dendrogram and Phenotype Heatmap"),
                        cex.dendroLabels = 0.4,
                        cex.rowText = 0.6)
    dev.off()
    
    #-----------------------------------------------------------------------------
    # 5.3 Sample Dendrogram with Metadata Heatmap - Phosphoproteomics
    #-----------------------------------------------------------------------------
    sampleTree2 = hclust(dist(phosphoproteinData), method = "average")
    
    customColorRamp <- colorRampPalette(c(phos_color, "white"))
    
    traitColorMatrix <- apply(metadata, 2, function(col) {
      col_numeric <- as.numeric(col)
      normalized <- (col_numeric - min(col_numeric, na.rm = TRUE)) /
        (max(col_numeric, na.rm = TRUE) - min(col_numeric, na.rm = TRUE))
      customColorRamp(100)[as.numeric(cut(normalized, breaks = 100, include.lowest = TRUE))]
    })
    
    traitColors <- as.data.frame(traitColorMatrix)
    rownames(traitColors) <- rownames(metadata)
    
    old_rownames <- names(metadata)
    grouplabels <- ifelse(
      nchar(old_rownames) > 15,
      paste0(substr(old_rownames, 1, 15), "..."),
      old_rownames
    )
    
    while (!is.null(dev.list())) dev.off()
    png(file.path(out_dir, paste0(omics2_name, "_Sample_Dendrogram_and_Metadata_Heatmap.png")),
        width = 12, height = 6, units = "in", res = 300)
    plotDendroAndColors(sampleTree2, traitColors,
                        groupLabels = grouplabels,
                        main = paste0(omics2_name, " Sample Dendrogram and Phenotype Heatmap"),
                        cex.dendroLabels = 0.4,
                        cex.rowText = 0.6)
    dev.off()
    
    #-----------------------------------------------------------------------------
    # 5.4 Module-Trait Correlation Analysis
    #-----------------------------------------------------------------------------
    
    # Check if metadata rownames match module eigengene rownames
    if (setequal(rownames(metadata), rownames(all_MEs))) {
      # Calculate correlations between modules and traits
      moduleTraitCor = cor(all_MEs, metadata, use = "p")
    } else {
      # If samples don't match, recalculate networks with matching samples
      if (length(rownames(proteinData)) < 30) {
        if (module_cor_method != 'kendall') {
          module_cor_method = 'kendall'
        }
      }
      
      # Recalculate soft thresholds and networks with subsetted data
      protein_SoftPower <- get_SoftPower_val(proteinData, protein_RsquareCut_val, 
                                             protein_networkType, paste0(omics1_name, "_metadata"))
      phosphoprotein_SoftPower <- get_SoftPower_val(phosphoproteinData, 
                                                    phosphoprotein_RsquareCut_val, 
                                                    phosphoprotein_networkType, 
                                                    paste0(omics2_name, "_metadata"))
      
      protein_network <- network_construction(proteinData, protein_SoftPower, 
                                              protein_networkType, protein_corFun_tmp,
                                              protein_corOptions_str, protein_corOptions_list,
                                              protein_cluster_method, protein_minModuleSize,
                                              protein_mergingThresh, "protein")
      
      phosphoprotein_network <- network_construction(phosphoproteinData, phosphoprotein_SoftPower,
                                                     phosphoprotein_networkType, phosphoprotein_corFun_tmp,
                                                     phosphoprotein_corOptions_str, phosphoprotein_corOptions_list,
                                                     phosphoprotein_cluster_method, phosphoprotein_minModuleSize,
                                                     phosphoprotein_mergingThresh, "phosphoprotein")
      
      # Extract results
      protein_moduleLabels1 <- protein_network$moduleLabels1
      protein_moduleColor <- protein_network$moduleColor
      protein_MEs <- protein_network$MEs
      
      phosphoprotein_moduleLabels1 <- phosphoprotein_network$moduleLabels1
      phosphoprotein_moduleColor <- phosphoprotein_network$moduleColor
      phosphoprotein_MEs <- phosphoprotein_network$MEs
      
      # Rename and combine module eigengenes
      colnames(protein_MEs) <- paste0(omics1_name, "_", colnames(protein_MEs))
      colnames(phosphoprotein_MEs) <- paste0(omics2_name, "_", colnames(phosphoprotein_MEs))
      all_MEs <- cbind(protein_MEs, phosphoprotein_MEs)
      
      # Calculate module-trait correlations
      moduleTraitCor = cor(all_MEs, metadata, use = "p")
    }
    
    # Calculate p-values for module-trait correlations
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(proteinData))
    
    # Truncate trait names for visualization
    truncated_yLabels <- ifelse(nchar(names(metadata)) > 18,
                                paste0(substr(names(metadata), 1, 15), "..."),
                                names(metadata))
    
    # Prepare matrix for heatmap
    plot_matrix <- t(moduleTraitCor)
    rownames(plot_matrix) <- truncated_yLabels
    colnames(plot_matrix) <- names(all_MEs)
    
    # Create text matrix with correlation and p-value
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    
    text_matrix <- matrix(sprintf("%.2f", plot_matrix),
                          nrow = nrow(plot_matrix))
    
    # Prepare annotation for module-trait heatmap
    pro_counts <- table(protein_moduleColor) / sum(table(protein_moduleColor))
    pro_module_names <- paste0(omics1_name, "_ME", names(pro_counts))
    pro_counts <- as.numeric(pro_counts)
    names(pro_counts) <- pro_module_names
    
    pro_color_fun <- colorRampPalette(c(pro_ME_color, "white"))
    pro_colors <- pro_color_fun(length(pro_counts))
    pro_colors <- pro_colors[rank(-pro_counts)]
    names(pro_colors) <- names(pro_counts)
    
    phospho_counts <- table(phosphoprotein_moduleColor) / sum(table(phosphoprotein_moduleColor))
    phospho_module_names <- paste0(omics2_name, "_ME", names(phospho_counts))
    phospho_counts <- as.numeric(phospho_counts)
    names(phospho_counts) <- phospho_module_names
    
    phospho_color_fun <- colorRampPalette(c(phos_ME_color, "white"))
    phospho_colors <- phospho_color_fun(length(phospho_counts))
    phospho_colors <- phospho_colors[rank(-phospho_counts)]
    names(phospho_colors) <- names(phospho_counts)
    
    # Create annotation data frame
    annotation_col <- data.frame(
      ModuleNumber = c(names(pro_counts), names(phospho_counts)),
      ModuleType = c(rep("Proteome", length(pro_counts)), rep("Phosphoproteome", length(phospho_counts))),
      Counts = c(pro_counts, phospho_counts),
      Color = c(pro_colors, phospho_colors),
      row.names = c(names(pro_counts), names(phospho_counts))
    )
    
    annotation_col <- annotation_col[colnames(plot_matrix), ]
    
    annotation_colors <- list(
      ModuleType = c(Proteome = pro_ME_color, Phosphoproteome = phos_ME_color),
      ModuleNumber = c(pro_colors, phospho_colors)
    )
    
    annotation_col <- annotation_col[, c('ModuleType', 'ModuleNumber')]
    
    # Find split point between proteomics and phosphoproteomics modules
    first_phos_index <- which(annotation_col$ModuleType == "Phosphoproteome")[1]
    if (is.na(first_phos_index)) {
      gaps_col <- 0  # All proteome
    } else if (first_phos_index == 1) {
      gaps_col <- 0  # No proteome
    } else {
      gaps_col <- first_phos_index - 1
    }
    
    # Calculate dynamic plot dimensions
    cell_height <- 12
    total_height_mm <- nrow(plot_matrix) * cell_height * 2 + 1
    height_inch <- total_height_mm / (25.4) + 0.5
    total_width_mm <- ncol(plot_matrix) * cell_height + 2 + 2
    widthinch <- total_width_mm / (25.4 * 2)
    
    # Save module-trait correlation heatmap
    png(file.path(out_dir, "module_metadata_correlation_pheatmap.png"),
        width = widthinch, height = height_inch, units = "in", res = 300)
    
    pheatmap(
      plot_matrix,
      color = colorRampPalette(pheatmap_color)(50),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      gaps_col = gaps_col,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      show_colnames = TRUE,
      angle_col = 45,
      angle_row = 45,
      main = paste0("Module vs. Metadata Relationships \n"),
      fontsize = 8,
      fontsize_row = 8,
      fontsize_col = 8,
      cellwidth = 12,
      cellheight = 12,
      zlim = c(-1, 1),
      border_color = NA,
      annotation_legend = FALSE
    )
    
    # Add custom legends
    grid.text("Module Type", x = 0.95, y = 0.95, just = "right", gp = gpar(fontface = "bold"))
    grid.rect(x = 0.95, y = 0.90, width = 0.04, height = 0.04, just = "right", 
              gp = gpar(fill = pro_ME_color, col = NA))
    grid.text(omics1_name, x = 0.90, y = 0.90, just = "right", gp = gpar(fontsize = 10))
    grid.rect(x = 0.95, y = 0.85, width = 0.04, height = 0.04, just = "right", 
              gp = gpar(fill = phos_ME_color, col = NA))
    grid.text(omics2_name, x = 0.90, y = 0.85, just = "right", gp = gpar(fontsize = 10))
    
    grid.text("Module Number", x = 0.95, y = 0.79, just = "right", gp = gpar(fontface = "bold"))
    grid.text(paste0(omics1_name, "_ME"), x = 0.86, y = 0.75, just = "right", gp = gpar(fontsize = 10))
    grid.text(paste0(omics2_name, "_ME"), x = 0.86, y = 0.69, just = "right", gp = gpar(fontsize = 10))
    
    # Add gradient legends
    gradient_legend(
      x = unit(0.87, "npc"),
      y = unit(0.74, "npc"),
      width = unit(0.08, "npc"),
      height = unit(0.02, "npc"),
      colors = c("white", pro_ME_color),
      num_element = ncol(proteinData)
    )
    
    gradient_legend(
      x = unit(0.87, "npc"),
      y = unit(0.68, "npc"),
      width = unit(0.08, "npc"),
      height = unit(0.02, "npc"),
      colors = c("white", phos_ME_color),
      num_element = ncol(phosphoproteinData)
    )
    
    dev.off()
  }
}

