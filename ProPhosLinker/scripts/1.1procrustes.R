#' Perform Procrustes Analysis and Visualization
#'
#' This function performs Procrustes analysis on two sets of multi-omics data
#' (e.g., proteomics and phosphoproteomics) to evaluate the consistency between
#' datasets in reduced dimensional space, and generates visualization plots.
#'
#' @param pro_path Character string, path to proteomics data file (TSV format).
#'                 Tab-delimited, first column contains sample names, first row contains protein names.
#' @param phos_path Character string, path to phosphoproteomics data file (TSV format).
#'                  Tab-delimited, first column contains sample names, first row contains phosphosite names.
#' @param group_path Character string, path to group information file (TSV format),
#'                   must contain 'sample' and 'group' columns.
#' @param out_dir Character string, path to output directory for results.
#' @param dim_rd_mtd Character string, dimensionality reduction method, either "PCA" or "PCoA".
#' @param omics1_label Character string, label for the first omics dataset (e.g., legend title).
#' @param omics2_label Character string, label for the second omics dataset (e.g., legend title).
#' @param point_size Numeric, size of points in the plot, default is 3.
#' @param PA_group1_color Character string, color for proteomics data points, default is "#a03c32".
#' @param PA_group2_color Character string, color for phosphoproteomics data points, default is "#1a5f6e".
#'
#' @return No return value, generates analysis result files and visualization plots.
#' @export
#'
#' @examples
#' \dontrun{
#' Procrustes_Analysis(
#'   pro_path = 'protein_data_normalized.tsv',
#'   phos_path = 'phosphopro_data_normalized.tsv',
#'   group_path = 'compare_groups.tsv',
#'   out_dir = './results/',
#'   dim_rd_mtd = 'PCA',
#'   omics1_label = "Pro",
#'   omics2_label = "Phos"
#' )
#' }

# Suppress package startup messages for cleaner output
options(warn = -1)
suppressWarnings(suppressPackageStartupMessages(library(readr)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(vegan)))
suppressWarnings(suppressPackageStartupMessages(library(ggrepel)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

#' Main function for Procrustes analysis between two omics datasets
#'
#' @param pro_path Path to proteomics data file
#' @param phos_path Path to phosphoproteomics data file
#' @param group_path Path to group information file
#' @param out_dir Output directory path
#' @param dim_rd_mtd Dimensionality reduction method ("PCA" or "PCoA")
#' @param omics1_label Label for first omics dataset
#' @param omics2_label Label for second omics dataset
#' @param point_size Point size in visualization
#' @param PA_group1_color Color for first group
#' @param PA_group2_color Color for second group
#'
#' @return None
Procrustes_Analysis <- function(pro_path, phos_path, group_path, out_dir, dim_rd_mtd, 
                                omics1_label, omics2_label, point_size = 3,
                                PA_group1_color = "#a03c32", PA_group2_color = "#1a5f6e") {
  
  #=============================================================================
  # 1. DATA LOADING
  #=============================================================================
  
  # Load proteomics data - samples in rows, proteins in columns
  pro_data <- read.csv(pro_path, sep = '\t', row.names = 1)
  
  # Load phosphoproteomics data - samples in rows, phosphosites in columns
  phosphopro_data <- read.csv(phos_path, sep = '\t', row.names = 1)
  
  # Load group information and select relevant columns
  group_data <- read.csv(group_path, sep = '\t')[c('sample', 'group')]
  
  #=============================================================================
  # 2. PROCRUSTES ANALYSIS FUNCTION
  #=============================================================================
  
  #' Perform Procrustes analysis on two omics datasets
  #'
  #' @param omics1 First omics dataset (samples in columns, features in rows)
  #' @param omics2 Second omics dataset (samples in columns, features in rows)
  #' @param dim_rd_mtd Dimensionality reduction method ("PCA" or "PCoA")
  #' @return Procrustes analysis result object from protest()
  procrutest <- function(omics1, omics2, dim_rd_mtd) {
    # Transpose to have samples in rows, features in columns (required for vegan functions)
    omics1 <- t(omics1)
    omics2 <- t(omics2)
    
    #' Principal Component Analysis (PCA) for continuous data
    #' @param x Data matrix with samples in rows, variables in columns
    #' @return PCA scores for first two components
    pca_metab <- function(x) {
      # Perform unconstrained PCA using vegan's rda function
      x_pca <- rda(x)
      # Extract sample coordinates for first two principal components
      x_pca <- scores(x_pca, choices = 1:2)
      return(x_pca)
    }
    
    #' Principal Coordinate Analysis (PCoA) for compositional data
    #' @param x Data matrix with samples in rows, variables in columns
    #' @return PCoA scores for first two coordinates
    PCoA_micro <- function(x) {
      # Calculate Bray-Curtis dissimilarity matrix
      micro_dist <- vegdist(x, method = "bray")
      # Perform weighted classical multidimensional scaling (PCoA)
      micro_pcoa <- wcmdscale(micro_dist)
      return(micro_pcoa)
    }
    
    # Apply selected dimensionality reduction method to both datasets
    if (dim_rd_mtd == "PCA") {
      matrix1 <- pca_metab(omics1)
      matrix2 <- pca_metab(omics2)
    } else {
      matrix1 <- PCoA_micro(omics1)
      matrix2 <- PCoA_micro(omics2)
    }
    
    # Perform Procrustes analysis with protest (significance test)
    pro_result <- protest(matrix1, matrix2)
    return(pro_result)
  }
  
  #=============================================================================
  # 3. PROCRUSTES VISUALIZATION FUNCTION
  #=============================================================================
  
  #' Create Procrustes analysis visualization plot
  #'
  #' @param pro Procrustes result object from protest()
  #' @param group Data frame with group information
  #' @param point_size Size of points in plot
  #' @param omics1 Label for first omics dataset
  #' @param omics2 Label for second omics dataset
  #' @return ggplot2 object
  plot_procrustes <- function(pro, group, point_size = 1, omics1, omics2) {
    # Extract first and last group labels for legend ordering
    group1 <- group$group[1]
    group2 <- group$group[length(group$group)]
    
    # Prepare data frame for visualization
    # Yrot: rotated coordinates of second dataset (target matrix)
    # X: original coordinates of first dataset (reference matrix)
    plot_data <- data.frame(
      rda1 = pro$Yrot[, 1],  # X-coordinate of rotated second dataset
      rda2 = pro$Yrot[, 2],  # Y-coordinate of rotated second dataset
      xrda1 = pro$X[, 1],    # X-coordinate of first dataset
      xrda2 = pro$X[, 2],    # Y-coordinate of first dataset
      group = group$group,    # Group information
      samplename = group$sample  # Sample names for labels
    )
    
    # Save Procrustes coordinates to file
    write_tsv(plot_data, file.path(out_dir, "procrustes_results.tsv"))
    
    # Create statistical annotation text
    pro_stat_text <- glue::glue(
      "Procrustes Analysis\nM\u00B2 = {round(pro$ss, digits = 3)}, r = {round(sqrt(1 - pro$ss), digits = 3)}, p = {round(pro$signif, digits = 3)}"
    )
    
    # Create ggplot visualization
    p <- ggplot(plot_data) +
      
      ## Plot first dataset points (reference: circles)
      geom_point(
        aes(x = rda1, y = rda2, colour = group, shape = omics1),
        size = point_size
      ) +
      
      ## Plot second dataset points (target: triangles)
      geom_point(
        aes(x = xrda1, y = xrda2, colour = group, shape = omics2),
        size = point_size
      ) +
      
      ## Add arrows showing the shift from first to second dataset
      geom_segment(
        aes(x = rda1, y = rda2, xend = xrda1, yend = xrda2, colour = group),
        linetype = "dashed",
        arrow = arrow(length = unit(0.2, "cm")),
        arrow.fill = NULL  # Arrow color follows group color
      ) +
      
      ## Add sample labels with repulsion to avoid overlap
      geom_text_repel(
        aes(x = rda1, y = rda2, label = samplename, colour = group),
        max.overlaps = Inf,          # Show all labels
        show.legend = FALSE,          # Don't show labels in legend
        size = 2,                     # Font size
        fontface = "bold",             # Bold for readability
        alpha = 0.4,                   # Transparency
        min.segment.length = 0,        # Always show connecting line
        segment.size = 0.3,             # Thin line
        segment.color = "grey50"        # Gray connecting line
      ) +
      
      ## Manual shape scale for legend
      scale_shape_manual(
        name = "Point Type",
        values = setNames(c(16, 17), c(omics1, omics2)),  # 16 = circle, 17 = triangle
        breaks = c(omics1, omics2)
      ) +
      
      ## Manual color scale for groups
      scale_colour_manual(
        name = "Group Type",
        values = setNames(c(PA_group1_color, PA_group2_color), c(group1, group2)),
        breaks = c(group1, group2)
      ) +
      
      ## Add Procrustes statistics annotation in top-left corner
      annotate(
        "text",
        x = -Inf,                     # Left-aligned
        y = Inf,                       # Top-aligned
        label = pro_stat_text,
        hjust = -0.05,                  # Slight inward offset
        vjust = 1.1,                    # Slight downward offset
        size = 6,                       # Font size
        colour = "black",
        fontface = "bold"
      ) +
      
      ## Labels and theme
      labs(
        colour = NULL,                   # Remove color legend title
        shape = NULL,                    # Remove shape legend title
        x = "Dimension 1",
        y = "Dimension 2"
      ) +
      
      ## Legend guides
      guides(
        colour = guide_legend(ncol = 1),  # Single column for color legend
        shape = guide_legend(ncol = 1, override.aes = list(size = point_size))
      ) +
      
      theme_classic()                     # Clean classic theme
    
    return(p)
  }
  
  #=============================================================================
  # 4. EXECUTE PROCRUSTES ANALYSIS AND VISUALIZATION
  #=============================================================================
  
  # Perform Procrustes analysis
  procrustes_result <- procrutest(pro_data, phosphopro_data, dim_rd_mtd)
  
  # Generate visualization plot
  p <- plot_procrustes(procrustes_result, group_data, 
                       point_size = point_size, 
                       omics1_label, omics2_label)
  
  # Save plot to file
  ggsave(file.path(out_dir, "procrustes.png"), p, width = 10, height = 6.5)
  
  # Extract and save statistical results
  pa_results <- c(M2 = procrustes_result$ss, p_value = procrustes_result$signif)
  write_tsv(as.data.frame(t(pa_results)), 
            file.path(out_dir, "procrustes_stat.tsv"), 
            col_names = TRUE)
}