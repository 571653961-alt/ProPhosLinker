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
suppressWarnings(suppressPackageStartupMessages(library(clusterProfiler)))  # Core enrichment analysis
suppressWarnings(suppressPackageStartupMessages(library(org.Hs.eg.db)))     # Human genome annotation
suppressWarnings(suppressPackageStartupMessages(library(enrichplot)))       # Visualization of enrichment results
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))          # Grammar of graphics
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))            # Data manipulation

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

#' #' Perform and visualize comparative GO enrichment analysis
#' #'
#' #' This function compares GO enrichment between two omics datasets (e.g.,
#' #' proteomics and phosphoproteomics) and generates a dot plot visualization
#' #' with faceting by GO ontology (BP, MF, CC).
#' #'
#' #' @param omics1_list Vector of feature IDs for first omics dataset
#' #' @param omics2_list Vector of feature IDs for second omics dataset
#' #' @param omics1_name Name for first omics dataset (for labeling)
#' #' @param omics2_name Name for second omics dataset (for labeling)
#' #' @param outdir Output directory for results
#' #' @param pvalueCutoff P-value cutoff for enrichment (default: 0.05)
#' #' @param showCategory Number of top categories to show per ontology (default: 15)
#' #' @param enrich_fromType ID type for enrichment ('UNIPROT' or 'SYMBOL', default: 'UNIPROT')
#' #' @param color_gradient_low Color for low end of p-value gradient (default: "#175663")
#' #' @param color_gradient_high Color for high end of p-value gradient (default: "#90362d")
#' two_omicses_GO_enrichment <- function(omics1_list = omics1_list, omics2_list = omics2_list,
#'                                       omics1_name = omics1_name, omics2_name = omics2_name,
#'                                       outdir = "./", pvalueCutoff = 0.05, showCategory = 15,
#'                                       enrich_fromType = 'UNIPROT',
#'                                       color_gradient_low = "#175663",
#'                                       color_gradient_high = "#90362d") {
#'   
#'   # Create a named list of gene clusters for comparison
#'   compare_lists <- list(omics1_list, omics2_list)
#'   names(compare_lists) <- c(omics1_name, omics2_name)
#'   
#'   # Perform comparative GO enrichment analysis
#'   go_compare <- compareCluster(
#'     geneClusters = compare_lists,
#'     fun = "enrichGO",
#'     OrgDb = org.Hs.eg.db,
#'     keyType = enrich_fromType,
#'     ont = "ALL",  # Include all three GO ontologies: BP, MF, CC
#'     pAdjustMethod = "BH",
#'     pvalueCutoff = pvalueCutoff,
#'     qvalueCutoff = 0.05
#'   )
#'   
#'   # Proceed only if enrichment results are not NULL
#'   if (!is.null(go_compare)) {
#'     
#'     # Extract and save enrichment results
#'     go_compare_data <- as.data.frame(go_compare)
#'     
#'     write.table(
#'       go_compare_data,
#'       file = paste0(outdir, "/", omics1_name, "_", omics2_name, "_GO_enrichment.tsv"),
#'       sep = "\t",
#'       row.names = FALSE,
#'       col.names = TRUE,
#'       quote = FALSE
#'     )
#'     
#'     #---------------------------------------------------------------------------
#'     # Data preprocessing for visualization
#'     #---------------------------------------------------------------------------
#'     go_compare_clean <- go_compare_data %>%
#'       # Calculate numeric GeneRatio and BgRatio for sorting
#'       mutate(
#'         GeneRatio_num = sapply(strsplit(GeneRatio, "/"), 
#'                                function(x) as.numeric(x[1]) / as.numeric(x[2])),
#'         BgRatio_num = sapply(strsplit(BgRatio, "/"), 
#'                              function(x) as.numeric(x[1]) / as.numeric(x[2])),
#'         FoldEnrichment = GeneRatio_num / BgRatio_num
#'       ) %>%
#'       # Group by ontology and sort by adjusted p-value
#'       group_by(ONTOLOGY) %>%
#'       arrange(p.adjust, .by_group = TRUE) %>%
#'       # Keep top categories per ontology
#'       slice_head(n = showCategory) %>%
#'       ungroup() %>%
#'       # Set factor levels for proper plotting order
#'       mutate(
#'         Description = factor(Description, levels = rev(unique(Description))),
#'         ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC"))
#'       ) %>%
#'       arrange(desc(GeneRatio_num), .by_group = TRUE) %>%
#'       mutate(
#'         Description = factor(Description, levels = rev(unique(Description))),
#'         ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC")),
#'         GeneRatio_factor = as.factor(GeneRatio_num)  # Optional factor version
#'       )
#'     
#'     #---------------------------------------------------------------------------
#'     # Create GO enrichment dot plot
#'     #---------------------------------------------------------------------------
#'     go_compare_GO <- ggplot(go_compare_clean, aes(x = Cluster, y = Description)) +
#'       
#'       # Add points with size = gene count, color = adjusted p-value
#'       geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
#'       
#'       # Facet by GO ontology with free y-axis scales
#'       facet_grid(ONTOLOGY ~ .,
#'                  scales = "free_y",
#'                  space = "free_y") +
#'       
#'       # Color gradient for p-values (log10 scale for better visualization)
#'       scale_color_gradient(low = color_gradient_low, high = color_gradient_high,
#'                            name = "Adjusted P-value",
#'                            trans = "log10") +
#'       
#'       # Size scale for point sizes
#'       scale_size_continuous(range = c(2, 8),
#'                             name = paste0(omics1_name, " / ", omics2_name, " Count"),
#'                             guide = guide_legend(
#'                               override.aes = list(color = "grey60")
#'                             )) +
#'       
#'       # Axis labels and title
#'       labs(x = "Protein Ratio",
#'            y = "",
#'            title = paste0("GO Enrichment Comparison: ", omics1_name, " vs ", omics2_name)) +
#'       
#'       # Theme customization
#'       theme_bw() +
#'       theme(
#'         # Title settings
#'         plot.title = element_text(hjust = 0.5, face = "bold", size = 20,
#'                                   margin = margin(b = 20)),
#'         
#'         # Facet strip settings
#'         strip.text = element_text(size = 10, face = "bold",
#'                                   margin = margin(t = 5, r = 5, b = 5, l = 5)),
#'         strip.background = element_rect(fill = "lightgray"),
#'         
#'         # Axis settings
#'         axis.text.x = element_text(size = 10, color = "grey20", face = "bold"),
#'         axis.text.y = element_text(size = 10, color = "grey20", face = "bold"),
#'         axis.title.x = element_text(size = 10, face = "bold",
#'                                     margin = margin(t = 15)),
#'         
#'         # Legend settings
#'         legend.title = element_text(size = 15, face = "bold"),
#'         legend.text = element_text(size = 10),
#'         legend.key.size = unit(0.8, "cm"),
#'         
#'         # Grid and margins
#'         panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
#'         panel.grid.minor = element_blank(),
#'         plot.margin = margin(1, 1, 1, 1, "cm")
#'       ) +
#'       
#'       # Legend customization
#'       guides(
#'         fill = guide_colorbar(barheight = unit(4, "cm")),
#'         shape = guide_legend(
#'           override.aes = list(
#'             fill = "grey80",
#'             size = 6,
#'             stroke = 0.5
#'           )
#'         )
#'       )
#'     
#'     # Save plot to file
#'     ggsave(
#'       paste0(outdir, "/", omics1_name, "_", omics2_name, "_GO_enrichment.png"),
#'       go_compare_GO,
#'       width = 12,
#'       height = 8,
#'       dpi = 300
#'     )
#'   }
#' }


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
# KEGG ENRICHMENT COMPARISON FUNCTION (Online version)
#===============================================================================

#' #' Perform and visualize comparative KEGG enrichment analysis (online)
#' #'
#' #' This function performs KEGG pathway enrichment comparison between two omics
#' #' datasets using the online KEGG database via clusterProfiler.
#' #'
#' #' @param omics1_list Vector of feature IDs for first omics dataset
#' #' @param omics2_list Vector of feature IDs for second omics dataset
#' #' @param omics1_name Name for first omics dataset
#' #' @param omics2_name Name for second omics dataset
#' #' @param outdir Output directory
#' #' @param pvalueCutoff P-value cutoff (default: 0.05)
#' #' @param showCategory Number of top categories to show (default: 15)
#' #' @param enrich_fromType ID type ('UNIPROT' or 'SYMBOL', default: 'UNIPROT')
#' #' @param color_gradient_low Low end color (default: "#175663")
#' #' @param color_gradient_high High end color (default: "#90362d")
#' two_omicses_KEGG_enrichment <- function(omics1_list = omics1_list, omics2_list = omics2_list,
#'                                         omics1_name = omics1_name, omics2_name = omics2_name,
#'                                         outdir = "./", pvalueCutoff = 0.05, showCategory = 15,
#'                                         enrich_fromType = 'UNIPROT',
#'                                         color_gradient_low = "#175663",
#'                                         color_gradient_high = "#90362d") {
#'   
#'   # Create named list of gene clusters
#'   compare_lists <- list(omics1_list, omics2_list)
#'   names(compare_lists) <- c(omics1_name, omics2_name)
#'   
#'   # Handle different ID types for KEGG enrichment
#'   if (enrich_fromType == 'UNIPROT') {
#'     # Direct KEGG enrichment with UniProt IDs
#'     kk_compare <- compareCluster(
#'       geneClusters = compare_lists,
#'       fun = "enrichKEGG",
#'       organism = "hsa",
#'       keyType = "uniprot",
#'       pvalueCutoff = pvalueCutoff
#'     )
#'   } else {
#'     # Convert SYMBOL to ENTREZID for KEGG enrichment
#'     omics1_list_entrez <- bitr(
#'       omics1_list,
#'       fromType = "SYMBOL",
#'       toType = "ENTREZID",
#'       OrgDb = org.Hs.eg.db
#'     )$ENTREZID
#'     
#'     omics2_list_entrez <- bitr(
#'       omics2_list,
#'       fromType = "SYMBOL",
#'       toType = "ENTREZID",
#'       OrgDb = org.Hs.eg.db
#'     )$ENTREZID
#'     
#'     # Update comparison lists with ENTREZ IDs
#'     compare_lists <- list(omics1_list_entrez, omics2_list_entrez)
#'     names(compare_lists) <- c(omics1_name, omics2_name)
#'     
#'     kk_compare <- compareCluster(
#'       geneClusters = compare_lists,
#'       fun = "enrichKEGG",
#'       organism = "hsa",
#'       pvalueCutoff = pvalueCutoff
#'     )
#'   }
#'   
#'   # Proceed only if enrichment results are not NULL
#'   if (!is.null(kk_compare)) {
#'     
#'     # Extract and save enrichment results
#'     kk_compare_data <- as.data.frame(kk_compare)
#'     
#'     if (nrow(kk_compare_data) > 0) {
#'       write.table(
#'         kk_compare_data,
#'         file = paste0(outdir, "/", omics1_name, "_", omics2_name, "_KEGG_enrichment.tsv"),
#'         sep = "\t",
#'         row.names = FALSE,
#'         col.names = TRUE,
#'         quote = FALSE
#'       )
#'       
#'       #---------------------------------------------------------------------------
#'       # Data preprocessing for visualization
#'       #---------------------------------------------------------------------------
#'       kk_compare_data_clean <- kk_compare_data %>%
#'         # Calculate numeric GeneRatio
#'         mutate(
#'           GeneRatio_num = sapply(strsplit(GeneRatio, "/"), 
#'                                  function(x) as.numeric(x[1]) / as.numeric(x[2])),
#'           BgRatio_num = sapply(strsplit(BgRatio, "/"), 
#'                                function(x) as.numeric(x[1]) / as.numeric(x[2])),
#'           FoldEnrichment = GeneRatio_num / BgRatio_num
#'         ) %>%
#'         # Sort and select top categories
#'         arrange(p.adjust, .by_group = TRUE) %>%
#'         slice_head(n = showCategory) %>%
#'         arrange(desc(Count), .by_group = TRUE) %>%
#'         # Set factor levels for proper ordering
#'         mutate(
#'           Description = factor(Description, levels = rev(unique(Description))),
#'           GeneRatio_factor = as.factor(GeneRatio_num)
#'         )
#'       
#'       #---------------------------------------------------------------------------
#'       # Create KEGG enrichment dot plot
#'       #---------------------------------------------------------------------------
#'       p_kk_compare <- ggplot(kk_compare_data_clean, aes(x = Cluster, y = Description)) +
#'         
#'         # Add points with size = gene count, color = adjusted p-value
#'         geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
#'         
#'         # Color gradient for p-values
#'         scale_color_gradient(low = color_gradient_low, high = color_gradient_high,
#'                              name = "Adjusted P-value",
#'                              trans = "log10") +
#'         
#'         # Size scale for points
#'         scale_size_continuous(range = c(2, 8),
#'                               name = paste0(omics1_name, " / ", omics2_name, " Count"),
#'                               guide = guide_legend(
#'                                 override.aes = list(color = "grey60")
#'                               )) +
#'         
#'         # Labels and title
#'         labs(x = "Protein Ratio",
#'              y = "",
#'              title = paste0("KEGG Enrichment Comparison: ", omics1_name, " vs ", omics2_name)) +
#'         
#'         # Theme customization
#'         theme_bw() +
#'         theme(
#'           plot.title = element_text(hjust = 0.5, face = "bold", size = 16,
#'                                     margin = margin(b = 20)),
#'           axis.text.x = element_text(size = 10, color = "grey20", face = "bold"),
#'           axis.text.y = element_text(size = 10, color = "grey20", face = "bold"),
#'           axis.title.x = element_text(size = 10, face = "bold",
#'                                       margin = margin(t = 15)),
#'           legend.title = element_text(size = 10, face = "bold"),
#'           legend.text = element_text(size = 8),
#'           legend.key.size = unit(0.8, "cm"),
#'           panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
#'           panel.grid.minor = element_blank(),
#'           plot.margin = margin(1, 1, 1, 1, "cm")
#'         ) +
#'         
#'         # Legend customization
#'         guides(
#'           fill = guide_colorbar(barheight = unit(4, "cm")),
#'           shape = guide_legend(
#'             override.aes = list(
#'               fill = "grey80",
#'               size = 6,
#'               stroke = 0.5
#'             )
#'           )
#'         )
#'       
#'       # Save plot to file
#'       ggsave(
#'         paste0(outdir, "/", omics1_name, "_", omics2_name, "_KEGG_enrichment.png"),
#'         p_kk_compare,
#'         width = 12,
#'         height = 8,
#'         dpi = 300
#'       )
#'     }
#'   }
#' }


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
# KEGG ENRICHMENT COMPARISON FUNCTION (Local database version)
#===============================================================================

#' Perform and visualize comparative KEGG enrichment analysis using local database
#'
#' This function performs KEGG pathway enrichment using a locally stored KEGG
#' database, which is useful when internet access is limited or for reproducibility.
#'
#' @param omics1_list Vector of feature IDs for first omics dataset
#' @param omics2_list Vector of feature IDs for second omics dataset
#' @param omics1_name Name for first omics dataset
#' @param omics2_name Name for second omics dataset
#' @param outdir Output directory
#' @param pvalueCutoff P-value cutoff (default: 0.05)
#' @param showCategory Number of top categories to show (default: 15)
#' @param enrich_fromType ID type ('UNIPROT' or 'SYMBOL', default: 'UNIPROT')
#' @param color_gradient_low Low end color (default: "#175663")
#' @param color_gradient_high High end color (default: "#90362d")
#' @param kegg_db_path Path to local KEGG database files
two_omicses_KEGG_enrichment_local <- function(omics1_list, omics2_list, 
                                              omics1_name, omics2_name,
                                              outdir = "./",
                                              pvalueCutoff = 0.05, 
                                              showCategory = 15,
                                              enrich_fromType = 'UNIPROT',
                                              color_gradient_low = "#175663",
                                              color_gradient_high = "#90362d",
                                              kegg_db_path = "/jdfsyt2/BC_PT/Old_Pipeline/FuctionalAnalysis/database/kegg/kegg/") {
  
  #---------------------------------------------------------------------------
  # Create comparison lists
  #---------------------------------------------------------------------------
  compare_lists <- list(omics1_list, omics2_list)
  names(compare_lists) <- c(omics1_name, omics2_name)
  
  #---------------------------------------------------------------------------
  # Convert IDs based on type
  #---------------------------------------------------------------------------
  if (enrich_fromType == 'UNIPROT') {
    # Convert UniProt to ENTREZID
    uniprot2entrez <- bitr(
      unique(unlist(compare_lists)),
      fromType = "UNIPROT",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    # Helper function to convert gene lists using mapping
    convert_cluster <- function(gene_list, map_df) {
      res <- map_df$ENTREZID[match(gene_list, map_df$UNIPROT)]
      res <- res[!is.na(res)]
      return(unique(res))
    }
    
    compare_lists_entrez <- lapply(compare_lists, convert_cluster, uniprot2entrez)
    
  } else {
    # Convert SYMBOL to ENTREZID
    symbol2entrez <- bitr(
      unique(unlist(compare_lists)),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    convert_cluster <- function(gene_list, map_df) {
      res <- map_df$ENTREZID[match(gene_list, map_df$SYMBOL)]
      res <- res[!is.na(res)]
      return(unique(res))
    }
    
    compare_lists_entrez <- lapply(compare_lists, convert_cluster, symbol2entrez)
  }
  
  #---------------------------------------------------------------------------
  # Load and process local KEGG database files
  #---------------------------------------------------------------------------
  
  # Load gene-KO (KEGG Orthology) mapping
  gene_ko_map <- read.table(
    file.path(kegg_db_path, "genes_ko.list"),
    sep = "\t",
    stringsAsFactors = FALSE
  )
  colnames(gene_ko_map) <- c("gene", "ko")
  
  # Keep only human genes (hsa prefix)
  gene_ko_map <- gene_ko_map[grepl("^hsa:", gene_ko_map$gene), ]
  gene_ko_map$gene <- sub("^hsa:", "", gene_ko_map$gene)
  gene_ko_map$ko <- sub("^ko:", "", gene_ko_map$ko)
  
  # Load KO to pathway mapping
  ko_raw <- read.table(
    file.path(kegg_db_path, "ko_map.tab"),
    sep = "\t",
    stringsAsFactors = FALSE,
    fill = TRUE,
    quote = ""
  )
  colnames(ko_raw) <- c("ko", "pathways")
  
  # Expand KO-pathway relationships (one row per KO-pathway pair)
  ko_pathway_map <- do.call(
    rbind,
    lapply(1:nrow(ko_raw), function(i) {
      data.frame(
        ko = ko_raw$ko[i],
        pathway = unlist(strsplit(ko_raw$pathways[i], " ")),
        stringsAsFactors = FALSE
      )
    })
  )
  
  # Load pathway names and descriptions
  pathway_names <- read.table(
    file.path(kegg_db_path, "map_title.tab"),
    sep = "\t",
    stringsAsFactors = FALSE,
    colClasses = c("character", "character", "character", "character")
  )
  colnames(pathway_names) <- c("pathway", "Description", "Category", "SubCategory")
  
  #---------------------------------------------------------------------------
  # Create term2gene and term2name for enrichment
  #---------------------------------------------------------------------------
  
  # Merge gene-KO and KO-pathway to get gene-pathway relationships
  term2gene <- merge(gene_ko_map, ko_pathway_map, by = "ko")
  term2gene <- term2gene[, c("pathway", "gene")]
  colnames(term2gene) <- c("term", "gene")
  
  # Create term to name mapping
  term2name <- pathway_names[, c("pathway", "Description")]
  colnames(term2name) <- c("term", "name")
  
  # Define universe as all genes in KEGG database
  kegg_universe <- unique(term2gene$gene)
  
  #---------------------------------------------------------------------------
  # Perform enrichment analysis using local data
  #---------------------------------------------------------------------------
  kk_compare <- compareCluster(
    geneClusters = compare_lists_entrez,
    fun = "enricher",
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    universe = kegg_universe,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  #---------------------------------------------------------------------------
  # Extract and visualize results
  #---------------------------------------------------------------------------
  kk_compare_data <- as.data.frame(kk_compare)
  
  # Save results
  write.table(
    kk_compare_data,
    file = paste0(outdir, "/", omics1_name, "_", omics2_name, "_KEGG_enrichment.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  # Create plot only if there are results
  if (nrow(kk_compare_data) != 0) {
    
    # Data preprocessing
    kk_compare_data_clean <- kk_compare_data %>%
      mutate(
        GeneRatio_num = sapply(strsplit(GeneRatio, "/"), 
                               function(x) as.numeric(x[1]) / as.numeric(x[2])),
        BgRatio_num = sapply(strsplit(BgRatio, "/"), 
                             function(x) as.numeric(x[1]) / as.numeric(x[2])),
        FoldEnrichment = GeneRatio_num / BgRatio_num
      ) %>%
      arrange(p.adjust, .by_group = TRUE) %>%
      slice_head(n = showCategory) %>%
      arrange(desc(Count), .by_group = TRUE) %>%
      mutate(
        Description = factor(Description, levels = rev(unique(Description))),
        GeneRatio_factor = as.factor(GeneRatio_num)
      )
    
    # Create enhanced dot plot with larger elements for better visibility
    p_kk_compare <- ggplot(kk_compare_data_clean, aes(x = Cluster, y = Description)) +
      geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
      scale_color_gradient(low = color_gradient_low, high = color_gradient_high,
                           name = "Adjusted P-value", trans = "log10") +
      scale_size_continuous(range = c(4, 12),
                            name = paste0(omics1_name, " / ", omics2_name, " Count"),
                            guide = guide_legend(override.aes = list(color = "grey60"))) +
      labs(x = "Protein Ratio", y = "",
           title = paste0("KEGG Enrichment Comparison: ", omics1_name, " vs ", omics2_name)) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24, 
                                  margin = margin(b = 20)),
        axis.text.x = element_text(size = 16, color = "grey20", face = "bold"),
        axis.text.y = element_text(size = 16, color = "grey20", face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "cm"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm")
      ) +
      guides(
        color = guide_colorbar(barheight = unit(6, "cm"),
                               title.theme = element_text(size = 16, face = "bold"),
                               label.theme = element_text(size = 14)),
        size = guide_legend(title.theme = element_text(size = 16, face = "bold"),
                            label.theme = element_text(size = 14),
                            override.aes = list(color = "grey60"))
      )
    
    # Save plot with larger dimensions
    ggsave(
      paste0(outdir, "/", omics1_name, "_", omics2_name, "_KEGG_enrichment.png"),
      p_kk_compare,
      width = 14,
      height = 10,
      dpi = 300
    )
    
  } else {
    message("❌ No KEGG enrichment found. Check gene overlap with KEGG database.")
  }
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
                                  color_gradient_high = "#90362d") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  #---------------------------------------------------------------------------
  # GO Enrichment Comparison
  #---------------------------------------------------------------------------
  # Only run if both lists have at least one gene
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
      color_gradient_high = color_gradient_high
    )
  }
  
  #---------------------------------------------------------------------------
  # KEGG Enrichment Comparison
  #---------------------------------------------------------------------------
  # Only run if both lists have at least two genes (for statistical power)
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
      color_gradient_high = color_gradient_high
    )
  }
}