# community_detection_plot
# Description: Detects community-related activities in network traffic.

suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(tidygraph)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggforce)))
suppressWarnings(suppressPackageStartupMessages(library(ggpubr)))
suppressWarnings(suppressPackageStartupMessages(library(tibble)))



#' Detect and Plot Community Structures in a Differential Network
#'
#' This function loads node and edge data from specified file paths, filters nodes based on membership criteria,
#' detects communities within the network, and generates visualizations of the top modules. The resulting plots
#' are saved to the specified output directory.
#'
#' @param node_file_path Character string specifying the path to the node data file (tab-separated).
#' @param edge_file_path Character string specifying the path to the edge data file (tab-separated).
#' @param outdir Character string specifying the directory where the output plots will be saved.
#' @param omics1_name Optional character string naming the first type of omics data for legend labels (default: 'Pro').
#' @param omics2_name Optional character string naming the second type of omics data for legend labels (default: 'Phos').
#' @param ModuleSize_show Numeric threshold for filtering out small modules by size (default: 0, no filtering).
#' @param top_module_num Integer specifying the number of top modules to display in the plot (default: 20).
#'
#' @return Invisibly returns NULL. Outputs are saved as PNG files in the specified output directory:
#'   - \code{network_by_top_modules.png}: A plot showing the entire network with highlighted top modules.
#'   - \code{network_by_separate_top_modules.png}: Separate plots for each of the top modules.
#'
#' @importFrom igraph graph_from_data_frame induced_subgraph centrality_degree
#' @importFrom tidygraph as_tbl_graph activate mutate
#' @importFrom ggraph ggraph geom_edge_link geom_mark_hull geom_node_point scale_color_manual scale_fill_manual scale_size_continuous theme_void labs
#' @importFrom ggplot2 aes element_text element_rect theme margin
#' @importFrom ggforce geom_mark_hull
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis
#' @export

diff_net_community_detection_plot <- function(node_file_path,edge_file_path,outdir, omics1_name = 'Pro',omics2_name = 'Phos',ModuleSize_show = 0,top_module_num = 20){
  
  # 0.Load node and edge data
  node_df <- read.csv(node_file_path, sep = "\t", header = TRUE)
  edge_df <- read.csv(edge_file_path, sep = "\t", header = TRUE)
 
  # 0.filter node and edge data
  node_df <- node_df[node_df$membership != "subnet_other",]
  filtered_nodes <- node_df$node
  edge_df <- edge_df[!(edge_df$cor_status %in% c("Non-significant")) & ((edge_df$from %in% filtered_nodes) & (edge_df$to %in% filtered_nodes)),]
  
  # 1. Create igraph object, and assign module IDs and Numbers to nodes [A whole graph]
  # Create igraph object from data frames
  g <- graph_from_data_frame(d = edge_df, vertices = node_df, directed = FALSE)
  g_tbl <- as_tbl_graph(g)
  
  # Visualize the entire graph with modules
  ### visualizing filter:: 1.ModuleSize_show, 2.top_module_num
  # ModuleSize_show <- 0 ##mask function
  
  modules_list <- table(vertex_attr(g_tbl)$membership)
  mask_modules <- names(modules_list[modules_list < ModuleSize_show])
  
  if(ModuleSize_show > 0){
    g_tbl <- g_tbl %>%
      activate(nodes) %>%
      mutate(
        membership = ifelse(
          membership %in% mask_modules, 
          paste0("other (num <", ModuleSize_show, ")"), 
          membership
        )
      )
  }
  
  # top_module_num <- 6
  modules_list <- table(vertex_attr(g_tbl)$membership)
  top_modules <- sort(modules_list,decreasing = TRUE)
  if(length(top_modules) < top_module_num){
    top_modules <- names(top_modules)
    top_module_num <- length(top_modules)
  }else{
    top_modules <- names(top_modules)[1:top_module_num]
  }

  if (length(top_modules) > 0) {
    colors_from_palettes <- c(
      RColorBrewer::brewer.pal(9, "Set1"),
      RColorBrewer::brewer.pal(8, "Dark2"),
      RColorBrewer::brewer.pal(8, "Set2"),
      RColorBrewer::brewer.pal(12, "Set3"),
      viridis::viridis(20)  
    )

    module_colors <- colors_from_palettes[1:length(top_modules)]
    color_mapping <- setNames(module_colors, top_modules)
  } else {
    color_mapping <- c()
  }
  
  top_module_node_df <- as.data.frame(vertex_attr(g_tbl))[vertex_attr(g_tbl)$membership %in% top_modules, ]
  top_module_nodes <- unique(top_module_node_df$name)

  top_module_g <- induced_subgraph(g_tbl, vids = V(g)$name %in% top_module_nodes)
  top_module_g_tbl <- as_tbl_graph(top_module_g)
  
  module_label_mapping <- table(top_module_node_df$membership)
  module_label_mapping <- module_label_mapping[names(module_label_mapping) %in% top_modules]
  
  module_label_mapping <- setNames(
    paste0(names(module_label_mapping), " (", module_label_mapping, ")"),
    names(module_label_mapping)
  )
  

  top_module_g_tbl <- top_module_g_tbl %>%
    activate(nodes) %>%
    mutate(node_color = membership,
           node_color = factor(node_color,
                               levels = c(top_modules[top_modules != paste0("other (num <", ModuleSize_show, ")")],paste0("other (num <", ModuleSize_show, ")"))),
           membership = as.factor(membership)
    )
  
  if(ModuleSize_show > 0 & (paste0("other (num <", ModuleSize_show, ")") %in% names(module_label_mapping))){
    color_mapping[paste0("other (num <", ModuleSize_show, ")")] <- "gray70"
  }
  full_color_mapping <- color_mapping
  

  layout_matrix <- as.matrix(node_df[node_df$node %in% top_module_nodes, c("x", "y")])

  top_module_p <- ggraph(top_module_g_tbl, layout = layout_matrix) +
    geom_edge_link(
      color = "grey50",
      alpha = 0.3,
      width = 0.6
    ) +
    ggforce::geom_mark_hull(
      aes(x, y,
          group = membership,
          fill = membership
      ),
      alpha = 0.1,
      concavity = 6,
      expand = unit(2, "mm"),
      radius = unit(3, "mm"),
      lwd = 0.05,
      color = 'grey90'
    ) +
    geom_node_point(
      aes(color = node_color,
          size = centrality_degree(),
          shape = Class),
      alpha = 0.9
    ) +
    scale_color_manual(
      name = "Module ID (Item Num.)",
      values = full_color_mapping,  
      drop = FALSE,
      labels = function(x) {
        sapply(x, function(id) {
          if (id %in% names(module_label_mapping)) {
            module_label_mapping[as.character(id)]
          } else {
            id
          }
        })
      }  
    ) +
    scale_fill_manual(values = full_color_mapping, na.value = "gray70",
                      guide = "none") +
    scale_size_continuous(name = "Nodes Centrality",
                          range = c(4, 8)) +
    theme_void() +
    theme(strip.text = element_text(face = "bold"),
          plot.background = element_rect(fill = "white", color = NA),
          legend.title = element_text(size = 30, face = "bold"),      
          legend.text = element_text(size = 25),                    
          plot.title = element_text(
            hjust = 0.5,           
            face = "bold",        
            size = 50,             
            margin = margin(t = 15, b = 15)  
          )) +
    labs(title = paste0("Differential Protein and Phosphorylation Sites Network with top ",top_module_num," Modules"))
  # top_module_p
  ggsave(paste0(outdir,"/network_by_top_modules.png"), plot = top_module_p, width = 10, height = 8, dpi = 300)
  
  top_module_p3 <- ggraph(top_module_g_tbl, layout = layout_matrix) +
    geom_edge_link(
      color = "grey30",
      alpha = 0.3) +
    geom_node_point(aes(color = membership,
                        size = centrality_degree(),
                        shape = Class)) +
    scale_shape_manual(
      name = "Node Type",
      values = setNames(c(16, 17), c(omics1_name, omics2_name))
    ) +
    facet_nodes(~ membership, scales = "free") +
    scale_color_manual(name = "Module ID (Item Num.)",
                       values = full_color_mapping,
                       labels = function(x) {
                         sapply(x, function(id) {
                           if (id %in% names(module_label_mapping)) {
                             module_label_mapping[as.character(id)]
                           } else {
                             id
                           }
                         })
                       },
                       na.value = "gray80") +
    scale_size_continuous(name = "Nodes Centrality",  
                          range = c(4, 8)) +
    theme_void() +
    theme( strip.text = element_text(
      face = "bold",
      size = 40,  
      margin = margin(b = 5)  
    ),
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 30, face = "bold"),      
    legend.text = element_text(size = 25),                     
    
    plot.title = element_text(
      hjust = 0.5,           
      face = "bold",         
      size = 60,             
      margin = margin(t = 15, b = 15)  
    )) +
    labs(title = paste0("Differential Protein and Phosphorylation Sites Network by top ",top_module_num," Modules"))
  # top_module_p3

  num_modules <- length(unique(vertex_attr(top_module_g_tbl, "membership")))
  expected_cols <- ceiling(sqrt(num_modules))

  ggsave(paste0(outdir,"/network_by_separate_top_modules.png"), plot = top_module_p3, width = 24, height = 6 *expected_cols, dpi = 300)
}
