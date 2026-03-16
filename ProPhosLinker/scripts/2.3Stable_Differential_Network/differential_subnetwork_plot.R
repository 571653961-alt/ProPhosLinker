#' Plot Differential Subnetworks for Omics Data.
#'
#' This function reads nodes and edges data from specified paths, processes the 
#' network based on correlation values between two groups (e.g., 'T' vs 'N'), 
#' and plots differential subnetworks. It generates three main plots: one for each 
#' group's specific network and a combined differential network plot. The plots 
#' are saved in PNG format at the given directory path.
#'
#' @param subnet_dir Character string specifying the directory containing nodes and edges files.
#' @param subnet_name Character string representing the name of the subnetwork.
#' @param group_name Character string indicating the comparison group names separated by "-vs-", default is "T-vs-N".
#' @param omics1_name Character string for the first type of omics data name, default is "Pro".
#' @param omics2_name Character string for the second type of omics data name, default is "Phos".
#' @param edge_color_pos Color hex code for positive correlation edges, default is "#9b6a65".
#' @param edge_color_neg Color hex code for negative correlation edges, default is "#5d8992".
#' @param Enhanced_in_N Color hex code for edges enhanced in N, default is "#5d8992".
#' @param Enhanced_in_T Color hex code for edges enhanced in T, default is "#9b6a65".
#' @param Only_in_N Color hex code for edges only present in N, default is "#0c2b32".
#' @param Only_in_T Color hex code for edges only present in T, default is "#381512".
#' @param Conflict_relation Color hex code for conflict relations, default is '#808080'.
#' @param fill_gradientn_color Vector of color hex codes for gradient fills, default is c("#175663", "#dce6e9", "#90362d").
#'
#' @importFrom igraph graph_from_data_frame delete_edges
#' @importFrom ggraph ggraph geom_edge_link geom_node_point scale_edge_color_manual scale_fill_gradientn scale_shape_manual theme_void theme ggtitle guides
#' @importFrom tidyverse "%>%"
#' @importFrom dplyr case_when
#' @importFrom grDevices colorRampPalette
#' @importFrom patchwork plot_layout
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' Differential_subnetwork_plot(subnet_dir = "/path/to/subnets",
#'                             subnet_name = "example_subnet",
#'                             group_name = "T-vs-N",
#'                             omics1_name = "Pro",
#'                             omics2_name = "Phos")
#' }
#'
#' @export


options(warn = -1)
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))

Differential_subnetwork_plot <- function(subnet_dir,subnet_name,group_name = 'T-vs-N', omics1_name = 'Pro',
                                         omics2_name = 'Phos',
                                         edge_color_pos = "#9b6a65",
                                         edge_color_neg = "#5d8992", 
                                         Enhanced_in_N = "#5d8992", 
                                         Enhanced_in_T = "#9b6a65",
                                         Only_in_N = "#0c2b32",
                                         Only_in_T = "#381512",
                                         Conflict_relation = '#808080',
                                         fill_gradientn_color = c("#175663", "#dce6e9", "#90362d")){
nodes_path <- file.path(subnet_dir,paste0("nodes_",subnet_name,"_",group_name,"_",subnet_name,".txt"))
edges_path <- file.path(subnet_dir,paste0("edges_",subnet_name,"_",group_name,"_",subnet_name,".txt"))

group_parts <- strsplit(group_name, "-vs-")[[1]]
group1 <- group_parts[1]
group2 <- group_parts[2]

nodes <- read.table(nodes_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
edges <- read.table(edges_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)


V(g)$case_norm_mean <- nodes$case_norm_mean
V(g)$omics_name <- nodes$omics_name


E(g)$cor_case <- edges$cor_case
E(g)$cor_T_type <- case_when(
  E(g)$cor_case < 0 ~ "< 0",
  E(g)$cor_case > 0 ~ "> 0",
  TRUE ~ "other"
)


g_filtered <- delete_edges(g, E(g)[cor_T_type == "other"])


layout_matrix <- as.matrix(nodes[, c("x", "y")])

current_mapping <- setNames(c(21, 24), c(omics1_name, omics2_name))

T_subnet <-ggraph(g_filtered, layout = layout_matrix) +
  geom_edge_link(
    aes(color = cor_T_type),
    width = 0.6,
    alpha = 0.7
  ) +
  geom_node_point(
    aes(fill = case_norm_mean,shape = Class),
    size = 10,
    colour = "transparent",
    stroke = 0
  ) +
  scale_edge_color_manual(
    name = "Correlation",
    values = c("< 0" = edge_color_neg , "> 0" = edge_color_pos, "other" = "transparent"),
    labels = c("< 0", "> 0", "Other")
  ) +
  scale_fill_gradientn(
    name = "Case Norm Mean",
    colours = colorRampPalette(fill_gradientn_color)(50)
  )+
  scale_shape_manual(
    name = "Omics Type",
    values =current_mapping, 
    guide = guide_legend(override.aes = list(fill = "grey50", size = 6)) 
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 60, face = "bold"),
    legend.title = element_text(size = 36, face = "bold"),     
    legend.text = element_text(size = 36),                   
  ) +
  ggtitle(paste0("Group ",group1," net")) +
  guides(
    fill = guide_colorbar(barheight = unit(4, "cm")),
   shape = guide_legend(
     override.aes = list(
       fill = "grey80",
       size = 6,
       stroke = 0.5
     ))
    )
T_subnet_lable <- T_subnet +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 14,
    fontface = "bold",
    color = "grey30",       
    family = "sans",        
    max.overlaps = 30,      
    box.padding = 2,        
    point.padding = 0.8,    
    force = 10,              
    min.segment.length = 0.1, 
    segment.color = "grey70", 
    segment.size = 0.3,     
    segment.alpha = 0.7     
  ) 
T_subnet_lable
ggsave(file.path(subnet_dir,paste0(subnet_name,"_",strsplit(group_name,"-vs-")[[1]][1],"_plot.png")),plot = T_subnet_lable, width = 10, height = 8, dpi = 300)


E(g)$cor_control <- edges$cor_control
E(g)$cor_N_type <- case_when(
  E(g)$cor_control < 0 ~ "< 0",
  E(g)$cor_control > 0 ~ "> 0",
  TRUE ~ "other"
)

g_filtered <- delete_edges(g, E(g)[cor_N_type == "other"])

N_subnet <-ggraph(g_filtered, layout = layout_matrix) +
  geom_edge_link(
    aes(color = cor_N_type),
    width = 0.6,
    alpha = 0.7
  ) +
  geom_node_point(
    aes(fill = control_norm_mean,shape = Class),
    size = 10,
    colour = "transparent", 
    stroke = 0
  ) +
  scale_edge_color_manual(
    name = "Correlation",
    # values = c("< 0" = "#487f8b", "> 0" = "#b3635b", "other" = "gray"),
    # values = c("< 0" = "#487f8b", "> 0" = "#b3635b", "other" = "transparent"),
    values = c("< 0" = edge_color_neg , "> 0" = edge_color_pos, "other" = "transparent"),
    labels = c("< 0", "> 0", "Other")
  ) +
  scale_fill_gradientn(
    name = "Case Norm Mean",
    colours = colorRampPalette(fill_gradientn_color)(50)
  )+
  scale_shape_manual(
    name = "Omics Type",
    values = current_mapping,  
    guide = guide_legend(override.aes = list(fill = "grey50", size = 4)) 
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 60, face = "bold"),
    legend.title = element_text(size = 36, face = "bold"),     
    legend.text = element_text(size = 36),             
  ) +
  ggtitle(paste0("Group ",group2," net")) +
  guides(
    fill = guide_colorbar(barheight = unit(4, "cm")),
    shape = guide_legend(
      override.aes = list(
        fill = "grey80",
        size = 6,
        stroke = 0.5
      ))
  )
N_subnet_lable <- N_subnet+
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 14,
    fontface = "bold",
    color = "grey30",       
    family = "sans",        
    max.overlaps = 30,     
    box.padding = 2,     
    point.padding = 0.8,   
    force = 10,          
    min.segment.length = 0.1, 
    segment.color = "grey70",
    segment.size = 0.3,     
    segment.alpha = 0.7    
  ) 
ggsave(file.path(subnet_dir,paste0(subnet_name,"_",strsplit(group_name,"-vs-")[[1]][2],"_plot.png")),plot = N_subnet_lable, width = 10, height = 8, dpi = 300)



g_filtered <- delete_edges(g, E(g)[cor_status == "Non-significant"])



edge_color_vector <- c(Enhanced_in_N, Enhanced_in_T, Only_in_N, Only_in_T, Conflict_relation)

names(edge_color_vector) <- c(
  paste0("Enhanced in ", group2),
  paste0("Enhanced in ", group1), 
  paste0("Only in ", group2),
  paste0("Only in ", group1),
  "Conflict relation"
)

edge_linetype_vector <- c( "dashed","dashed", "dotted", "solid","solid","dotted")
names(edge_linetype_vector) <- c(
  paste0("Enhanced in ", group2),
  paste0("Enhanced in ", group1), 
  "Non-significant",
  paste0("Only in ", group2),
  paste0("Only in ", group1),
  "Conflict relation"
)



diff_subnet <-ggraph(g_filtered, layout = layout_matrix) +
  geom_edge_link(
    aes(color = cor_status,
        linetype = cor_status),
    width = 0.6,
    alpha = 0.5
  ) +
  geom_node_point(
    aes(fill = Log2FC,shape = Class),
    size = 10,
    colour = "transparent",
    stroke = 0
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 16,
    fontface = "bold",
    color = "grey30",       
    family = "sans",       
    max.overlaps = 30,      
    box.padding = 2,      
    point.padding = 0.8,  
    force = 10,             
    min.segment.length = 0.1,
    segment.color = "grey70", 
    segment.size = 0.3,     
    segment.alpha = 0.7    
  ) +
  scale_edge_color_manual(
    name = "Correlation",
    values = edge_color_vector
  ) +
  scale_edge_linetype_manual(
    name = "Correlation",
    values = edge_linetype_vector
  ) +
  scale_shape_manual(
    name = "Omics Type",
    values = current_mapping, 
    guide = guide_legend(override.aes = list(fill = "grey50", size = 4))  
  ) +
  scale_fill_gradientn(
    name = paste0("Log2FC (",group1,":",group2,")"),
    colours = colorRampPalette(fill_gradientn_color)(50)
  )+
  theme_void() +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 60, face = "bold"),
    legend.title = element_text(size = 36, face = "bold"),    
    legend.text = element_text(size = 36),               
  ) +
  ggtitle(paste0(group_name," net")) +
  guides(
    fill = guide_colorbar(barheight = unit(4, "cm")),
    shape = guide_legend(
      override.aes = list(
        fill = "grey80",
        size = 6,
        stroke = 0.5
      ))
  )
ggsave(file.path(subnet_dir,paste0(subnet_name,"_diff_plot.png")),plot = diff_subnet, width = 10, height = 8, dpi = 300)


combined_plot <- T_subnet + N_subnet + diff_subnet+ plot_layout(ncol = 3)
ggsave(file.path(subnet_dir,paste0(subnet_name,"_",group_name,"_plot.png")), plot = combined_plot, width = 30, height = 8, dpi = 300)

}
