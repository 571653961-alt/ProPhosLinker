#' SubNetwork Class
#'
#' An S4 class to encapsulate the results of subnetwork clustering analysis.
#' It stores metadata about the comparison group, modularity of the clustered network,
#' color mappings for visualization, the full clustered network structure, and a list
#' of extracted subnetworks (each containing nodes, edges, and optional layout information).
#'
#' @slot group_name Character string indicating the group comparison label (e.g., "T-vs-N").
#' @slot modularity Numeric value representing the modularity score of the clustering result,
#'        which quantifies the strength of community structure in the network.
#' @slot colormapping A named vector or list mapping subnetwork identifiers to colors,
#'        used for consistent coloring in downstream visualizations.
#' @slot overall_cluster_network A list containing two data frames: \code{nodes} and \code{edges},
#'        representing the full network after cluster assignment.
#' @slot subnetworks A list where each element corresponds to a significant subnetwork,
#'        typically including \code{nodes}, \code{edges}, and optionally \code{plot_layout}.
#'
#' @name SubNetwork-class
#' @rdname SubNetwork-class
#' @keywords internal

setClass("SubNetwork",
         slots = c(
           group_name= "character",
           modularity="numeric",
           colormapping = "ANY",
           overall_cluster_network = "ANY",
           subnetworks = "ANY"
         )
)

#' Cluster and Analyze Subnetworks from Nodes and Edges Data.
#'
#' This function clusters the input network data using either fast greedy or Louvain method,
#' based on the cluster size threshold. It further refines subnetworks, ensuring each contains
#' at least 5 nodes, and returns a detailed SubNetwork object encapsulating clustering results.
#'
#' @param nodes A data frame representing nodes in the network.
#' @param edges A data frame representing edges in the network.
#' @param cfg_t1 The clustering result of the network.
#' @param group_name Character string indicating the name of the group for which the analysis is performed.
#' @param diffmessage Optional parameter indicating differential message type ("NULL", "diff", or "multi").
#'
#' @return An S4 object of class "SubNetwork" containing group name, modularity score, color mapping,
#' overall cluster network details, and subnetworks information.
#'
#' @import igraph dplyr RColorBrewer bootnet qgraph
#' @export

run_subnet_cluster<-function(nodes=NULL,edges=NULL,clustersize=30){
  cor_igraph<-run_igraph(nodes,edges)
  cfg_t1 <- igraph::cluster_fast_greedy(igraph::as.undirected(cor_igraph))
  if(max(sort(table(igraph::membership(cfg_t1)), decreasing = TRUE)) > clustersize) {
    resolution <- 0.1
    repeat {
      cfg_t1 <- igraph::cluster_louvain(cor_igraph, resolution = resolution)
      max_cluster_size <- max(sort(table(igraph::membership(cfg_t1)), decreasing = TRUE))
      if (max_cluster_size <= clustersize) {
        break
      }
      resolution <- resolution + 0.1
    }
  }
  return(cfg_t1)
}


run_add_cluster<-function(cfg_t1=NULL,nodes=NULL,edges=NULL,group_name=NULL,
                          diffmessage="NULL"){
  if(length(group_name)>1){
    precor_group_name<-paste(group_name, collapse = "-vs-")
  }else{
    precor_group_name<-gsub(":","-vs-", group_name)
  }
  modularity <- igraph::modularity(cfg_t1)
  cor_igraph<-run_igraph(nodes,edges)
  t_modules_1 <- sort(table(igraph::membership(cfg_t1)),decr=T)
  t_modules<-table(igraph::membership(cfg_t1))
  membership_df <- data.frame(
    old_membership =as.character(igraph::membership(cfg_t1)),
    node = names(igraph::membership(cfg_t1))
  )
  node_names <- names(igraph::V(cor_igraph))

  old_membership_values <- membership_df$old_membership[match(node_names, membership_df$node)]

  igraph::V(cor_igraph)$membership <- old_membership_values
  
  subgraphs1 <- lapply(names(t_modules_1), function(i) {
    modelcut(igraph::subgraph(cor_igraph, cfg_t1$membership == i))
  })
  sum_cor <- lapply(seq(1,length(subgraphs1)), function(i) {
    if(diffmessage %in% c("diff","multi")){
      cor_case <- igraph::E(subgraphs1[[i]])$cor_case
      cor_control <- igraph::E(subgraphs1[[i]])$cor_control
      sum(!is.na(cor_case) | !is.na(cor_control))
    }else{
      sum(!is.na(igraph::E(subgraphs1[[i]])$cor))
    }
   
  })
  sum_cor<-unlist(sum_cor)
  sum_node <- lapply(seq(1,length(subgraphs1)), function(i) {
    sum(length(names(igraph::V(subgraphs1[[i]]))))
  })
  sum_node<-unlist(sum_node)
  non_zero_positions <- which(sum_cor != 0 & sum_node>4)
  if(length(non_zero_positions)==0){
    return(NULL)
    message("No subnetwork containing at least 5 nodes, the network clustering analysis is terminated.")
  }else{
  subgraphs1 <- subgraphs1[non_zero_positions]
  t_modules_top <- t_modules_1[non_zero_positions]
  cor_igraph_top<-cor_igraph
  membership_df<-membership_df |>
    dplyr::mutate(membership_top=ifelse(old_membership %in% names(t_modules_top),old_membership,"other"))
  bootnet_cluster_edges<-edges
  bootnet_cluster_nodes<-nodes
  membership_mapping<-data.frame(membership_top=as.character(names(t_modules_top)),membership_top_new =as.character(c(1:length(non_zero_positions))))
  membership_df <- membership_df |>
    dplyr::left_join(membership_mapping,by=c("membership_top"="membership_top")) |>
    dplyr::mutate(membership_top_new = dplyr::coalesce(membership_top_new, "other"))
  
  membership_values_top <- membership_df$membership_top_new[match(bootnet_cluster_nodes$node, membership_df$node)]
  bootnet_cluster_nodes$membership<-paste0("subnet_",membership_values_top)

  names(subgraphs1)<-paste0("subnet_",membership_mapping$membership_top_new)
  subnet_bootnet_list <- lapply(subgraphs1, function(subgraphs) {
    sub_edge<-igrah_to_edgedata(cor_ppi_igraph_f=subgraphs,diffmessage=diffmessage)
    sub_node<-data.frame(node = union(sub_edge$from,sub_edge$to)) |>
      dplyr::left_join(bootnet_cluster_nodes,by=c("node"="node"))
    x<-sub_node$node
    if(diffmessage %in% c("diff","multi")){
      case_edges<-sub_edge |>
        dplyr::mutate(cor=cor_case)
      control_edges<-sub_edge |>
        dplyr::mutate(cor=cor_control)
      case_igraph <- run_igraph(nodes=sub_node,edges=case_edges)
      control_igraph <- run_igraph(nodes=sub_node,edges=control_edges)
      subplot_layout<- qgraph::averageLayout(case_igraph,control_igraph)
      list(nodes = sub_node,edges = sub_edge,plot_layout=subplot_layout)
    }else{
      list(nodes = sub_node,edges = sub_edge)
    }
  })

  unique_classes<-as.factor(paste0("subnet_",unique(membership_values_top)))
  unique_classes <- unique_classes[!(unique_classes %in% "subnet_other")]
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  if(length(unique_classes)>length(cluster_Palette)){
    cluster_Palette <-grDevices::colorRampPalette(cluster_Palette)(length(unique_classes))
  }
  color_mapping_cluster <- setNames(c(cluster_Palette[1:length(unique_classes)],"grey"), c(as.character(unique_classes),"subnet_other"))
  SubNetwork<-new("SubNetwork",
                  group_name=precor_group_name,
                  modularity=modularity,
                  colormapping = color_mapping_cluster,
                  overall_cluster_network = list(
                    nodes=bootnet_cluster_nodes,
                    edges=bootnet_cluster_edges
                  ),
                  subnetworks = subnet_bootnet_list
  )
  return(SubNetwork)
 }
  
} 

#' Convert an igraph Object to Edge Data Frame.
#'
#' Converts an igraph object to a data frame of edges with optional additional columns depending on the diffmessage parameter.
#' This function extracts edge attributes such as correlation values, p-values, etc., into a structured format.
#'
#' @param cor_ppi_igraph_f An igraph object representing the network.
#' @param diffmessage Optional parameter indicating differential message type ("NULL", "diff", or "multi").
#'
#' @return A data frame containing edge information with columns like 'from', 'to', 'cor', etc.
#'
#' @import igraph
#' @export

igrah_to_edgedata<-function(cor_ppi_igraph_f,diffmessage="NULL"){
  edges <- igraph::as_edgelist(cor_ppi_igraph_f, names = TRUE)
  if(diffmessage=="diff"){
    edge_data <- data.frame(
      from = edges[, 1],
      to = edges[, 2],
      cor = igraph::E(cor_ppi_igraph_f)$cor,
      cor_case = igraph::E(cor_ppi_igraph_f)$cor_case,
      cor_control = igraph::E(cor_ppi_igraph_f)$cor_control,
      from_to = igraph::E(cor_ppi_igraph_f)$from_to,
      cor_FC = igraph::E(cor_ppi_igraph_f)$cor_FC,
      cor_p_value = igraph::E(cor_ppi_igraph_f)$cor_p_value,
      cor_status = igraph::E(cor_ppi_igraph_f)$cor_status,
      CIrange_case = igraph::E(cor_ppi_igraph_f)$CIrange_case,
      CIrange_control = igraph::E(cor_ppi_igraph_f)$CIrange_control,
      p_adjust_case = igraph::E(cor_ppi_igraph_f)$p_adjust_case,
      p_adjust_control = igraph::E(cor_ppi_igraph_f)$p_adjust_control
    )
  }else if(diffmessage=="multi"){
    edge_data <- data.frame(
      from = edges[, 1],
      to = edges[, 2],
      cor = igraph::E(cor_ppi_igraph_f)$cor,
      cor_case = igraph::E(cor_ppi_igraph_f)$cor_case,
      cor_control = igraph::E(cor_ppi_igraph_f)$cor_control,
      from_to = igraph::E(cor_ppi_igraph_f)$from_to,
      cor_FC = igraph::E(cor_ppi_igraph_f)$cor_FC,
      cor_p_value = igraph::E(cor_ppi_igraph_f)$cor_p_value,
      cor_status = igraph::E(cor_ppi_igraph_f)$cor_status,
      multiplex_status = igraph::E(cor_ppi_igraph_f)$multiplex_status,
      CIrange_case = igraph::E(cor_ppi_igraph_f)$CIrange_case,
      CIrange_control = igraph::E(cor_ppi_igraph_f)$CIrange_control,
      p_adjust_case = igraph::E(cor_ppi_igraph_f)$p_adjust_case,
      p_adjust_control = igraph::E(cor_ppi_igraph_f)$p_adjust_control
    )
  }else{
    edge_data <- data.frame(
      from = edges[, 1],
      to = edges[, 2],
      cor = igraph::E(cor_ppi_igraph_f)$cor,
      p_adjust= igraph::E(cor_ppi_igraph_f)$p_adjust
    ) 
  }
 
  return(edge_data)
}

#' Filter Network Components Based on Size.
#'
#' Filters out components of a network that are smaller than a specified size (default is 5).
#' Removes isolated vertices and returns the largest connected component(s) of the network.
#'
#' @param net An igraph object representing the network.
#'
#' @return An igraph object representing the filtered network.
#'
#' @import igraph
#' @export

modelcut<-function(net){
  comp <- igraph::components(net)
  if(max(comp$csize)>5){
    net <- igraph::delete_vertices(net, igraph::degree(net) == 0)
    threshold <- 5
    large_comp_indices <- which(comp$csize >= threshold)
    large_comp_nodes <- unlist(comp$membership[comp$membership %in% large_comp_indices])
    g_filtered <- igraph::induced_subgraph(net, names(large_comp_nodes))
    return(g_filtered)
  }else{
    return(net)
  }
}


