#' SubNetwork S4 Class
#'
#' An S4 class that stores subnetwork clustering results from network analysis.
#' Contains modularity information, cluster color mapping, and individual
#' subnetwork structures for visualization and downstream analysis.
#'
#' @slot group_name character, group label identifying the network source
#' @slot modularity numeric, modularity score of the network clustering
#' @slot colormapping ANY, named color vector mapping cluster IDs to colors
#' @slot overall_culster_network ANY, list containing nodes and edges for the
#'                               full clustered network
#' @slot subnetworks ANY, list of individual subnetworks with their nodes,
#'                    edges, and optional layout information
#'
#' @name SubNetwork-class
#' @rdname SubNetwork-class
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a SubNetwork object
#' subnet_obj <- new("SubNetwork",
#'   group_name = "Treatment",
#'   modularity = 0.35,
#'   colormapping = c("subnet_1" = "#FF6B6B", "subnet_2" = "#4ECDC4"),
#'   overall_culster_network = list(nodes = node_df, edges = edge_df),
#'   subnetworks = subnet_list
#' )
#' 
#' # Access components
#' subnet_obj@modularity
#' subnet_obj@colormapping
#' }
setClass("SubNetwork",
         slots = c(
           group_name= "character",
           modularity="numeric",
           colormapping = "ANY",
          # Conditional_network_layout="ANY",
           overall_culster_network = "ANY",
           subnetworks = "ANY"
         )
)

#' Perform subnetwork clustering on a correlation network
#'
#' This function applies community detection algorithms to identify subnetworks
#' (clusters) within a correlation network. It iteratively adjusts resolution
#' to ensure no cluster exceeds the specified maximum size.
#'
#' @param nodes data.frame, node information containing at least a 'node' column
#'              with node identifiers.
#' @param edges data.frame, edge information containing 'from' and 'to' columns
#'              for network edges.
#' @param clustersize integer, maximum allowed size for any single cluster.
#'                    Default 25.
#'
#' @return A community object from igraph containing clustering results with
#'         membership assignments for each node.
#'
#' @details
#' Clustering algorithm:
#' 1. Convert nodes and edges to igraph object using run_igraph.
#' 2. Attempt fast greedy clustering initially.
#' 3. If any cluster exceeds clustersize, switch to Louvain algorithm.
#' 4. Incrementally increase resolution parameter until all clusters are
#'    ≤ clustersize.
#' 5. Return final community structure.
#'
#' @importFrom igraph cluster_fast_greedy as.undirected membership cluster_louvain
#'
#' @examples
#' \dontrun{
#' # Cluster network with maximum cluster size of 20
#' clusters <- run_subnet_cluster(
#'   nodes = node_data,
#'   edges = edge_data,
#'   clustersize = 20
#' )
#' 
#' # View cluster membership
#' table(igraph::membership(clusters))
#' }
#'
#' @export
run_subnet_cluster<-function(nodes=NULL,edges=NULL,clustersize=25){
  cor_igraph<-run_igraph(nodes,edges)#network_show.R
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


#' Add cluster information to network and extract subnetworks
#'
#' This function processes clustering results to create a structured SubNetwork
#' object. It filters out small or empty subnetworks, assigns colors to clusters,
#' and generates individual subnetwork objects with their associated data.
#'
#' @param cfg_t1 igraph community object, clustering results from run_subnet_cluster
#' @param nodes data.frame, node information containing at least a 'node' column
#' @param edges data.frame, edge information with 'from', 'to', and correlation columns
#' @param group_name character, group label(s) for identification. Multiple groups
#'                   are joined with "-vs-".
#' @param diffmessage character, type of differential analysis context:
#'                    "diff" for differential networks, "multi" for multiplex networks,
#'                    or "NULL" for standard networks.
#'
#' @return A SubNetwork S4 object containing:
#'   \item{group_name}{Processed group label}
#'   \item{modularity}{Network modularity score}
#'   \item{colormapping}{Named color vector for cluster visualization}
#'   \item{overall_culster_network}{List with nodes and edges for full clustered network}
#'   \item{subnetworks}{List of individual subnetwork objects, each containing:
#'                      nodes, edges, and optionally plot_layout for differential networks}
#'
#' @details
#' Processing workflow:
#' 1. Calculate modularity from clustering results.
#' 2. Build igraph network and assign membership to nodes.
#' 3. Extract individual subgraphs for each cluster.
#' 4. Filter out clusters with insufficient edges (< 1) or nodes (≤ 4).
#' 5. Map original cluster IDs to sequential numbering (subnet_1, subnet_2, etc.).
#' 6. Generate color palette for clusters using RColorBrewer.
#' 7. Create overall clustered network with membership annotations.
#' 8. For each subnetwork, extract node and edge data, optionally including
#'    layout information for differential/multiplex comparisons.
#' 9. Return structured SubNetwork object.
#'
#' @importFrom dplyr mutate left_join coalesce filter
#' @importFrom igraph V E modularity membership as.undirected subgraph
#' @importFrom qgraph averageLayout
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' \dontrun{
#' # Add cluster information to a standard network
#' subnet_obj <- run_add_cluster(
#'   cfg_t1 = cluster_results,
#'   nodes = node_data,
#'   edges = edge_data,
#'   group_name = "Treatment",
#'   diffmessage = "NULL"
#' )
#' 
#' # Access individual subnetworks
#' first_subnet <- subnet_obj@subnetworks[[1]]
#' first_subnet$nodes
#' first_subnet$edges
#' }
#'
#' @export
run_add_cluster<-function(cfg_t1=NULL,nodes=NULL,edges=NULL,group_name=NULL,
                          diffmessage="NULL"
                       
){
 
  if(length(group_name)>1){
    precor_group_name<-paste(group_name, collapse = "-vs-")
  }else{
    precor_group_name<-gsub(":","-vs-", group_name)#diff
  }
 
  modularity <- igraph::modularity(cfg_t1) 
  ########################membership
  cor_igraph<-run_igraph(nodes,edges)
  t_modules_1 <- sort(table(igraph::membership(cfg_t1)),decr=T)#table(igraph::membership(cfg_t1))#
  t_modules<-table(igraph::membership(cfg_t1))
 
  membership_df <- data.frame(
    old_membership =as.character(igraph::membership(cfg_t1)),#membership_arrange,
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
  # membership_values_top <- membership_df$membership_top[match(node_names, membership_df$node)]
  # V(cor_igraph_top)$membership <- membership_values_top
  ######################################overall_culster_network
  bootnet_cluster_edges<-edges
  bootnet_cluster_nodes<-nodes
  membership_mapping<-data.frame(membership_top=as.character(names(t_modules_top)),membership_top_new =as.character(c(1:length(non_zero_positions))))
  membership_df <- membership_df |>
    dplyr::left_join(membership_mapping,by=c("membership_top"="membership_top")) |>
    dplyr::mutate(membership_top_new = dplyr::coalesce(membership_top_new, "other"))
  
  membership_values_top <- membership_df$membership_top_new[match(bootnet_cluster_nodes$node, membership_df$node)]
  bootnet_cluster_nodes$membership<-paste0("subnet_",membership_values_top)
  #######################################subnetworks
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
  ########################color_mapping_cluster
  unique_classes<-as.factor(paste0("subnet_",unique(membership_values_top)))
  unique_classes <- unique_classes[!(unique_classes %in% "subnet_other")]
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  if(length(unique_classes)>length(cluster_Palette)){
    cluster_Palette <-grDevices::colorRampPalette(cluster_Palette)(length(unique_classes))#避免颜色不够
  }
  color_mapping_cluster <- setNames(cluster_Palette[1:length(unique_classes)], as.character(unique_classes))
  bootnet_cluster_nodes<-bootnet_cluster_nodes |>
     dplyr::filter(membership != "subnet_other")
  bootnet_cluster_edges<-bootnet_cluster_edges |>
      dplyr::filter(from %in% bootnet_cluster_nodes$node) |>
      dplyr::filter(to %in% bootnet_cluster_nodes$node)

  ##################################
  SubNetwork<-new("SubNetwork",
                  group_name=precor_group_name,
                  modularity=modularity,
                  colormapping = color_mapping_cluster,
                  #Conditional_network_layout=Conditional_network_layout,
                  overall_culster_network = list(
                    nodes=bootnet_cluster_nodes,
                    edges=bootnet_cluster_edges
                  ),
                  subnetworks = subnet_bootnet_list
                  
  )
  return(SubNetwork)
 }
  
} 

#' Convert igraph object to edge data frame
#'
#' Internal function that extracts edge attributes from an igraph object and
#' converts them to a structured data frame. Handles different network types
#' including standard, differential, and multiplex networks.
#'
#' @param cor_ppi_igraph_f igraph object, network to convert
#' @param diffmessage character, type of network:
#'                    "diff" for differential networks,
#'                    "multi" for multiplex networks,
#'                    "NULL" for standard networks.
#'
#' @return A data.frame of edges with columns appropriate to network type:
#'   \item{Standard}{from, to, cor, p_adjust, CIrange}
#'   \item{Differential}{from, to, cor, cor_case, cor_control, from_to, cor_FC,
#'                       cor_p_value, cor_status, CIrange_case, CIrange_control,
#'                       p_adjust_case, p_adjust_control}
#'   \item{Multiplex}{Same as differential plus multiplex_status}
#'
#' @details
#' Extracted attributes include:
#' \itemize{
#'   \item from, to: Edge endpoints
#'   \item cor: Correlation coefficient(s)
#'   \item p_adjust: Adjusted p-values
#'   \item CIrange: Confidence interval range
#'   \item For differential/multiplex: Case/control specific values
#'   \item from_to: Combined node identifier
#'   \item cor_FC: Fold change of correlations
#'   \item cor_status: Edge status (present/absent/divergent)
#' }
#'
#' @keywords internal
#'
#' @importFrom igraph as_edgelist E
#'
#' @noRd
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
  }else if(diffmessage=="multi"){   ##################################
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
      p_adjust= igraph::E(cor_ppi_igraph_f)$p_adjust,      
      CIrange=igraph::E(cor_ppi_igraph_f)$CIrange
    ) 
  }
 
  return(edge_data)
}

#' Filter network to retain only connected components above size threshold
#'
#' Internal function that removes isolated nodes and filters network components
#' to only those with at least 5 nodes. Ensures subnetworks have sufficient
#' connectivity for meaningful analysis.
#'
#' @param net igraph object, network to filter
#'
#' @return An igraph object containing only:
#'   \itemize{
#'     \item Nodes with degree > 0 (non-isolated)
#'     \item Connected components with size ≥ 5
#'   }
#'
#' @details
#' Processing steps:
#' 1. Identify connected components in the network.
#' 2. Remove nodes with degree 0.
#' 3. Filter to keep only components with ≥ 5 nodes.
#' 4. Return induced subgraph of retained nodes.
#'
#' If no component meets the size threshold, the original network is returned
#' (though likely small or empty).
#'
#' @keywords internal
#'
#' @importFrom igraph components degree delete_vertices induced_subgraph
#'
#' @noRd
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


