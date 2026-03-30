#' Network_show S4 Class
#'
#' S4 class for storing network visualization results including data, igraph objects,
#' and ggplot2 plots.
#'
#' @slot data list, contains network data including nodes, edges, and layout information
#' @slot igraph ANY, igraph object(s) representing the network(s)
#' @slot plot ANY, ggplot2 plot object(s) for visualization
#'
#' @name Network_show-class
#' @rdname Network_show-class
#' @export
setClass("Network_show", slots = c(
  data="list",
  igraph = "ANY",
  plot="ANY"
))

#' Enrichment_show S4 Class
#'
#' S4 class for storing functional enrichment results including data and plots.
#'
#' @slot data list, contains enrichment analysis results
#' @slot plot ANY, ggplot2 plot object(s) for enrichment visualization
#'
#' @name Enrichment_show-class
#' @rdname Enrichment_show-class
#' @export
setClass("Enrichment_show", slots = c(
  data="list",
  plot="ANY"
))

# 定义 show 方法
setMethod("show", "Network_show", function(object) {
  cat("=== Network_show ===\n")
  cat("- Data: object@data\n")
  cat("- Network (igraph): object@igraph\n")
  cat("- Plot: object@plot\n")
  print(object@plot)
})
setMethod("show", "Enrichment_show", function(object) {
  cat("=== Enrichment_show ===\n")
  cat("- Data: object@data\n")
  cat("- Plot: object@plot\n")
  print(object@plot)
})



#' Visualize networks from various network analysis objects
#'
#' This comprehensive visualization function generates network plots for different
#' network types including stability tests, overall networks, clustered networks,
#' subnetworks, differential networks, and multiplex networks. It supports
#' various layout algorithms, node coloring schemes, edge styling, and interactive
#' visualizations.
#'
#' @param Network S4 object, network result object (StableNetwork, Conditional_network,
#'                SubNetwork, Differential_network, etc.)
#' @param plot_type character, type of network to plot. Options:
#'   \itemize{
#'     \item "stable_test", "case_stable_test", "control_stable_test": Bootstrap stability networks
#'     \item "overall_network", "case_overall_network", "control_overall_network": Overall correlation networks
#'     \item "overall_culster_network": Clustered overall network
#'     \item "sub_network": Individual subnetworks
#'     \item "diff_network": Differential network
#'     \item "diff_overall_cluster_network": Clustered differential network
#'     \item "diff_subnetwork": Differential subnetworks
#'     \item "differential_network": Combined case/control/differential networks
#'     \item "differential_subnetwork": Combined case/control/differential subnetworks
#'     \item "interaction_network": Interaction/multiplex networks
#'     \item "case_multi_network", "control_multi_network": Multiplex networks by group
#'   }
#' @param stable_num integer, number of stable networks to display in stability tests.
#'                   Default 9.
#' @param subnetwork_name character vector, names of subnetworks to display.
#'                        "all" displays all subnetworks. Default "all".
#' @param layout_type character, layout algorithm for network visualization.
#'                     Options: "fr" (Fruchterman-Reingold), "kk" (Kamada-Kawai),
#'                     "nicely" (automatic). Default "fr".
#' @param input_layout matrix, optional custom layout matrix for network coordinates.
#' @param richfactor_threshold numeric, threshold for enrichment factor filtering.
#'                              Default 0.
#' @param node_colortype character, node coloring scheme. Options:
#'                        "Class" (by omics type), "membership" (by cluster),
#'                        "FC", "Log2FC", "Normalized mean". Default "Class".
#' @param R_threshold numeric, correlation threshold for edge filtering. Default 0.
#' @param focus character vector, focus on specific edge types for differential networks.
#'               Default "all".
#' @param node_name_size numeric, font size for node labels. Default 2.
#' @param node_name_type character, node label type. "id" for feature ID,
#'                       "name" for feature name. Default "id".
#' @param node_size numeric, fixed node size. If NULL, size is scaled by centrality.
#'                  Default NULL.
#' @param image_margin_size numeric, margin around network plot. Default 0.5.
#' @param show_node_name logical, whether to display node labels. Default FALSE.
#' @param show_node_legend logical, whether to display node legend. Default FALSE.
#' @param show_edge_legend logical, whether to display edge legend. Default FALSE.
#' @param plot_title_size numeric, font size for plot titles. Default 12.
#' @param axis_title_size numeric, font size for axis titles. Default 8.
#' @param text_size numeric, general text size. Default 8.
#' @param legend_title_size numeric, font size for legend titles. Default 8.
#' @param legend_text_size numeric, font size for legend text. Default 8.
#' @param font_family character, font family for text. Default "Arial".
#' @param add_enrichement logical, whether to add enrichment annotations. Default FALSE.
#' @param add_Centrality character vector, centrality measures to add for node sizing.
#'                       Options: "betweenness", "degree", "eigenvector". Default NULL.
#' @param interactive logical, whether to generate interactive plots with ggiraph.
#'                    Default FALSE.
#' @param outdir character, output directory for saving results. Default "./".
#' @param centrality_scatterplot logical, whether to generate centrality scatterplots.
#'                                Default TRUE.
#' @param datasave logical, whether to save node/edge data to files. Default TRUE.
#' @param alpha numeric, transparency for edges. Default 1.
#'
#' @return A Network_show S4 object containing:
#'   \item{data}{List of network data including nodes, edges, and layout coordinates}
#'   \item{igraph}{List of igraph objects}
#'   \item{plot}{Combined ggplot2 plot(s)}
#'
#' @details
#' The function handles multiple network types:
#'
#' **Stability Tests**:
#' - Displays multiple bootstrap networks
#' - Uses average layout for consistent positioning
#'
#' **Overall Networks**:
#' - Single network visualization with correlation edges
#' - Supports node coloring by class, FC, or normalized expression
#'
#' **Subnetworks**:
#' - Individual subnetwork plots with consistent layouts
#' - Optional centrality measures for node sizing
#'
#' **Differential Networks**:
#' - Combined visualization of case, control, and differential networks
#' - Edge coloring by correlation status (enhanced/only/conflict)
#'
#' **Layout Options**:
#' - "fr": Fruchterman-Reingold force-directed layout
#' - "kk": Kamada-Kawai spring-embedded layout
#' - "nicely": Automatic layout selection
#'
#' @importFrom dplyr mutate filter arrange
#' @importFrom igraph graph_from_data_frame V E degree betweenness eigen_centrality
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggplot2 scale_fill_manual scale_fill_gradientn scale_edge_color_manual
#' @importFrom patchwork wrap_plots
#' @importFrom ggiraph girafe
#' @importFrom qgraph averageLayout
#'
#' @examples
#' \dontrun{
#' # Plot overall network with class-based node coloring
#' net_show <- network_show(
#'   Network = stable_network,
#'   plot_type = "overall_network",
#'   show_node_name = TRUE,
#'   show_node_legend = TRUE
#' )
#' 
#' # Plot differential network with focus on specific edge types
#' diff_show <- network_show(
#'   Network = differential_network,
#'   plot_type = "differential_network",
#'   focus = c("Only in Treatment", "Enhanced in Treatment")
#' )
#' 
#' # Plot top 5 subnetworks with centrality-based node sizing
#' sub_show <- network_show(
#'   Network = subnetwork_object,
#'   plot_type = "sub_network",
#'   subnetwork_name = 1:5,
#'   add_Centrality = "betweenness",
#'   show_node_name = TRUE
#' )
#' }
#'
#' @export
network_show<-function(Network=NULL,plot_type="stable_test",stable_num=9,
                       subnetwork_name=c("all"),layout_type="fr",input_layout=NULL,richfactor_threshold=0,#diff_subnetwork_name=c("all")
                       node_colortype="Class",R_threshold=0,focus=c("all"),node_name_size=2,node_name_type="id",node_size=NULL,
                       image_margin_size=0.5,
                       show_node_name=FALSE,show_node_legend=FALSE,show_edge_legend=FALSE,
                       plot_title_size=12,
                       axis_title_size=8,
                       text_size=8,
                       legend_title_size=8,
                       legend_text_size=8,
                       font_family="Arial",
                       add_enrichement=FALSE,add_Centrality=NULL,interactive=FALSE,outdir="./",
                       centrality_scatterplot=TRUE,datasave=TRUE,
                       alpha=1){ 
  
  # showtext::showtext_auto()
  set.seed(123)
  nodes<-NULL
  LAYOUT<-NULL
  edges<-NULL
  plot<-NULL
  plotlist<-NULL
  is_related_node<-NULL
  group_name=Network@group_name
  split_result <- strsplit(group_name, "-vs-")[[1]]
  if(length(split_result)>1){
    casename <- split_result[1] 
    controlname <- split_result[2] 
  }else{
    casename <- "case" 
    controlname <- "control"
  }
  Only_in_control<-paste0("Only in ",controlname)
  Only_in_case<-paste0("Only in ",casename)
  Enhanced_in_case<-paste0("Enhanced in ",casename)
  Enhanced_in_control<-paste0("Enhanced in ",controlname) 
  case_control<-paste0("(",casename,":",controlname,")") 
  
  edge_status <- c("positive", "negative", "other")
  edgecolormapping <- setNames(c("red", "blue", "gray"), edge_status)
  linetype_map <-setNames(c("longdash","longdash","solid"), c("interaction","cor","interaction_cor"))
  
    if(plot_type %in% c("stable_test","case_stable_test","control_stable_test")){
      if(plot_type =="stable_test"){
        if(is(Network, "Stable_SubNetwork")) {
          net_list<-Network@StableNetwork@bootnet_result@bootnet_list
        }else if(is(Network, "StableNetwork")){
          net_list<-Network@bootnet_result@bootnet_list
        }else{
          stop(paste0("‘stable_test‘ does not applies to the target object :",class(Network)))
        } 
      }
      if(plot_type =="case_stable_test"){
        group_name<-unlist(strsplit(group_name,split = "-vs-"))[1]
        if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork") ) {
          net_list<-Network@Conditional_network@network_case@bootnet_result@bootnet_list
        }else if(is(Network, "Conditional_network")){
          net_list<-Network@network_case@bootnet_result@bootnet_list 
        }else{
          stop(paste0("‘case_stable_test‘ does not applies to the target object :",class(Network)))
        } 
      }
      if(plot_type =="control_stable_test"){
        group_name<-unlist(strsplit(group_name,split = "-vs-"))[2]
        if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork") ) {
          net_list<-Network@Conditional_network@network_control@bootnet_result@bootnet_list
        }else if(is(Network, "Conditional_network")){
          net_list<-Network@network_control@bootnet_result@bootnet_list 
        }else{
          stop(paste0("‘control_stable_test‘ does not applies to the target object :",class(Network)))
        } 
      } 
      nodes<-net_list[[1]]$nodes
      ################################
    
      }else if(plot_type %in% c("overall_network","case_overall_network","control_overall_network")){
      if(plot_type =="overall_network"){
        if(is(Network, "Stable_SubNetwork")) {
          nodes<-Network@StableNetwork@bootnet_result_filter@bootnet_node
          edges<-Network@StableNetwork@bootnet_result_filter@bootnet_edge |>
            dplyr::filter(abs(cor)>R_threshold)
        }else if(is(Network, "StableNetwork")){
          nodes<-Network@bootnet_result_filter@bootnet_node
          edges<-Network@bootnet_result_filter@bootnet_edge |>
            dplyr::filter(abs(cor)>R_threshold)
        }else{
          stop(paste0("‘overall_network‘ does not applies to the target object :",class(Network)))
        }
      }else if (plot_type =="case_overall_network"){
        group_name<-unlist(strsplit(group_name,split = "-vs-"))[1]
        if ((is(Network, "Stable_DifferentialNetwork")) || is(Network, "Stable_MultiplexNetwork")){
          nodes<-Network@Conditional_network@network_case@bootnet_result_filter@bootnet_node |>
            dplyr::mutate(`Normalized mean` = case_norm_mean)
          edges<-Network@Conditional_network@network_case@bootnet_result_filter@bootnet_edge |>
            
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else if(is(Network, "Conditional_network")){
          nodes<-Network@network_case@bootnet_result_filter@bootnet_node |>
            dplyr::mutate(`Normalized mean` = case_norm_mean)
          edges<-Network@network_case@bootnet_result_filter@bootnet_edge |>
            
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else{
          stop(paste0("‘case_overall_network‘ does not applies to the target object :",class(Network)))
        } 
      }else if (plot_type =="control_overall_network"){
        group_name<-unlist(strsplit(group_name,split = "-vs-"))[2]
        if ((is(Network, "Stable_DifferentialNetwork")) || is(Network, "Stable_MultiplexNetwork")){
          nodes<-Network@Conditional_network@network_control@bootnet_result_filter@bootnet_node |>
            dplyr::mutate(`Normalized mean` = control_norm_mean)
          edges<-Network@Conditional_network@network_control@bootnet_result_filter@bootnet_edge |>
            
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else if(is(Network, "Conditional_network")){
          nodes<-Network@network_control@bootnet_result_filter@bootnet_node |>
            dplyr::mutate(`Normalized mean` = control_norm_mean)
          edges<-Network@network_control@bootnet_result_filter@bootnet_edge |>
            
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else{
          stop(paste0("‘control_overall_network‘ does not applies to the target object :",class(Network)))
        } 
      }
      ################################
    
        }else if(plot_type =="overall_culster_network"){
      if(is(Network, "Stable_SubNetwork")) {
        nodes<-Network@SubNetwork@overall_culster_network$nodes
        edges<-Network@SubNetwork@overall_culster_network$edges |>
          dplyr::filter(abs(cor)>R_threshold)
      }else if(is(Network, "SubNetwork")){
        nodes<-Network@overall_culster_network$nodes
        edges<-Network@overall_culster_network$edges |>
          dplyr::filter(abs(cor)>R_threshold)
      }else{
        stop(paste0("‘overall_culster_network‘ does not applies to the target object :",class(Network)))
      }
      ####################################
    
          }else if(plot_type =="sub_network"){
      if(is(Network, "Stable_SubNetwork")) {
        net_list<-Network@SubNetwork@subnetworks
      }else if(is(Network, "SubNetwork")){
        net_list<-Network@subnetworks
      }else{
        stop(paste0("‘sub_network‘ does not applies to the target object :",class(Network)))
      }
      ################################
    
            }else if(plot_type =="interaction_network"){
      if(is(Network, "Stable_MultiplexNetwork")){
        # LAYOUT<-Network@Interaction_network_layout
        nodes<-Network@Interaction_network@nodes
        edges<-Network@Interaction_network@edges 
      }else if(is(Network, "Interaction_network")){
        nodes<-Network@nodes
        edges<-Network@edges
      }else{
        stop(paste0("‘interaction_network‘ does not applies to the target object :",class(Network)))
      }
      ################################
      
    
              }else if(plot_type %in% c("case_multi_network","control_multi_network")){
      #LAYOUT<-Network@Interaction_network_layout
      # LAYOUT<-"kk"#Network@Interaction_network_layout
      if(plot_type =="case_multi_network"){
        group_name<-unlist(strsplit(group_name,split = "-vs-"))[1]
        if(is(Network, "Stable_MultiplexNetwork")){ 
          nodes<-Network@Conditional_multiplexnetwork@network_case@nodes |>
            dplyr::mutate(`Normalized mean` = case_norm_mean)
          edges<-Network@Conditional_multiplexnetwork@network_case@edges |>
            
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else if(is(Network, "Conditional_multiplexnetwork")){
          nodes<-Network@network_case@nodes |>
            dplyr::mutate(`Normalized mean` = case_norm_mean)
          edges<-Network@network_case@edges |>
            
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else{
          stop(paste0("‘case_multi_network‘ does not applies to the target object :",class(Network)))
        }
      }else if(plot_type =="control_multi_network"){
        group_name<-unlist(strsplit(group_name,split = "-vs-"))[2]
        if(is(Network, "Stable_MultiplexNetwork")){ 
          nodes<-Network@Conditional_multiplexnetwork@network_control@nodes |>
            dplyr::mutate(`Normalized mean` = control_norm_mean)
          edges<-Network@Conditional_multiplexnetwork@network_control@edges |>
            
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else if(is(Network, "Conditional_multiplexnetwork")){
          nodes<-Network@network_control@nodes |>
            dplyr::mutate(`Normalized mean` = control_norm_mean)
          edges<-Network@network_control@edges |>
            
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else{
          stop(paste0("‘control_multi_network‘ does not applies to the target object :",class(Network)))
        }
      }
      
      ################################
    
                }else if(plot_type %in% c("diff_network")){
      if(is(Network, "Stable_DifferentialNetwork")){
        nodes<-Network@Differential_network@diff_nodes
        edges<-Network@Differential_network@diff_edges 
      }else if(is(Network, "Differential_network")){
        nodes<-Network@diff_nodes
        edges<-Network@diff_edges 
      }else if(is(Network, "Stable_MultiplexNetwork")){
        #LAYOUT<-Network@Interaction_network_layout
        nodes<-Network@Differential_multiplexnetwork@diff_nodes
        edges<-Network@Differential_multiplexnetwork@diff_edges
      }else if(is(Network, "Differential_multiplexnetwork")){
        nodes<-Network@diff_nodes
        edges<-Network@diff_edges
      }else{
        stop(paste0("‘diff_network‘ does not applies to the target object :",class(Network)))
      }
    }else if(plot_type %in% c("diff_overall_cluster_network")){
      if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork")){
        nodes<-Network@Differential_subnetwork@overall_culster_network$nodes
        edges<-Network@Differential_subnetwork@overall_culster_network$edges
      }else if(is(Network, "SubNetwork")){
        nodes<-Network@overall_culster_network$nodes
        edges<-Network@overall_culster_network$edges
      }else{
        stop(paste0("‘diff_overall_cluster_network‘ does not applies to the target object :",class(Network)))
      }
    }else if(plot_type %in% c("diff_subnetwork")){
      if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork")){
        net_list<-Network@Differential_subnetwork@subnetworks
      }else if(is(Network, "SubNetwork")){
        net_list<-Network@subnetworks
       
      }else{
        stop(paste0("‘diff_subnetwork‘ does not applies to the target object :",class(Network)))
      }
    
      
      }else if(plot_type %in% c("case_subnetwork","control_subnetwork")){
      if(plot_type %in% c("case_subnetwork")){
        group_name<-unlist(strsplit(group_name,split = "-vs-"))[1]
        if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork")){
          net_list<-Network@Differential_subnetwork@subnetworks
          # net_name<-names(Network@Differential_subnetwork@subnetworks)
          # nodes<-Network@Differential_subnetwork@subnetworks[[net_name]]$nodes
          # edges<-Network@Differential_subnetwork@subnetworks[[net_name]]$edges |>
          #   mutate(cor=cor_case)
        }else if(is(Network, "SubNetwork")){
          net_list<-Network@subnetworks
          # net_name<-names(Network@subnetworks)
          # nodes<-Network@subnetworks[[net_name]]$nodes
          # edges<-Network@subnetworks[[net_name]]$edges |>
          #   mutate(cor=cor_case)
        }else{
          stop(paste0("‘case_subnetwork‘ does not applies to the target object :",class(Network)))
        }
        net_list <- lapply(net_list, function(subnet) {
          subnet$edges <- subnet$edges |>
            dplyr::mutate(cor = cor_case)
          subnet$nodes<-subnet$nodes |>
            dplyr::mutate(`Normalized mean` = case_norm_mean)
          subnet
        })
      }else if(plot_type %in% c("control_subnetwork")){
        group_name<-unlist(strsplit(group_name,split = "-vs-"))[2]
        if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork")){
          net_list<-Network@Differential_subnetwork@subnetworks
          # net_name<-names(Network@Differential_subnetwork@subnetworks)
          # nodes<-Network@Differential_subnetwork@subnetworks[[net_name]]$nodes
          # edges<-Network@Differential_subnetwork@subnetworks[[net_name]]$edges |>
          #   mutate(cor=cor_control)
        }else if(is(Network, "SubNetwork")){
          net_list<-Network@subnetworks
          # net_name<-names(Network@subnetworks)
          # nodes<-Network@subnetworks[[net_name]]$nodes
          # edges<-Network@subnetworks[[net_name]]$edges |>
          #   mutate(cor=cor_control)
        }else{
          stop(paste0("‘control_subnetwork‘ does not applies to the target object :",class(Network)))
        }
        net_list <- lapply(net_list, function(subnet) {
          subnet$edges <- subnet$edges |>
            dplyr::mutate(cor = cor_control)
          subnet$nodes<-subnet$nodes |>
            dplyr::mutate(`Normalized mean` = control_norm_mean)
          subnet
        })
      }
      #####################test
    }else if(plot_type %in% c("differential_network")){
      if ((is(Network, "Stable_DifferentialNetwork"))){
        case_nodes<-Network@Conditional_network@network_case@bootnet_result_filter@bootnet_node |>
          dplyr::mutate(`Normalized mean` = case_norm_mean)
        case_edges<-Network@Conditional_network@network_case@bootnet_result_filter@bootnet_edge |>
          
          dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        control_nodes<-Network@Conditional_network@network_control@bootnet_result_filter@bootnet_node |>
          dplyr::mutate(`Normalized mean` = control_norm_mean)
        control_edges<-Network@Conditional_network@network_control@bootnet_result_filter@bootnet_edge |>
          
          dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        nodes<-Network@Differential_network@diff_nodes
        edges<-Network@Differential_network@diff_edges 
      }else if(is(Network, "Stable_MultiplexNetwork")){ 
        case_nodes<-Network@Conditional_multiplexnetwork@network_case@nodes |>
          dplyr::mutate(`Normalized mean` = case_norm_mean)
        case_edges<-Network@Conditional_multiplexnetwork@network_case@edges |>
          
          dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        control_nodes<-Network@Conditional_multiplexnetwork@network_control@nodes |>
          dplyr::mutate(`Normalized mean` = control_norm_mean)
        control_edges<-Network@Conditional_multiplexnetwork@network_control@edges |>
          
          dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor)) 
        nodes<-Network@Differential_multiplexnetwork@diff_nodes
        edges<-Network@Differential_multiplexnetwork@diff_edges
      }else{
        stop(paste0("‘differential_network‘ does not applies to the target object :",class(Network)))
      }
    }else if(plot_type %in% c("differential_subnetwork")){
      if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork")){
        casecon_net_list<-Network@Differential_subnetwork@subnetworks
        case_net_list <- lapply(casecon_net_list, function(subnet) {
          subnet$edges <- subnet$edges |>
            dplyr::mutate(cor = cor_case)
          subnet$nodes<-subnet$nodes |>
            dplyr::mutate(`Normalized mean` = case_norm_mean)
          subnet
        })
        control_net_list <- lapply(casecon_net_list, function(subnet) {
          subnet$edges <- subnet$edges |>
            dplyr::mutate(cor = cor_control)
          subnet$nodes<-subnet$nodes |>
            dplyr::mutate(`Normalized mean` = control_norm_mean)
          subnet
        })
      }else{
        stop(paste0("‘differential_subnetwork‘ does not applies to the target object :",class(Network)))
      }
    }

    ##################################################### 
    
    if(plot_type %in% c("interaction_network","overall_network","overall_culster_network",
                        "case_multi_network","control_multi_network",
                        "case_overall_network","control_overall_network",
                        "diff_network","diff_overall_cluster_network")){
      igraph <- run_igraph(nodes=nodes,edges=edges)
      if(layout_type=="kk"){
        LAYOUT<-igraph::layout_with_kk(igraph)
      }else if(layout_type=="nicely"){
        LAYOUT<-igraph::layout_nicely(igraph)
      }else if(layout_type=="fr"){
        LAYOUT<-igraph::layout_with_fr(igraph)
      }else{
        stop("The layout_type parameter must be either ‘kk’, ‘nicely’, or ‘fr’.")
      }
      
      igraphlist<-list(net=igraph)
      net_name<-"net"
      names(igraphlist)=net_name
    }else if(plot_type %in% c("stable_test","case_stable_test","control_stable_test")){
      if(stable_num<=length(net_list)){
        net_name<- paste0("boot ",c(1:stable_num))
      }else{
        net_name<-paste0("boot ",c(1:length(net_list)))
      }
      igraphlist<-lapply(net_name,function(x){
        nodes<-net_list[[x]]$nodes
        edges<-net_list[[x]]$edges |>
          dplyr::filter(abs(cor)>R_threshold)
        run_igraph(nodes=nodes,edges=edges)
      }) 
      names(igraphlist)=net_name
      #####
      LAYOUT<-qgraph::averageLayout(igraphlist)
      net_list<-lapply(net_name, function(subnet) {
        net_list[[subnet]]$plot_layout <- LAYOUT
        net_list[[subnet]]
      })
      names(net_list)=net_name
    }else if(plot_type %in% c("sub_network")){
      if("all" %in% subnetwork_name){
        net_name<-names(net_list)
      }else{
        #  net_name<-paste0("subnet_",subnetwork_name) 
        net_name<- ifelse(grepl("^subnet_", subnetwork_name), subnetwork_name, paste0("subnet_", subnetwork_name))#paste0("subnet_",subnetwork_name) 
        
      }
      igraphlist <-lapply(net_name, function(x) {
        nodes<-net_list[[x]]$nodes
        edges<-net_list[[x]]$edges |>
          dplyr::filter(abs(cor)>R_threshold)
        igraph <- run_igraph(nodes=nodes,edges=edges)
        if(!is.null(add_Centrality)){
          igraph<-run_add_Centrality(igraph=igraph,add_Centrality=add_Centrality)
          if(centrality_scatterplot){
            Centralitylist<-run_plot_Centrality(igraph=igraph,add_Centrality=add_Centrality,
                                                nodes=nodes,node_name_type=node_name_type)
            Centrality_data<-Centralitylist[["Centrality_data"]]
            Centrality_plot<-Centralitylist[["Centrality_plot"]]
            Centrality_plotname<-paste0("scatterplot_", stringr::str_c(add_Centrality, collapse = "_"),"_",
                                        x,"_",group_name)
            # readr::write_delim(Centrality_data,file.path(outdir,paste0(Centrality_plotname,".txt")),delim="\t")
            ggplot2::ggsave(file.path(outdir,paste0(Centrality_plotname,".png")),
                            plot = Centrality_plot, width = 6, height = 6, dpi = 300,bg = "white")
          }
        }
        igraph
      })
      names(igraphlist)=net_name
      LAYOUT_list<-lapply(net_name, function(x) { 
        igraph<-igraphlist[[x]]
        if(layout_type=="kk"){
          igraph::layout_with_kk(igraph)
        }else if(layout_type=="nicely"){
          igraph::layout_nicely(igraph)
        }else if(layout_type=="fr"){
          igraph::layout_with_fr(igraph)
        }else{
          stop("The layout_type parameter must be either ‘kk’, ‘nicely’, or ‘fr’.")
        }
      })
      names(LAYOUT_list)=net_name
      net_list<-lapply(net_name, function(subnet) {
        net_list[[subnet]]$plot_layout <- LAYOUT_list[[subnet]]
        net_list[[subnet]]
      })
      names(net_list)=net_name
      
    }else if(plot_type %in% c("diff_subnetwork","control_subnetwork","case_subnetwork")){
      if("all" %in% subnetwork_name){
        net_name<-names(net_list)
      }else{
        #net_name<-paste0("subnet_",subnetwork_name) 
        net_name<- ifelse(grepl("^subnet_", subnetwork_name), subnetwork_name, paste0("subnet_", subnetwork_name))#paste0("subnet_",subnetwork_name) 
      }
      igraphlist <-lapply(net_name, function(x) {
        nodes<-net_list[[x]]$nodes
        edges<-net_list[[x]]$edges 
        run_igraph(nodes=nodes,edges=edges)
        
      })
      names(igraphlist)=net_name
      LAYOUT_list<-lapply(net_name, function(x) {
        net_list[[x]]$plot_layout
      })
      names(LAYOUT_list)=net_name
     #######test 
    }else if(plot_type %in% c("differential_network")){
      case_igraph <- run_igraph(nodes=case_nodes,edges=case_edges)
      control_igraph <- run_igraph(nodes=control_nodes,edges=control_edges)
      igraph <- run_igraph(nodes=nodes,edges=edges)
      if(layout_type=="kk"){
        LAYOUT<-igraph::layout_with_kk(case_igraph)
      }else if(layout_type=="nicely"){
        LAYOUT<-igraph::layout_nicely(case_igraph)
      }else if(layout_type=="fr"){
        LAYOUT<-igraph::layout_with_fr(case_igraph)
      }else{
        stop("The layout_type parameter must be either ‘kk’, ‘nicely’, or ‘fr’.")
      }
      casename=unlist(strsplit(group_name,split = "-vs-"))[1]
      controlname=unlist(strsplit(group_name,split = "-vs-"))[2]
      igraphlist<-list(case_igraph,control_igraph,igraph)
      net_name<-c(casename,controlname,group_name)
      names(igraphlist)=net_name
      net_list<-list(list(nodes=case_nodes,edges=case_edges,plot_layout=LAYOUT),
                     list(nodes=control_nodes,edges=control_edges,plot_layout=LAYOUT),
                     list(nodes=nodes,edges=edges,plot_layout=LAYOUT))
      names(net_list)=net_name
      LAYOUT_list <-list(LAYOUT,LAYOUT,LAYOUT)
      names(LAYOUT_list)=net_name
    }else if(plot_type %in% c("differential_subnetwork")){
      if("all" %in% subnetwork_name){
        net_name<-names(casecon_net_list)
      }else{
        #net_name<-paste0("subnet_",subnetwork_name) 
        net_name<- ifelse(grepl("^subnet_", subnetwork_name), subnetwork_name, paste0("subnet_", subnetwork_name))#paste0("subnet_",subnetwork_name) 
      }
      casecon_igraphlist <-lapply(net_name, function(x) {
        nodes<-casecon_net_list[[x]]$nodes
        edges<-casecon_net_list[[x]]$edges 
        run_igraph(nodes=nodes,edges=edges)
      })
      case_igraphlist <-lapply(net_name, function(x) {
        nodes<-case_net_list[[x]]$nodes
        edges<-case_net_list[[x]]$edges 
        run_igraph(nodes=nodes,edges=edges)
        
      })
      control_igraphlist <-lapply(net_name, function(x) {
        nodes<-control_net_list[[x]]$nodes
        edges<-control_net_list[[x]]$edges 
        run_igraph(nodes=nodes,edges=edges)
        
      })
      igraphlist <-list()
      # Loop through each index and append the sublists in the required order
      for (i in 1:length(net_name)) {
        igraphlist <- append(igraphlist, list(case_igraphlist[[i]], control_igraphlist[[i]], casecon_igraphlist[[i]]))
      }
      casename=unlist(strsplit(group_name,split = "-vs-"))[1]
      controlname=unlist(strsplit(group_name,split = "-vs-"))[2]
      prefixes <- c(casename, controlname, group_name)
      # Generate combinations
      combine_net_name <- character(0)
      for (name in net_name) {
        combine_net_name <- c(combine_net_name, paste0(prefixes, "_", name))
      }
      names(igraphlist)=combine_net_name
      net_list <- list()
      # Loop through each index and append the sublists in the required order
      for (i in net_name) {
        net_list <- append(net_list, list(case_net_list[[i]], control_net_list[[i]], casecon_net_list[[i]]))
      }
      names(net_list)=combine_net_name
      LAYOUT_list0<-lapply(net_name, function(x) {
        casecon_net_list[[x]]$plot_layout
      })
      LAYOUT_list <-list()
      # Loop through each index and append the sublists in the required order
      for (i in 1:length(net_name)) {
        LAYOUT_list <- append(LAYOUT_list, list(LAYOUT_list0[[i]], LAYOUT_list0[[i]], LAYOUT_list0[[i]]))
      }
      names(LAYOUT_list)=combine_net_name
      net_name=combine_net_name

    }
    
    
    ##################################################### 
    if(plot_type %in% c("stable_test","case_stable_test","control_stable_test","interaction_network",
                        "overall_network")){
      color_type="Class"
    }else if(plot_type %in% c("case_overall_network","control_overall_network",
                              "case_multi_network","control_multi_network",
                              "control_subnetwork","case_subnetwork")){
      color_type=node_colortype
      if(!(color_type %in% c("Class","FC","Log2FC","Normalized mean"))){
        stop("The node_colortype parameter of ‘sub_network’ must be either ‘Class’, ‘FC’ , ‘Log2FC’, or ‘Normalized mean’.")
      }
    }else if(plot_type %in% c("overall_culster_network","diff_overall_cluster_network")){
      color_type="membership"
    }else if(plot_type %in% c("sub_network")){
      color_type=node_colortype
      if(!(color_type %in% c("Class","membership"))){
        stop("The node_colortype parameter of ‘sub_network’ must be either ‘Class’ or ‘membership’.")
      }
    }else if(plot_type %in% c("diff_network","diff_subnetwork")){
      color_type=node_colortype
      if(!(color_type %in% c("Class","FC","Log2FC"))){
        stop("The node_colortype parameter of ‘sub_network’ must be either ‘Class’, ‘FC’ ,or ‘Log2FC’.")
      }
    }else if(plot_type %in% c("differential_subnetwork")){
      color_typelist=rep(c("Normalized mean","Normalized mean","Log2FC"),length(case_igraphlist))
      color_type <-NULL
    }else if(plot_type %in% c("differential_network")){
      color_typelist=rep(c("Normalized mean","Normalized mean","Log2FC"),1)
      color_type <-NULL
    }else{
      color_type <-NULL
    }
    #####################################################
    if(is(Network, "Stable_SubNetwork") || is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork") ) {
      color_mapping<-Network@colormapping
    }else{
     # color_mapping<-color_mapping#f1
      color_mapping<-run_color(annotation_table=nodes)
    }
    if( c("membership") %in% color_type){
    if(is(Network, "Stable_SubNetwork")){
      color_mapping_cluster<-Network@SubNetwork@colormapping
    }else if(is(Network, "SubNetwork")){
      color_mapping_cluster<-Network@colormapping
    }else if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork")){
      color_mapping_cluster<-Network@Differential_subnetwork@colormapping
    }else{
      color_mapping_cluster<-NULL
    }
    }
    ##################################################### 
    if(plot_type %in% c("diff_network","diff_overall_cluster_network","diff_subnetwork")){
      edge_color_type="temp_cor_status"
    }else if(plot_type %in% c("differential_subnetwork")){
      edge_color_typelist=rep(c("color","color","temp_cor_status"),length(case_igraphlist))
    }else if(plot_type %in% c("differential_network")){
      edge_color_typelist=rep(c("color","color","temp_cor_status"),1)
    }else{
      edge_color_type="color"
    }
    
    ##################################################### 
    if(plot_type %in% c("case_multi_network","control_multi_network",
                        "diff_network","diff_overall_cluster_network"
    )){
      if(length(unique(igraph::E(igraph)$multiplex_status))>0){
        edgelinetypemap <- linetype_map[unique(igraph::E(igraph)$multiplex_status)]
      }else{
        edgelinetypemap<-NULL
      }
    }else if(plot_type %in% c("control_subnetwork","case_subnetwork","diff_subnetwork")){
      
      edgelinetypemaplist<-lapply(igraphlist,function(x){
        if(length(unique(igraph::E(x)$multiplex_status))>0){
          linetype_map[unique(igraph::E(x)$multiplex_status)]
        }else{NULL}
      })
      names(edgelinetypemaplist)<-net_name
    }else if(plot_type %in% c("differential_network","differential_subnetwork")){
      edgelinetypemaplist<-lapply(igraphlist,function(x){
        if(length(unique(igraph::E(x)$multiplex_status))>0){
          linetype_map[unique(igraph::E(x)$multiplex_status)]
        }else{NULL}
      })
      names(edgelinetypemaplist)<-net_name

    }else{
      edgelinetypemap<-NULL
    }
    
    ##################################################### 

    plotlist<-list()
    is_related_node<-list()
    for (k in 1:length(net_name)) { 
      network_name=net_name[k]
      igraph<-igraphlist[[network_name]]
      if(plot_type %in% c("differential_network","differential_subnetwork")){
        color_type<-color_typelist[k]
        edge_color_type<-edge_color_typelist[k]
      }
      if(color_type=="Class"){
        colormap <- color_mapping[unique(igraph::V(igraph)$Class)]
      }else if(color_type=="membership"){
        colormap <- color_mapping_cluster[unique(igraph::V(igraph)$membership)]
      }else if(color_type=="FC"){
        colormap <-grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
      }else if(color_type=="Log2FC"){
        colormap <-grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
      }else if(color_type=="Normalized mean"){
        colormap <-grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
        
      }
      if(plot_type %in% c("sub_network","diff_subnetwork","control_subnetwork","case_subnetwork","differential_subnetwork")){
        LAYOUT<-LAYOUT_list[[network_name]]
        if(exists("input_layout")){
          if(!is.null(input_layout)){
            LAYOUT=input_layout
            LAYOUT_list[[network_name]]<-LAYOUT
          }
        }

      }else{
        if(exists("input_layout")){
        if(!is.null(input_layout)){
          LAYOUT=input_layout
        }}
      }
      
      if(plot_type %in% c("diff_network","diff_overall_cluster_network","diff_subnetwork") ||
        ((k %% 3 == 0)  && (plot_type %in% c("differential_subnetwork","differential_network")))){
    
          igraph::E(igraph)$edge_width<-0.5
        
        focus_name1=NULL
        if("all" %in% focus){
          focus_name1<-unique(igraph::E(igraph)$cor_status)
        }else{
          focus_name1<-focus
        }
        

        igraph::E(igraph)$temp_cor_status<-ifelse(igraph::E(igraph)$cor_status %in% focus_name1, igraph::E(igraph)$cor_status, "Other")
        diff_edge_status <- c(Only_in_case,Enhanced_in_case, Only_in_control,Enhanced_in_control, "Conflict relation","Non-significant","Neither group","Other")
        diff_cor_status_colormap <- setNames(c("#d31b16","#FF00FF","#0000FF", "#00FFFF","#f28500","#E0E08C" ,"#cfc0bb","gray"), diff_edge_status)
        
        related_nodes <- unique(unlist(igraph::ends(igraph, which(igraph::E(igraph)$cor_status %in% focus_name1))))
        is_related_node[[k]] <- igraph::V(igraph)$name %in% related_nodes
      }
      if(plot_type %in% c("control_subnetwork","case_subnetwork","diff_subnetwork",
                          "differential_subnetwork","differential_network")){
        edgelinetypemap<-edgelinetypemaplist[[network_name]]
      }
      
      if(plot_type %in% c("diff_network","diff_overall_cluster_network","diff_subnetwork") ||
         ((k %% 3 == 0)  && (plot_type %in% c("differential_subnetwork","differential_network")))){
        edgecolormap<-diff_cor_status_colormap[unique(igraph::E(igraph)$temp_cor_status)]
      }else{
        edgecolormap <- edgecolormapping[unique(igraph::E(igraph)$color)]
      }

      xlim_max=max(LAYOUT[,1])+ image_margin_size
      xlim_min=min(LAYOUT[,1])- image_margin_size
      ylim_max=max(LAYOUT[,2])+ image_margin_size
      ylim_min=min(LAYOUT[,2])- image_margin_size
      plot <- run_ggraph_plot(igraph=igraph,colormap=colormap,edgecolormap=edgecolormap, color_type =color_type,plot_title_size=plot_title_size,
                              plot_layout=LAYOUT,show_node_legend=show_node_legend,show_edge_legend=show_edge_legend,alpha=alpha,
                              edgelinetypemap=edgelinetypemap,edge_color_type=edge_color_type,case_control=case_control) 
      plot <-plot + xlim(xlim_min, xlim_max)+ 
        ylim(ylim_min, ylim_max)
      
      # #节点大小调整
      if(length(unique(igraph::V(igraph)$size))>1){
        if(plot_type=="sub_network"){
          plot <- plot + ggplot2:: scale_size_continuous(range = c(4, 12))
        }else{
          plot <- plot + ggplot2:: scale_size_continuous(range = c(3, 6))
        }
      }else{
        if(is.null(node_size)){
          plot <- plot + ggplot2:: scale_size(range = c(6, 6))
        }else{
          plot <- plot + ggplot2:: scale_size(range = c(as.numeric(node_size), as.numeric(node_size)))
        }
        
      }
      
      #线宽调整
      if(!(plot_type %in% c("diff_network","diff_overall_cluster_network","diff_subnetwork")) ||
         !((k %% 3 == 0)  && (plot_type %in% c("differential_subnetwork","differential_network")))){
       plot <- plot + ggraph::scale_edge_width(range = c(0.2, 1),limits = c(0, 1))
      }
      if(plot_type %in% c("stable_test","case_stable_test","control_stable_test")){
        plotlist[[network_name]]<- plot +  ggplot2::ggtitle(network_name)+# ggplot2::ggtitle(paste0("boot_",network_name))+
          ggplot2::theme(plot.background = ggplot2::element_rect(color = "black"))
      }else if(plot_type %in% c("sub_network","diff_subnetwork","control_subnetwork","case_subnetwork",
                                "differential_subnetwork","differential_network" )){
        plotlist[[network_name]]<- plot + ggplot2::ggtitle(network_name)
      }else{
        plotlist[[network_name]]<- plot
      }
    }
    if(plot_type %in% c("diff_network","diff_overall_cluster_network","diff_subnetwork","differential_subnetwork","differential_network")){
      names(is_related_node)<-net_name
    }
    #######################################################
    if(datasave){
      if(plot_type %in% c("stable_test","case_stable_test","control_stable_test")){
        lapply(net_name,function(net_name){
          nodes<-net_list[[net_name]]$nodes |>
            dplyr::mutate(x=LAYOUT[,1],y=LAYOUT[,2]) 
          edges<-net_list[[net_name]]$edges |>
            dplyr::filter(abs(cor)>R_threshold)
          net_name <- gsub(" ", "", net_name)
          # readr::write_delim(nodes,file.path(outdir,paste0("nodes_",net_name,"_",group_name,".txt")),delim="\t")
          # readr::write_delim(edges,file.path(outdir,paste0("edges_",net_name,"_",group_name,".txt")),delim="\t")
        })
      }else if(plot_type %in% c("overall_network","overall_culster_network",
                                "case_overall_network","control_overall_network","interaction_network",
                                "case_multi_network","control_multi_network","diff_network","diff_overall_cluster_network"
      )){
        nodes<-nodes |>
          dplyr::mutate(x=LAYOUT[,1],y=LAYOUT[,2]) 
        # readr::write_delim(nodes,file.path(outdir,paste0("nodes_",group_name,".txt")),delim="\t")
        # readr::write_delim(edges,file.path(outdir,paste0("edges_",group_name,".txt")),delim="\t")
      }else if(plot_type %in% c("sub_network")){
        lapply(net_name,function(net_name){
          LAYOUT<-LAYOUT_list[[net_name]]
          nodes<-net_list[[net_name]]$nodes |>
            dplyr::mutate(x=LAYOUT[,1],y=LAYOUT[,2]) 
          edges<-net_list[[net_name]]$edges |>
            dplyr::filter(abs(cor)>R_threshold)
          # readr::write_delim(nodes,file.path(outdir,paste0("nodes_",net_name,"_",group_name,".txt")),delim="\t")
          # readr::write_delim(edges,file.path(outdir,paste0("edges_",net_name,"_",group_name,".txt")),delim="\t") 
        })
      }else if(plot_type %in% c("diff_subnetwork","control_subnetwork","case_subnetwork")){
        lapply(net_name,function(net_name){
          LAYOUT<-LAYOUT_list[[net_name]]
          nodes<-net_list[[net_name]]$nodes |>
            dplyr::mutate(x=LAYOUT[,1],y=LAYOUT[,2]) 
          edges<-net_list[[net_name]]$edges
          # readr::write_delim(nodes,file.path(outdir,paste0("nodes_",net_name,"_",group_name,".txt")),delim="\t")
          # readr::write_delim(edges,file.path(outdir,paste0("edges_",net_name,"_",group_name,".txt")),delim="\t") 
        })
        
      }else if(plot_type %in% c("differential_subnetwork","differential_network")){
        lapply(net_name,function(net_name1){
          LAYOUT<-LAYOUT_list[[net_name1]]
          nodes<-net_list[[net_name1]]$nodes |>
            dplyr::mutate(x=LAYOUT[,1],y=LAYOUT[,2]) 
          edges<-net_list[[net_name1]]$edges
          # readr::write_delim(nodes,file.path(outdir,paste0("nodes_",net_name1,".txt")),delim="\t")
          # readr::write_delim(edges,file.path(outdir,paste0("edges_",net_name1,".txt")),delim="\t") 
        })
      }
    }
    
    if(add_enrichement){
      if(plot_type %in% c("overall_culster_network")){
        if(is(Network, "Stable_SubNetwork")) {
          add_anno<-run_add_annotation(Network=Network,plot_layout=LAYOUT)
          if(datasave){
          # readr::write_delim(add_anno,file.path(outdir,paste0("add_anno_",group_name,".txt")),delim="\t")
          }
          plotlist <- lapply(plotlist, function(plot) {

            plot + ggrepel::geom_text_repel(
              data = other,
              ggplot2::aes(x = x, y = y, label = enrichment),
              size = 4,
              fontface = "bold",
              color = "black",
              bg.color = "white",
              bg.r = 0.15,
              seed = 123
            )
          })
        }
      }
    }
    
   if(node_name_type=="name"){
      node_name_vector="feature_name"
    }else{
      node_name_vector="name"
    }
    if (show_node_name){
      # if(plot_type!="enrichment"){
      # if(all(!is.null(is_related_node))){
      if(plot_type %in% c("diff_network","diff_overall_cluster_network","diff_subnetwork")){
        plotlist <- lapply(names(plotlist), function(plot) {
          plotlist[[plot]] + ggraph::geom_node_text(ggplot2::aes(label = ifelse(is_related_node[[plot]], 
                              !!rlang::sym(node_name_vector), "")), vjust = -0.5,size = node_name_size)
        })  
        
      }else{
          plotlist <- lapply(plotlist, function(plot) {
            plot + ggraph::geom_node_text(ggplot2::aes(label = !!rlang::sym(node_name_vector)), vjust = -0.5,size = node_name_size)
          })  

      }
    
    }
    
    if (interactive){
      # if(plot_type!="enrichment"){
      plotlist <- lapply(plotlist, function(plot) {
        plot <- ggiraph::girafe(ggobj = plot)
      })  
    }
    
    ################################
  


  
  

  plotlist <- na.omit(plotlist)
  plotlist <- plotlist[!sapply(plotlist, is.null)]
  if(length(plotlist)>1){
    if(plot_type=="differential_network"){
      plot<-patchwork::wrap_plots(plotlist, ncol = 3, nrow = 1)
    }else if(plot_type=="differential_subnetwork"){
      plot<-patchwork::wrap_plots(plotlist, ncol = 3, nrow = length(plotlist)/3)
    }else{
      # plot<-patchwork::wrap_plots(plotlist, ncol = 4, nrow = 2)#f1
      npar<-ceiling(sqrt(length(plotlist)))
      plot<-patchwork::wrap_plots(plotlist, ncol = npar, nrow = npar)
    }
   
  }else if(length(plotlist)==1){
    if(is.list(plotlist)){
      plot<-plotlist[[1]]
    }else{
      plot<-plotlist
    }
  }else{
    message(paste("Enrichment plot is empty."))
  }
  plot
  
  if(!(plot_type %in% c("enrichment","differential_subnetwork_enrichment","subnetwork_enrichment"))){
    if(exists("net_list")){
      return(new("Network_show",
                 data=net_list,
                 igraph=igraphlist,
                 plot = plot
      )) 
    }else{
      return(new("Network_show",
                 data=list(net=list(nodes=nodes,edges=edges,plot_layout=LAYOUT)),
                 igraph=igraphlist,
                 plot = plot
      ))
    }
  }else if(plot_type %in% c("enrichment","differential_subnetwork_enrichment","subnetwork_enrichment")){
    return(new("Enrichment_show",
               data=annotations_filterlist,
               plot = plot
    ))
  }
  
}




#' Generate color mapping for node classes
#'
#' Internal function that creates a named color vector for node classes using
#' RColorBrewer qualitative palettes. Automatically expands the palette if more
#' colors are needed than available.
#'
#' @param annotation_table data.frame, must contain a 'Class' column with node
#'                         class annotations
#'
#' @return A named character vector mapping class names to colors
#'
#' @details
#' Color selection:
#' 1. Uses RColorBrewer qualitative palettes (Set1, Dark2, Set2, Set3)
#' 2. If unique classes exceed available colors, uses colorRampPalette to
#'    interpolate additional colors
#' 3. Returns colors in order of unique class appearance
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices colorRampPalette
#'
#' @keywords internal
#'
#' @noRd
run_color<-function(annotation_table=NULL){
  unique_classes<-unique(annotation_table$Class)
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  if(length(unique_classes)>length(cluster_Palette)){
    cluster_Palette <-grDevices::colorRampPalette(cluster_Palette)(length(unique_classes))
  }
  color_mapping <- setNames(cluster_Palette[1:length(unique_classes)], unique_classes)
  return(color_mapping)
}


#' Convert node and edge data frames to igraph object
#'
#' Internal function that creates an igraph object from node and edge data frames,
#' calculates centrality measures, and assigns edge colors based on correlation
#' sign.
#'
#' @param nodes data.frame, node information with at least a 'node' column
#' @param edges data.frame, edge information with 'from', 'to', and optionally
#'              'cor' columns
#'
#' @return An igraph object with:
#'   \item{color}{Edge attribute: "positive", "negative", or "other"}
#'   \item{edge_width}{Edge attribute: absolute correlation value or 0.5 if NA}
#'   \item{size}{Node attribute: initial size set to 1}
#'   \item{betweenness}{Node centrality measure}
#'   \item{degree}{Node degree}
#'   \item{eigenvector}{Node eigenvector centrality}
#'
#' @importFrom igraph graph_from_data_frame E V betweenness degree eigen_centrality
#'
#' @keywords internal
#'
#' @noRd
run_igraph <- function(nodes=NULL,edges=NULL){
  net1 <- igraph::graph_from_data_frame(d=as.data.frame(edges),vertices=as.data.frame(nodes),directed = F)

  igraph::E(net1)$color <- ifelse(is.na(igraph::E(net1)$cor), "other",
                                  ifelse(igraph::E(net1)$cor > 0,"positive","negative"))
  igraph::E(net1)$edge_width<-ifelse(is.na(igraph::E(net1)$cor), 0.5, igraph::E(net1)$cor)
  igraph::E(net1)$edge_width<-abs(igraph::E(net1)$edge_width)
  igraph::V(net1)$size<-1
  
  igraph::V(net1)$betweenness <- igraph::betweenness(
    graph = net1,
    v = igraph::V(net1),
    directed = FALSE
  )
  igraph::V(net1)$degree <- igraph::degree(
    graph = net1,
    v= igraph::V(net1)
  )
  igraph::V(net1)$eigenvector <- igraph::eigen_centrality(
    graph = net1
    , directed = FALSE)$vector
  return(net1)
}


#' Create ggplot2 network visualization
#'
#' Internal function that generates a ggplot2/ggraph network plot with customizable
#' node and edge aesthetics.
#'
#' @param igraph igraph object, network to visualize
#' @param colormap function or named vector, color mapping for nodes
#' @param edgecolormap named vector, color mapping for edges
#' @param color_type character, node coloring scheme
#' @param plot_layout matrix, layout coordinates or layout algorithm name
#' @param edge_color_type character, edge coloring scheme ("color" or custom)
#' @param edgelinetypemap named vector, line type mapping for edges
#' @param plot_title_size numeric, font size for plot title
#' @param alpha numeric, edge transparency
#' @param show_edge_legend logical, whether to show edge legend
#' @param show_node_legend logical, whether to show node legend
#' @param case_control character, label for case-control comparison
#'
#' @return A ggplot2 object with the network visualization
#'
#' @importFrom ggraph ggraph geom_edge_link scale_edge_color_manual
#' @importFrom ggplot2 aes scale_fill_manual scale_fill_gradientn theme_void
#' @importFrom ggiraph geom_point_interactive
#'
#' @keywords internal
#'
#' @noRd
run_ggraph_plot <- function(igraph=NULL, colormap=NULL,edgecolormap=NULL, color_type = "Class",plot_layout="fr",
                            edge_color_type="color",edgelinetypemap=NULL,plot_title_size=1,
                            alpha = 1,show_edge_legend=FALSE,show_node_legend=FALSE,case_control="(case:control)") {
  plot <- ggraph::ggraph(igraph, layout = plot_layout)
  if(all(!is.null(edgelinetypemap))){
    plot <- plot + ggraph::geom_edge_link(
      ggplot2::aes(edge_color =!!rlang::sym(edge_color_type), edge_width = edge_width, linetype = multiplex_status),
      alpha = alpha,
      show.legend = c(edge_color = show_edge_legend, edge_width = FALSE, linetype = TRUE)
    ) +
      ggraph::scale_edge_linetype_manual(values = edgelinetypemap, name = "Type",
                                         guide = ggplot2::guide_legend(order=1)) 
    ##############################
  }else{
    plot <- plot + ggraph::geom_edge_link(
      ggplot2::aes(
        edge_color =!!rlang::sym(edge_color_type), edge_width = edge_width),
      alpha = alpha,
      show.legend = c(edge_color = show_edge_legend, edge_width = FALSE, linetype = FALSE)
      #show.legend = FALSE
    )
  }
  if(edge_color_type=="color"){
    plot <- plot + ggraph::scale_edge_color_manual(values = edgecolormap, name = "Correlation",
                                                   guide = ggplot2::guide_legend(order=2))
  }else{
    plot <- plot + ggraph::scale_edge_color_manual(values = edgecolormap, name = "Correlation status",
                                                   guide = ggplot2::guide_legend(order=2))
  }
  ######################################
  
  
  omics_num<-length(unique(igraph::V(igraph)$omics_name))
  if(omics_num>1){
    nodeshape<-c(21,24,22,23,25)
    omics_levels <- sort(unique(igraph::V(igraph)$omics_name))
    omics_num <- min(length(omics_levels), length(nodeshape))
    shape_values <- setNames(nodeshape[1:omics_num], omics_levels)
    plot <- plot + ggiraph::geom_point_interactive(
      ggplot2::aes(x=x,y=y,size = size,tooltip = name,data_id = name,
                   fill = !!rlang::sym(color_type),shape=omics_name),
      # shape = 21,              
      colour = "black",
      show.legend = c(size = FALSE, fill = show_node_legend,shape=show_node_legend)
    ) +
      ggplot2::scale_shape_manual(name="Omics",values =shape_values) 
  }else{
    plot <- plot + ggiraph::geom_point_interactive(
      ggplot2::aes(x=x,y=y,tooltip = name,data_id = name,
                   size = size, fill = !!rlang::sym(color_type)),
      shape = 21,             
      colour = "black",
      show.legend = c(size = FALSE, fill = show_node_legend)
    ) 
  }
  
  if(is.function(colormap)){
    if(color_type=="Normalized mean"){
      plot <- plot + ggplot2::scale_fill_gradientn(colours = colormap(100),limits = c(0, 1),#scale_color_gradientn
                                                   name = color_type,
                                                   guide = ggplot2::guide_colorbar(order=3)) #
    }else{
      plot <- plot + ggplot2::scale_fill_gradientn(colours = colormap(100), name = paste0(color_type,case_control),
                                                   guide = ggplot2::guide_colorbar(order=3)) #scale_color_gradientn
    }
  }else{
    plot <- plot + ggplot2::scale_fill_manual(values = colormap, name = color_type,
                                              guide = ggplot2::guide_legend(order=3)) #scale_color_manual
  }
  #主题
  plot <- plot +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = plot_title_size)) +
    coord_equal()

  return(plot)
}


#' Add centrality measures as node sizes
#'
#' Internal function that sets node sizes based on specified centrality measures
#' for use in network visualizations.
#'
#' @param igraph igraph object, network to modify
#' @param add_Centrality character vector, centrality measures to use for sizing.
#'                       Options: "betweenness", "degree", "eigenvector".
#'                       Only the first element is used.
#'
#' @return Modified igraph object with node 'size' attribute set to centrality values
#'
#' @keywords internal
#'
#' @noRd
run_add_Centrality<-function(igraph=NULL,add_Centrality="betweenness"){
  if(add_Centrality[1]=="betweenness"){
    node_center<-igraph::V(igraph)$betweenness
  }else if(add_Centrality[1]=="degree"){
    node_center<-igraph::V(igraph)$degree
  }else if(add_Centrality[1]=="eigenvector"){
    node_center<-igraph::V(igraph)$eigenvector
  }
  node_size <-node_center
  igraph::V(igraph)$size <- node_size
  return(igraph)
}


#' Generate centrality scatterplot
#'
#' Internal function that creates a scatterplot showing centrality measures for
#' network nodes, useful for identifying influential nodes.
#'
#' @param igraph igraph object, network with calculated centrality measures
#' @param add_Centrality character vector, centrality measures to plot
#' @param nodes data.frame, node annotation table
#' @param node_name_type character, "id" or "name" for node labeling
#'
#' @return A list containing:
#'   \item{Centrality_plot}{ggplot2 scatterplot}
#'   \item{Centrality_data}{Data frame with centrality values in long format}
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter mutate arrange
#' @importFrom ggplot2 ggplot aes geom_path geom_point theme_bw
#'
#' @keywords internal
#'
#' @noRd
run_plot_Centrality<-function(igraph=NULL,add_Centrality="betweenness",
                              nodes=nodes,node_name_type="id"){
  if(node_name_type=="name"){
    Centrality_data_long<-data.frame(id=nodes$feature_name,eigenvector=igraph::V(igraph)$eigenvector,
                                     degree=igraph::V(igraph)$degree,
                                     betweenness=igraph::V(igraph)$betweenness)
  }else{
    Centrality_data_long<-data.frame(id=names(igraph::V(igraph)),eigenvector=igraph::V(igraph)$eigenvector,
                                     degree=igraph::V(igraph)$degree,
                                     betweenness=igraph::V(igraph)$betweenness)
  }
  Centrality_data_long<-Centrality_data_long |>
    tidyr::pivot_longer(
      cols = c("betweenness", "degree", "eigenvector"), 
      names_to = "Centrality",
      values_to = "value"
    ) |>
    dplyr::filter(Centrality %in% add_Centrality)
  id_Centrality_level= names(sort(tapply(Centrality_data_long$value,Centrality_data_long$id, mean),decreasing = FALSE))
  Centrality_data<-Centrality_data_long |>
    dplyr::mutate(id = factor(id, levels = unique(id_Centrality_level))) |>
    dplyr::mutate(Centrality = factor(Centrality, levels = unique(Centrality)))   |>
    dplyr::arrange(Centrality, id)
  Centrality_plot<-ggplot2::ggplot(Centrality_data, ggplot2::aes(x = value, y =id, color = Centrality, group = Centrality)) +
    ggplot2::geom_path() + ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::geom_point()+
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 15, face="bold"),
      plot.title = ggplot2::element_text(hjust = 0.5))+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  return(list(Centrality_plot=Centrality_plot,Centrality_data=Centrality_data))
}


#' Create comprehensive subnetwork visualization for differential networks
#'
#' This function generates three network visualizations for a subnetwork:
#' case group network, control group network, and differential network,
#' saved as separate files and a combined plot.
#'
#' @param subnet_dir character, directory containing node and edge files
#' @param subnet_name character, name of the subnetwork
#' @param group_name character, comparison label (e.g., "T-vs-N")
#' @param omics1_name character, name of first omics type (e.g., 'Pro')
#' @param omics2_name character, name of second omics type (e.g., 'Phos')
#' @param edge_color_pos character, color for positive correlations
#' @param edge_color_neg character, color for negative correlations
#' @param Enhanced_in_N character, color for edges enhanced in control group
#' @param Enhanced_in_T character, color for edges enhanced in treatment group
#' @param Only_in_N character, color for edges only in control group
#' @param Only_in_T character, color for edges only in treatment group
#' @param Conflict_relation character, color for conflicting edges
#' @param fill_gradientn_color character vector, colors for continuous gradients
#' @param plot_title_size numeric, font size for plot titles
#' @param legend_title_size numeric, font size for legend titles
#' @param legend_text_size numeric, font size for legend text
#' @param node_label_size numeric, font size for node labels
#' @param point_size numeric, size of node points
#'
#' @return No return value. Saves three PNG files to the subnet_dir:
#'   \item{subnet_name_group1_plot.png}{Case group network visualization}
#'   \item{subnet_name_group2_plot.png}{Control group network visualization}
#'   \item{subnet_name_diff_plot.png}{Differential network visualization}
#'   \item{subnet_name_group_name_plot.png}{Combined three-panel plot}
#'
#' @importFrom igraph graph_from_data_frame V E delete_edges
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggplot2 aes scale_edge_color_manual scale_fill_gradientn
#' @importFrom patchwork plot_layout
#'
#' @export
Differential_subnetwork_plot <- function(subnet_dir,subnet_name,group_name = 'T-vs-N', omics1_name = 'Pro',
                                         omics2_name = 'Phos',
                                         edge_color_pos = "#9b6a65",
                                         edge_color_neg = "#5d8992", 
                                         Enhanced_in_N = "#5d8992", 
                                         Enhanced_in_T = "#9b6a65",
                                         Only_in_N = "#0c2b32",
                                         Only_in_T = "#381512",
                                         Conflict_relation = '#808080',
                                         fill_gradientn_color = c("#175663", "#dce6e9", "#90362d"),
                                         plot_title_size = 20,        
                                         legend_title_size = 12,      
                                         legend_text_size = 12,       
                                         node_label_size = 6,         
                                         point_size = 10){
  nodes_path <- file.path(subnet_dir,paste0("nodes_",subnet_name,"_",group_name,"_",subnet_name,".tsv"))
  edges_path <- file.path(subnet_dir,paste0("edges_",subnet_name,"_",group_name,"_",subnet_name,".tsv"))
  
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
      aes(fill = case_norm_mean, shape = Class),
      size = point_size,      
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
      plot.title = element_text(hjust = 0.5, size = plot_title_size, face = "bold"),  
      legend.title = element_text(size = legend_title_size, face = "bold"),           
      legend.text = element_text(size = legend_text_size),                            
    ) +
    ggtitle(paste0("Group ",group1," net")) +
    guides(
      fill = guide_colorbar(barheight = unit(4, "cm")),
      shape = guide_legend(
        override.aes = list(
          fill = "grey80",
          size = 4,
          stroke = 0.5
        ))
    )
  T_subnet_lable <- T_subnet +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      # size = 6,
      size = node_label_size,    
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
      size = point_size,         
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
      values = current_mapping, 
      guide = guide_legend(override.aes = list(fill = "grey50", size = 4))
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = plot_title_size, face = "bold"),  
      legend.title = element_text(size = legend_title_size, face = "bold"),          
      legend.text = element_text(size = legend_text_size),                           
    ) +
    ggtitle(paste0("Group ",group2," net")) +
    guides(
      fill = guide_colorbar(barheight = unit(4, "cm")),
      shape = guide_legend(
        override.aes = list(
          fill = "grey80",
          size = 4,
          stroke = 0.5
        ))
    )
  N_subnet_lable <- N_subnet+
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      # size = 6,
      size = node_label_size,   
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
      aes(fill = case_norm_mean, shape = Class),
      size = point_size,         
      colour = "transparent",
      stroke = 0
    ) +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = node_label_size,   
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
      plot.title = element_text(hjust = 0.5, size = plot_title_size, face = "bold"),  
      legend.title = element_text(size = legend_title_size, face = "bold"),           
      legend.text = element_text(size = legend_text_size),                            
    ) +
    ggtitle(paste0(group_name," net")) +
    guides(
      fill = guide_colorbar(barheight = unit(4, "cm")),
      shape = guide_legend(
        override.aes = list(
          fill = "grey80",
          size = 4,
          stroke = 0.5
        ))
    )
  ggsave(file.path(subnet_dir,paste0(subnet_name,"_diff_plot.png")),plot = diff_subnet, width = 10, height = 8, dpi = 300)
  
  combined_plot <- T_subnet + N_subnet + diff_subnet+ plot_layout(ncol = 3)
  ggsave(file.path(subnet_dir,paste0(subnet_name,"_",group_name,"_plot.png")), plot = combined_plot, width = 30, height = 8, dpi = 300)
  
}


#' Create community detection plots for differential networks
#'
#' This function generates network visualizations showing module/community
#' structures in differential networks. Creates both an overall module plot
#' and a faceted plot with separate panels for each module.
#'
#' @param node_file_path character, path to node table file (TSV format)
#' @param edge_file_path character, path to edge table file (TSV format)
#' @param outdir character, output directory for saving plots
#' @param omics1_name character, name of first omics type (e.g., 'Pro')
#' @param omics2_name character, name of second omics type (e.g., 'Phos')
#' @param ModuleSize_show integer, minimum module size to display. Modules
#'                        smaller than this are grouped as "other". Default 0.
#' @param top_module_num integer, number of top modules to display. Default 20.
#' @param title_size numeric, font size for plot title. Default 20.
#' @param legend_title_size numeric, font size for legend titles. Default 15.
#' @param legend_text_size numeric, font size for legend text. Default 13.
#' @param strip_text_size numeric, font size for facet strip text. Default 20.
#' @param node_label_size numeric, font size for node labels. Default 3.
#' @param point_size_range numeric vector, range for node point sizes (min, max).
#'                         Default c(4, 8).
#' @param edge_alpha numeric, transparency for edges. Default 0.3.
#'
#' @return No return value. Saves two PNG files to outdir:
#'   \item{network_by_top_modules.png}{Overall network with module coloring}
#'   \item{network_by_separate_top_modules.png}{Faceted plot with separate panels for each module}
#'
#' @details
#' Processing steps:
#' 1. Load node and edge data from files
#' 2. Filter out nodes with membership "subnet_other"
#' 3. Remove edges with non-significant correlation status
#' 4. Build igraph object and convert to tidygraph
#' 5. Filter modules below ModuleSize_show threshold
#' 6. Select top_module_num modules by size
#' 7. Generate color mapping using RColorBrewer and viridis palettes
#' 8. Create overall module plot with convex hulls around modules
#' 9. Create faceted plot with separate panels per module
#'
#' @importFrom igraph graph_from_data_frame induced_subgraph V
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggplot2 aes scale_color_manual scale_fill_manual theme_void
#' @importFrom ggforce geom_mark_hull
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis
#'
#' @export
diff_net_community_detection_plot <- function(node_file_path,edge_file_path,outdir, omics1_name = 'Pro',omics2_name = 'Phos',ModuleSize_show = 0,top_module_num = 20,
                                              title_size = 20,
                                              legend_title_size = 15,
                                              legend_text_size = 13,
                                              strip_text_size = 20,
                                              node_label_size = 3,
                                              point_size_range = c(4, 8),
                                              edge_alpha = 0.3) {
  
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
    geom_edge_link(color = "grey50", alpha = edge_alpha, width = 0.6)+ 
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
    theme(
      strip.text = element_text(face = "bold", size = strip_text_size), 
      plot.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = legend_title_size, face = "bold"),
      legend.text = element_text(size = legend_text_size),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = title_size,
        margin = margin(t = 15, b = 15)
      )
    )+
    labs(title = paste0("Differential Protein and Phosphorylation Sites Network with top ",top_module_num," Modules"))
  # top_module_p
  ggsave(paste0(outdir,"/network_by_top_modules.png"), plot = top_module_p, width = 16, height = 10, dpi = 300)
  
  top_module_p3 <- ggraph(top_module_g_tbl, layout = layout_matrix) +
    geom_edge_link(color = "grey50", alpha = edge_alpha, width = 0.6)+ 
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
    theme(
      strip.text = element_text(
        face = "bold",
        size = strip_text_size,  # 原来是40
        margin = margin(b = 5)
      ),
      plot.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = legend_title_size, face = "bold"),  
      legend.text = element_text(size = legend_text_size),  # 原来是25
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = title_size,
        margin = margin(t = 15, b = 15)
      )
    )+
    labs(title = paste0("Differential Protein and Phosphorylation Sites Network by top ",top_module_num," Modules"))
  num_modules <- length(unique(vertex_attr(top_module_g_tbl, "membership")))
  expected_cols <- ceiling(sqrt(num_modules))
  ggsave(paste0(outdir,"/network_by_separate_top_modules.png"), plot = top_module_p3, width = 24, height = 6 *expected_cols, dpi = 300)
}




