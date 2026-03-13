#depends:ggplot2 dplyr igraph qgraph stringr readr RColorBrewer rlang ggiraph patchwork ggraph tidyr tibble  


#' Network visualization and analysis class
#'
#' This class stores network data, igraph objects, and visualization plots
#' for various types of network analyses.
#'
#' @slot data List containing network data components
#' @slot igraph ANY type storing igraph object(s)
#' @slot plot ANY type storing visualization plot(s)
#' @slot picture_width Numeric width for plot output
#' @slot picture_height Numeric height for plot output
#' @slot other ANY type for additional results (e.g., topological analysis)
#'
#' @name Network_show-class
#' @rdname Network_show-class
#' @exportClass Network_show
#'
#' @examples
#' # Create a simple Network_show object for demonstration
#' demo_net <- new("Network_show",
#'                 data = list(nodes = data.frame(node = c("A", "B"), 
#'                                               Class = c("Type1", "Type2")),
#'                            edges = data.frame(from = "A", to = "B", cor = 0.8,cor=0.5,color="red")),
#'                 igraph = igraph::graph_from_data_frame(
#'                   data.frame(from = "A", to = "B",cor=0.5,color="red"),vertices=data.frame(name=c("A","B"),Class = c("Type1", "Type2")), 
#'                   directed = FALSE),
#'                 plot = ggraph::ggraph(igraph::graph_from_data_frame(
#'                 data.frame(from = "A", to = "B",cor=0.5,color="red"),vertices=data.frame(name=c("A","B"),Class = c("Type1", "Type2")), directed = FALSE)) + 
#'                 ggraph::geom_edge_link(ggplot2::aes(edge_color =color , edge_width = cor)) + 
#'                 ggiraph::geom_point_interactive(ggplot2::aes(x=x,y=y,color =Class ,tooltip = name,data_id = name,size = 1)) +
#'                 ggplot2::theme_void(),
#'                 picture_width = 10,
#'                 picture_height = 8,
#'                 other = NULL)
#' demo_net
setClass("Network_show", slots = c(
  data="list",
  igraph = "ANY",
  plot="ANY",
  picture_width="numeric",
  picture_height="numeric",
  other= "ANY"
))

#' Enrichment visualization class
#'
#' This class stores enrichment data and visualization plots.
#'
#' @slot data List containing enrichment data
#' @slot plot ANY type storing visualization plot(s)
#' @slot picture_width Numeric width for plot output
#' @slot picture_height Numeric height for plot output
#'
#' @name Enrichment_show-class
#' @rdname Enrichment_show-class
#' @exportClass Enrichment_show
#'
#' @examples
#' # Create a simple Enrichment_show object for demonstration
#' demo_enrich <- new("Enrichment_show",
#'                   data = list(data.frame(Pathway = c("Pathway1", "Pathway2"),
#'                              Pvalue = c(0.001, 0.01),
#'                              RichFactor = c(0.5, 0.3),Count=c(5,3))),
#'                   plot = ggplot2::ggplot(data.frame(Pathway = c("Pathway1", "Pathway2"),Pvalue = c(0.001, 0.01),RichFactor = c(0.5, 0.3),Count=c(5,3))) + 
#'                   ggplot2::geom_point(ggplot2::aes(x=RichFactor,y=Pathway,size=Count,colour=Pvalue)),
#'                   picture_width = 10,
#'                   picture_height = 8)
#' demo_enrich
setClass("Enrichment_show", slots = c(
  data="list",
  plot="ANY"  ,
  picture_width="numeric",
  picture_height="numeric"
))

#' Topological analysis result class
#'
#' This class stores topological analysis results and plots.
#'
#' @slot data List containing topological data
#' @slot plot ANY type storing visualization plot(s)
#' @slot picture_width Numeric width for plot output
#' @slot picture_height Numeric height for plot output
#'
#' @name Topological_result-class
#' @rdname Topological_result-class
#' @exportClass Topological_result
#'
#' @examples
#' # Create a simple Topological_result object for demonstration
#' demo_topological <- new("Topological_result",
#'                        data = list(data.frame(id = c("a","b","c" ,"d"),Centrality=c("betweenness","betweenness","degree","degree"), value = c(1.2,1.4,5, 8))),
#'                        plot = ggplot2::ggplot(data = data.frame(id = c("a","b","c" ,"d"),Centrality=c("betweenness","betweenness","degree","degree"), value = c(1.2,1.4,5, 8)), 
#'                        ggplot2::aes(x = value, y =id, color = Centrality, group = Centrality)) +
#'                        ggplot2::geom_path() + ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::geom_point(),
#'                        picture_width = 10,
#'                        picture_height = 8)
#' demo_topological
setClass("Topological_result", slots = c(
  data="list",
  plot="ANY"  ,
  picture_width="numeric",
  picture_height="numeric"
))

#' 网络与富集对象的内部显示方法
#'
#' 用于 Network_show 和 Enrichment_show S4 类的内部显示方法。
#' 这些方法仅供包内部使用，不对外导出，也不适合用户直接交互。
#'
#' @param object Network_show 或 Enrichment_show 类的对象
#'
#' @details
#' 这些内部显示方法为调试和开发目的提供格式化的控制台输出。
#' 它们不是公共 API 的一部分，可能会在不通知的情况下发生变化。
#'
#' @name internal-show-methods
#' @aliases show,Network_show-method show,Enrichment_show-method
#' @keywords internal
NULL

#' @rdname internal-show-methods
#' @description
#' Network_show 对象的显示方法。展示关于网络数据结构和可视化的基本信息。
#' 
#' @slot data 底层网络数据（通常是包含节点和边信息的数据框或列表）
#' @slot igraph 表示网络结构的 igraph 对象
#' @slot plot 包含网络可视化内容的 ggplot2 对象
setMethod("show", "Network_show", function(object) {
  cat("=== Network_show ===\n")
  cat("- 数据: object@data\n")
  cat("- 网络 (igraph): object@igraph\n")
  cat("- 图形: object@plot\n")
  print(object@plot)
})

#' @rdname internal-show-methods
#' @description
#' Enrichment_show 对象的显示方法。展示关于富集分析结果和相关可视化的信息。
#' 
#' @slot data 富集分析结果（通常是包含术语、p值、富集倍数等的数据框）
#' @slot plot 可视化富集结果的 ggplot2 对象（如条形图、显著术语的点图等）
setMethod("show", "Enrichment_show", function(object) {
  cat("=== Enrichment_show ===\n")
  cat("- 数据: object@data\n")
  cat("- 图形: object@plot\n")
  print(object@plot)
})

#' @rdname internal-show-methods
#' @description
#' Topological_result 对象的显示方法。展示关于网络拓扑分析结果和可视化的信息。
#' 
#' @slot data 网络拓扑指标（如度分布、介数中心性、聚类系数等）
#' @slot plot 可视化网络拓扑属性的 ggplot2 对象
setMethod("show", "Topological_result", function(object) {
  cat("=== Topological_result ===\n")
  cat("- 数据: object@data\n")
  cat("- 图形: object@plot\n")
  print(object@plot)
})
#' Network visualization and analysis function
#'
#' This function creates various types of network visualizations and analyses
#' from different network objects.
#'
#' @param Network A network object (e.g., Stable_SubNetwork, StableNetwork, etc.)
#' @param plot_type Type of plot to generate. Options include: "stable_test", 
#'        "overall_network", "sub_network", "enrichment", etc.
#' @param stable_num Number of bootstrap networks to display (for stable_test)
#' @param subnetwork_name Names of subnetworks to display (default: "all")
#' @param layout_type Layout algorithm for network visualization. 
#'        Options: "kk", "nicely", "fr"
#' @param input_layout Custom layout matrix (optional)
#' @param richfactor_threshold Threshold for enrichment rich factor (default: 0)
#' @param node_colortype Variable to use for node coloring. Options: 
#'        "Class", "FC", "Log2FC", "Normalized mean", "membership"
#' @param R_threshold Correlation threshold for edge filtering (default: 0)
#' @param focus Which correlation status to focus on (default: "all")
#' @param node_name_size Size of node labels (default: 2)
#' @param node_name_type Type of node labels. Options: "id", "name" (default: "id")
#' @param node_size Custom node size (optional)
#' @param image_margin_size Margin size around plot (default: 0.3)
#' @param show_node_name Whether to show node names (default: FALSE)
#' @param show_node_legend Whether to show node legend (default: FALSE)
#' @param show_edge_legend Whether to show edge legend (default: FALSE)
#' @param plot_title_size Size of plot title (default: 12)
#' @param axis_title_size Size of axis titles (default: 8)
#' @param text_size Size of text elements (default: 8)
#' @param legend_title_size Size of legend title (default: 8)
#' @param legend_text_size Size of legend text (default: 8)
#' @param font_family Font family for text (default: "Arial")
#' @param add_enrichement Whether to add enrichment annotations (default: FALSE)
#' @param add_Centrality Centrality measure to add. Options: "betweenness", 
#'        "degree", "eigenvector" (default: NULL)
#' @param interactive Whether to create interactive plot (default: FALSE)
#' @param centrality_scatterplot Whether to create centrality scatterplot (default: TRUE)
#' @param picture_width Width of output plot (default: 12)
#' @param picture_height Height of output plot (default: 8)
#' @param alpha Transparency level (default: 1)
#'
#' @return Returns a Network_show or Enrichment_show object containing 
#'         visualization and data
#' @export
#'
#' @examples
#' \donttest{
#' # Load example data
#' data("stable_subnetwork_result")
#' 
#' # Create stability test plot
#' network_result <- network_show(
#'   Network = stable_subnetwork_result,
#'   plot_type = "stable_test",
#'   stable_num = 4,
#'   R_threshold = 0.6
#' )
#' 
#' # Display the result
#' network_result@plot
#' 
#' # Create enrichment plot
#' enrichment_result <- network_show(
#'   Network = stable_subnetwork_result,
#'   plot_type = "enrichment",
#'   richfactor_threshold = 0.1
#' )
#' 
#' # Display enrichment result
#' enrichment_result@plot
#' }
network_show<-function(Network=NULL,plot_type="stable_test",stable_num=9,
                       subnetwork_name=c("all"),layout_type="fr",input_layout=NULL,richfactor_threshold=0,#diff_subnetwork_name=c("all")
                       node_colortype="Class",R_threshold=0,focus=c("all"),node_name_size=2,node_name_type="id",node_size=NULL,
                       image_margin_size=0.3,
                       show_node_name=FALSE,show_node_legend=FALSE,show_edge_legend=FALSE,
                       plot_title_size=12,#图标题的大小
                       axis_title_size=8,#轴标题大小
                       text_size=8,#轴刻度文本大小
                       legend_title_size=8,#图例标题大小
                       legend_text_size=8,#图例文字大小
                       font_family="Arial",
                       add_enrichement=FALSE,add_Centrality=NULL,interactive=FALSE,
                       centrality_scatterplot=TRUE,
                       picture_width = 12, picture_height = 8, 
                       alpha=1){ 
  showtext::showtext_auto()#加载字体
  set.seed(123)
  nodes<-NULL
  LAYOUT<-NULL
  edges<-NULL
  plot<-NULL
  plotlist<-NULL
  other<-NULL
  is_related_node<-NULL
  group_name=Network@group_name
  split_result <- strsplit(group_name, "-vs-")[[1]]
  Networkclass=class(Network)[1]
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

  
  if(!(plot_type %in% c("enrichment","differential_subnetwork_enrichment","subnetwork_enrichment"))){ 
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

      net_list <- lapply(net_list, function(bootnet) {
        bootnet$edges <- dplyr::filter(bootnet$edges, abs(cor) > R_threshold)
        return(bootnet)
      })

      
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
            #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else if(is(Network, "Conditional_network")){
          nodes<-Network@network_case@bootnet_result_filter@bootnet_node |>
            dplyr::mutate(`Normalized mean` = case_norm_mean)
          edges<-Network@network_case@bootnet_result_filter@bootnet_edge |>
            #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
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
            #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else if(is(Network, "Conditional_network")){
          nodes<-Network@network_control@bootnet_result_filter@bootnet_node |>
            dplyr::mutate(`Normalized mean` = control_norm_mean)
          edges<-Network@network_control@bootnet_result_filter@bootnet_edge |>
            #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else{
          stop(paste0("‘control_overall_network‘ does not applies to the target object :",class(Network)))
        } 
      }
      ################################
    }else if(plot_type =="overall_cluster_network"){
      if(is(Network, "Stable_SubNetwork")) {
        nodes<-Network@SubNetwork@overall_cluster_network$nodes
        edges<-Network@SubNetwork@overall_cluster_network$edges |>
          dplyr::filter(abs(cor)>R_threshold)
      }else if(is(Network, "SubNetwork")){
        nodes<-Network@overall_cluster_network$nodes
        edges<-Network@overall_cluster_network$edges |>
          dplyr::filter(abs(cor)>R_threshold)
      }else{
        stop(paste0("‘overall_cluster_network‘ does not applies to the target object :",class(Network)))
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
      net_list <- lapply(net_list, function(subnet) {
        subnet$edges <- dplyr::filter(subnet$edges, abs(cor) > R_threshold)
        return(subnet)
      })
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
            #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else if(is(Network, "Conditional_multiplexnetwork")){
          nodes<-Network@network_case@nodes |>
            dplyr::mutate(`Normalized mean` = case_norm_mean)
          edges<-Network@network_case@edges |>
            #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
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
            #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
            dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        }else if(is(Network, "Conditional_multiplexnetwork")){
          nodes<-Network@network_control@nodes |>
            dplyr::mutate(`Normalized mean` = control_norm_mean)
          edges<-Network@network_control@edges |>
            #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
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
        nodes<-Network@Differential_subnetwork@overall_cluster_network$nodes
        edges<-Network@Differential_subnetwork@overall_cluster_network$edges
      }else if(is(Network, "SubNetwork")){
        nodes<-Network@overall_cluster_network$nodes
        edges<-Network@overall_cluster_network$edges
      }else{
        stop(paste0("‘diff_overall_cluster_network‘ does not applies to the target object :",class(Network)))
      }
    }else if(plot_type %in% c("diff_subnetwork")){
      if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork")){
        net_list<-Network@Differential_subnetwork@subnetworks
        #   net_name<-names(Network@Differential_subnetwork@subnetworks)
        # nodes<-Network@Differential_subnetwork@subnetworks[[net_name]]$nodes
        # edges<-Network@Differential_subnetwork@subnetworks[[net_name]]$edges
      }else if(is(Network, "SubNetwork")){
        net_list<-Network@subnetworks
        # net_name<-names(Network@subnetworks)
        # nodes<-Network@subnetworks[[net_name]]$nodes
        # edges<-Network@subnetworks[[net_name]]$edges
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
          #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
          dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        control_nodes<-Network@Conditional_network@network_control@bootnet_result_filter@bootnet_node |>
          dplyr::mutate(`Normalized mean` = control_norm_mean)
        control_edges<-Network@Conditional_network@network_control@bootnet_result_filter@bootnet_edge |>
          #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
          dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        nodes<-Network@Differential_network@diff_nodes
        edges<-Network@Differential_network@diff_edges 
      }else if(is(Network, "Stable_MultiplexNetwork")){ 
        case_nodes<-Network@Conditional_multiplexnetwork@network_case@nodes |>
          dplyr::mutate(`Normalized mean` = case_norm_mean)
        case_edges<-Network@Conditional_multiplexnetwork@network_case@edges |>
          #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
          dplyr::mutate(cor = ifelse(!is.na(cor) & abs(cor) <= R_threshold, NA, cor))
        control_nodes<-Network@Conditional_multiplexnetwork@network_control@nodes |>
          dplyr::mutate(`Normalized mean` = control_norm_mean)
        control_edges<-Network@Conditional_multiplexnetwork@network_control@edges |>
          #如果 cor 不为 NA 且绝对值小于等于 R_threshold，则将其替换为 NA，否则保留原值。
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
    
    if(plot_type %in% c("interaction_network","overall_network","overall_cluster_network",
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
        edges<-net_list[[x]]$edges# |>
      #    dplyr::filter(abs(cor)>R_threshold)
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
      Centrality_data_list <- list()
      Centrality_plot_list <- list()
      
      igraphlist <- lapply(net_name, function(x) {
        nodes <- net_list[[x]]$nodes
        edges <- net_list[[x]]$edges# |>
#          dplyr::filter(abs(cor) > R_threshold)
        igraph <- run_igraph(nodes = nodes, edges = edges)
        
        if(!is.null(add_Centrality)){
          igraph <- run_add_Centrality(igraph = igraph, add_Centrality = add_Centrality)
          
          if(centrality_scatterplot){
            Centralitylist <- run_plot_Centrality(
              igraph = igraph, 
              add_Centrality = add_Centrality,
              nodes = nodes, 
              node_name_type = node_name_type,
              text_size=text_size,
              legend_title_size=legend_title_size,
              legend_text_size=legend_text_size,
              font_family=font_family
            )
            
            Centrality_data <- Centralitylist[["Centrality_data"]]
            Centrality_plot <- Centralitylist[["Centrality_plot"]]
            
            # 将数据存储到列表中，使用网络名称作为键
            Centrality_data_list[[x]] <<- Centrality_data
            Centrality_plot_list[[x]] <<- Centrality_plot
            # 注释掉的保存代码保持不变
            # Centrality_plotname<-paste0("scatterplot_", stringr::str_c(add_Centrality, collapse = "_"),"_",
            #                             x,"_",group_name)
            # if(datasave){
            #   outdir1<-file.path(outdir,Networkclass,"SubNetWork",x,"Topological_analysis")
            #   if (!dir.exists(outdir1)) {
            #     dir.create(outdir1, recursive = TRUE)
            #   }
            #   readr::write_delim(Centrality_data,file.path(outdir1,paste0(Centrality_plotname,".txt")),delim="\t")
            # }
            # outdir1<-file.path(outdir,Networkclass,x,"Topological_analysis")
            # if (!dir.exists(outdir1)) {
            #   dir.create(outdir1, recursive = TRUE)
            # }
            # pdf(file = file.path(outdir1,paste0(Centrality_plotname,".pdf")), 
            #     width = picture_width, height = picture_height)
            # print(Centrality_plot)
            # dev.off()
          }
        }
        igraph
      })
      names(igraphlist)=net_name
      if(!is.null(add_Centrality)){
      for(i in names(igraphlist)) {
        # 获取当前bootnet的nodes数据
        current_nodes <- net_list[[i]]$nodes
        # 获取对应的Centrality数据
        centralityigraph <- igraphlist[[i]]
        if(add_Centrality[1]=="betweenness"){
          node_center<-data.frame(node=igraph::V(centralityigraph)$name,betweenness=igraph::V(centralityigraph)$betweenness)
         
        }else if(add_Centrality[1]=="degree"){
          node_center<-data.frame(node=igraph::V(centralityigraph)$name,degree=igraph::V(centralityigraph)$degree)
        }else if(add_Centrality[1]=="eigenvector"){
          node_center<-data.frame(node=igraph::V(centralityigraph)$name,eigenvector=igraph::V(centralityigraph)$eigenvector)
        }
        merged_nodes<-current_nodes |>
          dplyr::left_join(node_center,by=c("node"="node"))
        # 将合并后的数据保存回net_list，保持数据结构不变
        net_list[[i]]$nodes <- merged_nodes
      }
      }
      other<-new("Topological_result",
          data=Centrality_data_list,
          plot = Centrality_plot_list,
          picture_width=picture_width,
          picture_height=picture_height
      )

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
    }else if(plot_type %in% c("overall_cluster_network","diff_overall_cluster_network")){
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
      if(all(is.null(nodes))){
      list_of_nodes <- lapply(net_list, function(x) x$nodes)
        nodes <- do.call("rbind", list_of_nodes) |>
        dplyr::distinct(node, .keep_all = TRUE)   # 根据node列去重，保留所有其他列
      }
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
        linetypes<-unique(igraph::E(igraph)$multiplex_status)
        linetype_map <-setNames(ifelse(grepl("_cor$", linetypes), "solid", ifelse(linetypes == "cor", "longdash", "longdash")), linetypes)
        edgelinetypemap <- linetype_map[linetypes]
      }else{
        edgelinetypemap<-NULL
      }
    }else if(plot_type %in% c("control_subnetwork","case_subnetwork","diff_subnetwork")){
      
      edgelinetypemaplist<-lapply(igraphlist,function(x){
        if(length(unique(igraph::E(x)$multiplex_status))>0){
          linetypes<-unique(igraph::E(x)$multiplex_status)
          linetype_map <-setNames(ifelse(grepl("_cor$", linetypes), "solid", ifelse(linetypes == "cor", "longdash", "longdash")), linetypes)
          linetype_map[linetypes]
        }else{NULL}
      })
      names(edgelinetypemaplist)<-net_name
    }else if(plot_type %in% c("differential_network","differential_subnetwork")){
      edgelinetypemaplist<-lapply(igraphlist,function(x){
        if(length(unique(igraph::E(x)$multiplex_status))>0){
          linetypes<-unique(igraph::E(x)$multiplex_status)
          linetype_map <-setNames(ifelse(grepl("_cor$", linetypes), "solid", ifelse(linetypes == "cor", "longdash", "longdash")), linetypes)
          linetype_map[linetypes]
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
        #igraph::V(igraph)$Log2FC<-log2(igraph::V(igraph)$FC)
      }else if(color_type=="Normalized mean"){
        colormap <-grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))
        # if(plot_type %in% c("case_overall_network","case_multi_network","case_subnetwork")){
        #   igraph::V(igraph)$`Normalized mean`<-igraph::V(igraph)$case_norm_mean
        # }else if(plot_type %in% c("control_overall_network","control_multi_network","control_subnetwork")){
        #   igraph::V(igraph)$`Normalized mean`<-igraph::V(igraph)$control_norm_mean
        # }
      }
      if(plot_type %in% c("sub_network","diff_subnetwork","control_subnetwork","case_subnetwork","differential_subnetwork")){
        LAYOUT<-LAYOUT_list[[network_name]]
        if(exists("input_layout")){
          if(!is.null(input_layout)){#########################test
            LAYOUT=input_layout
            LAYOUT_list[[network_name]]<-LAYOUT
          }
        }

      }else{
        if(exists("input_layout")){
        if(!is.null(input_layout)){#########################test
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
        
        #将线形和线颜色图例结合显示
        igraph::E(igraph)$temp_cor_status<-ifelse(igraph::E(igraph)$cor_status %in% focus_name1, igraph::E(igraph)$cor_status, "Other")
        diff_edge_status <- c(Only_in_case,Enhanced_in_case, Only_in_control,Enhanced_in_control, "Conflict relation","Non-significant","Neither group","Other")
        #diff_cor_status_colormap <- setNames(c("#d31b16","#FF00FF","#228b22", "#02fda4", "#87CEFA","#E0E08C" ,"#cfc0bb","gray"), diff_edge_status)
        diff_cor_status_colormap <- setNames(c("#d31b16","#FF00FF","#0000FF", "#00FFFF","#f28500","#E0E08C" ,"#cfc0bb","gray"), diff_edge_status)
        
        #指定focus节点
        related_nodes <- unique(unlist(igraph::ends(igraph, which(igraph::E(igraph)$cor_status %in% focus_name1))))
        # 创建一个逻辑向量，指示哪些节点是相关的
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
      plot <- run_ggraph_plot(igraph=igraph,colormap=colormap,edgecolormap=edgecolormap, color_type =color_type,
                              plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                              plot_layout=LAYOUT,show_node_legend=show_node_legend,show_edge_legend=show_edge_legend,alpha=alpha,
                              edgelinetypemap=edgelinetypemap,edge_color_type=edge_color_type,case_control=case_control) 
      plot <-plot + ggplot2::xlim(xlim_min, xlim_max)+ 
        ggplot2::ylim(ylim_min, ylim_max) 

      
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
          ggplot2::theme(plot.background = ggplot2::element_rect(color = "black"))#图边框
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
    
    ############################增加富集的标注
    if(add_enrichement){
      if(plot_type %in% c("overall_cluster_network")){
        if(is(Network, "Stable_SubNetwork")) {
          other<-run_add_annotation(Network=Network,plot_layout=LAYOUT)
          
          # if(datasave){
          #   outdir1<-file.path(outdir,Networkclass,"NetworkClustering")
          #   if (!dir.exists(outdir1)) {
          #     dir.create(outdir1, recursive = TRUE)
          #   }
          #   readr::write_delim(add_anno,file.path(outdir1,paste0("add_anno_",group_name,".txt")),delim="\t")
          # }
          plotlist <- lapply(plotlist, function(plot) {
            plot + ggplot2::geom_text(
              data = other,
              ggplot2::aes(x = x, y = y, label = enrichment),
              size = 4,
              fontface = "bold",
              color = "black"
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
    #显示节点名称
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
      
      #        }
    }
    
    
    #节点悬停显示
    if (interactive){
      # if(plot_type!="enrichment"){
      plotlist <- lapply(plotlist, function(plot) {
        plot <- ggiraph::girafe(ggobj = plot)
      })  
    }
    
    ################################
  }else if(plot_type %in% c("enrichment","differential_subnetwork_enrichment","subnetwork_enrichment")){
    sub_name1=NULL
    focus_name1=NULL
    if(plot_type =="enrichment"){
      if(is(Network, "Stable_SubNetwork") || is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork") ) {
        plotdatalist<-Network@Enrichment@network
      }else if(is(Network, "Enrichment")){
        plotdatalist<-Network@network
      }else{
        stop(paste0("‘enrichment‘ does not applies to the target object :",class(Network)))
      }
    }else if(plot_type =="subnetwork_enrichment"){
      if(is(Network, "Stable_SubNetwork")) {
        plotdatalist<-Network@Subnet_Enrichment@network
      }else if(is(Network, "Enrichment")){
        plotdatalist<-Network@network
      }else{
        stop(paste0("‘subnetwork_enrichment‘ does not applies to the target object :",class(Network)))
      }
    }else if(plot_type =="differential_subnetwork_enrichment"){
      if(is(Network, "Stable_DifferentialNetwork") || is(Network, "Stable_MultiplexNetwork") ) {
        plotdatalist<-Network@DiffSubnet_Enrichment@network
      }else if(is(Network, "Enrichment")){
        plotdatalist<-Network@network
      }else{
        stop(paste0("‘differential_subnetwork_enrichment‘ does not applies to the target object :",class(Network)))
      }
    }
    if(plot_type =="enrichment"){
      # if(all(grepl("^subnet", names(plotdatalist)))){
      #   enrichplot_type="subnet"
      #   if("all" %in% subnetwork_name){
      #     sub_name1<-names(plotdatalist)
      #   }else{
      #     sub_name1<- ifelse(grepl("^subnet_", subnetwork_name), subnetwork_name, paste0("subnet_", subnetwork_name))#paste0("subnet_",subnetwork_name) 
      #   }
      if(all(grepl("^network", names(plotdatalist)))){#条件特异性总网络富集plotdatalist$network$annotations_filter
        enrichplot_type="network"
        sub_name1="network"
      }else{#差异总网络富集plotdatalist$annotations_filter
        enrichplot_type="diffnet"
        sub_name1="net"
        if("all" %in% focus){
          focus_name1<-unique(plotdatalist$annotations_filter$cor_status)
        }else{
          focus_name1<-focus
        }
        plotdatalist<-list(net=plotdatalist)
      }
    }else if(plot_type =="subnetwork_enrichment"){
      if(all(grepl("^subnet", names(plotdatalist)))){
        enrichplot_type="subnet"
        if("all" %in% subnetwork_name){
          sub_name1<-names(plotdatalist)
        }else{
          sub_name1<- ifelse(grepl("^subnet_", subnetwork_name), subnetwork_name, paste0("subnet_", subnetwork_name))#paste0("subnet_",subnetwork_name) 
        }
        
      }
    }else if(plot_type =="differential_subnetwork_enrichment"){
      enrichplot_type="diffsubnet"
      if("all" %in% subnetwork_name){
        sub_name1<-names(plotdatalist)
      }else{
        sub_name1<-ifelse(grepl("^subnet_", subnetwork_name), subnetwork_name, paste0("subnet_", subnetwork_name))#paste0("subnet_",subnetwork_name) 
      }
      if("all" %in% focus){
        focus_name1<- c(Only_in_case,Enhanced_in_case, Only_in_control,Enhanced_in_control, "Conflict relation","Non-significant","Neither group","Other")
      }else{
        focus_name1<-focus
      }
    }

  annotations_filterlist <- lapply(sub_name1, function(sub_name1) {
    tryCatch({
      plotdata <- as.data.frame(plotdatalist[[sub_name1]]$annotations_filter)
      if(enrichplot_type %in% c("diffnet","diffsubnet")){
        plotdata<- plotdata |>
          dplyr::filter(cor_status %in% focus_name1)
      }
    }, error = function(e) {
      message("Error in plot: ", e$message)
      plotdata <- NULL  # 出错时返回 NULL
    })
    
    if (!is.null(plotdata) && nrow(plotdata) > 0) {
      # if (datasave) {
      #   if(sub_name1=="net"){
      #     outdir1<-file.path(outdir,Networkclass,"OverallNetwork")  
      #   }else{
      #     outdir1<-file.path(outdir,Networkclass,"SubNetwork",sub_name1)
      #   }
      #   if (!dir.exists(outdir1)) {
      #     dir.create(outdir1, recursive = TRUE)
      #   }
      #   readr::write_delim(plotdata, 
      #                      file.path(outdir1, paste0("BubbleDiagram_", gsub(" ", "_", sub_name1), ".txt")), 
      #                      delim = "\t")
      # }
      return(plotdata)
    } else {
      message(paste("Enrichment plot data is empty:",sub_name1))
      return(NULL)  # 如果为空，返回 NULL
    }
  })
  names(annotations_filterlist)<-sub_name1
    # 移除列表中为 NULL 的项，并保留有效项的名称
    annotations_filterlist <- Filter(Negate(is.null), annotations_filterlist)
    sub_name1<-names(annotations_filterlist)
    
    plotlist <-lapply(sub_name1, function(sub_name1) {
      plotdata<-annotations_filterlist[[sub_name1]]
        run_enrichmet_plot(plotdata,richfactor_threshold=richfactor_threshold,plot_title_size=plot_title_size,
                           axis_title_size=axis_title_size,
                           text_size=text_size,
                           legend_title_size=legend_title_size,
                           legend_text_size=legend_text_size,
                           font_family=font_family) +
          ggplot2::ggtitle(sub_name1)


    })
  }else{
    message(paste("Please select the correct ‘plot_type’ parameter."))
  }
  
  

  plotlist <- na.omit(plotlist)
  plotlist <- plotlist[!sapply(plotlist, is.null)]
  if(length(plotlist)>1){
    #图布局
    if(plot_type=="differential_network"){
      plot<-patchwork::wrap_plots(plotlist, ncol = 3, nrow = 1)
    }else if(plot_type=="differential_subnetwork"){
      n_rows <- ceiling(length(plotlist) / 3)
      plot<-patchwork::wrap_plots(plotlist, ncol = 3, nrow = n_rows)
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
  

    if (length(plotlist) > 1) {
      if (plot_type == "differential_network") {
        fig_picture_width <- picture_width * 3
        fig_picture_height <- picture_height
      } else if (plot_type == "differential_subnetwork") {
        n_rows <- ceiling(length(plotlist) / 3)
        fig_picture_width <- picture_width * 3
        fig_picture_height <- picture_height * n_rows
      } else {
        npar <- ceiling(sqrt(length(plotlist)))
        fig_picture_width <- picture_width * npar
        fig_picture_height <- picture_height * npar

      }

    } else if (length(plotlist) == 1) {
      # 单个图的尺寸
      fig_picture_width <- picture_width
      fig_picture_height <- picture_height
    }
  plot
  if(!(plot_type %in% c("enrichment","differential_subnetwork_enrichment","subnetwork_enrichment"))){
    if(exists("net_list")){
      net_list <- replace_colnames_recursive(net_list, casename, controlname)
    }else{
      net_list=list(net=list(nodes=nodes,edges=edges,plot_layout=LAYOUT))
      net_list <- replace_colnames_recursive(net_list, casename, controlname)

    }
    return(new("Network_show",
               data=net_list,
               igraph=igraphlist,
               plot = plot,
               picture_width=fig_picture_width,
               picture_height=fig_picture_height,
               other=other
    )) 
  }else if(plot_type %in% c("enrichment","differential_subnetwork_enrichment","subnetwork_enrichment")){
    return(new("Enrichment_show",
               data=annotations_filterlist,
               plot = plot,
               picture_width=fig_picture_width,
               picture_height=fig_picture_height
    ))
  }
  #return(plot)
}

#' Create enrichment plot
#'
#' This function creates an enrichment plot from enrichment data.
#'
#' @param plotdata Data frame containing enrichment data
#' @param richfactor_threshold Threshold for rich factor (default: 0)
#' @param plot_title_size Size of plot title (default: 12)
#' @param axis_title_size Size of axis titles (default: 8)
#' @param text_size Size of text elements (default: 8)
#' @param legend_title_size Size of legend title (default: 8)
#' @param legend_text_size Size of legend text (default: 8)
#' @param font_family Font family for text (default: "Arial")
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' example_data <- create_example_network_data()
#' enrich_plot <- run_enrichmet_plot(plotdata = example_data$enrichment_data)
#' print(enrich_plot)
run_enrichmet_plot <- function(plotdata=NULL,richfactor_threshold=0,
                               plot_title_size=12,#图标题的大小
                               axis_title_size=8,#轴标题大小
                               text_size=8,#轴刻度文本大小
                               legend_title_size=8,#图例标题大小
                               legend_text_size=8,#图例文字大小
                               font_family="Arial"){
  #top20
  if("cor_status" %in% colnames(plotdata)){
    # 按Types拆分为列表
    data_list <- split(plotdata, plotdata$Types)
    # 处理每个数据集
    processed_list <- purrr::map(data_list, function(df) {
      current_n_cor <- length(unique(df$cor_status))
      # 确定每组应取的行数
      top_n_value <- ifelse(current_n_cor < 4, 15, 5)
      # 按cor_status分组处理
      df <-df |>
      dplyr::filter(RichFactor > richfactor_threshold) |>
      dplyr::group_by(cor_status) |>
      dplyr::arrange(Pvalue, .by_group = TRUE) |>
      dplyr::slice_head(n = top_n_value) |>
      dplyr::ungroup()
      df
    })
    # 合并所有处理后的数据集
    final_data <- dplyr::bind_rows(processed_list)
    # 获取唯一通路列表
    TopPathway <- unique(final_data$Pathway)
    x="cor_status"
    size="RichFactor"
  }else{
    TopPathway<-plotdata |>
      dplyr::filter(RichFactor > richfactor_threshold) |>
      dplyr::group_by(Types) |>#组学
      dplyr::arrange(Pvalue, .by_group = TRUE) |>
      dplyr::distinct(Pathway, .keep_all = FALSE) |>
      dplyr::slice_head(n = 30) |>
      dplyr::pull(Pathway)
    x="RichFactor"
    size="Count"
  }
  plotdata<-plotdata |>
    dplyr::filter(Pathway %in% TopPathway)
  
  plot <- ggplot2::ggplot(plotdata)

#  if("Types" %in% colnames(plotdata)){#多组学
#    plot <- plot + ggplot2::geom_point(ggplot2::aes(x=!!rlang::sym(x),y=Pathway,
#                                                    size=!!rlang::sym(size),colour=Pvalue, shape = Types))+
#      ggplot2::scale_shape(guide = ggplot2::guide_legend(order = 1))#图例顺序
#  }else{
    plot <- plot + ggplot2::geom_point(ggplot2::aes(x=!!rlang::sym(x),y=Pathway,size=!!rlang::sym(size),colour=Pvalue))
    plot <- plot + ggplot2::guides(
    shape = "none"  # 隐藏 shape 图例
  )

#  }
  if(!("cor_status" %in% colnames(plotdata))){
  maxnum <- max(plotdata$Count)
  if(maxnum < 5){
    plot <- plot + ggplot2::scale_size_continuous(limits=c(1,5),range=c(1,5),guide = ggplot2::guide_legend(order = 1))
  }else{
    plot <- plot + ggplot2::scale_size_continuous(range=c(2,8),guide = ggplot2::guide_legend(order = 1))
  }
  }else{
    plot <- plot + ggplot2::scale_size_continuous(limits=c(richfactor_threshold,1),range=c(1,4),guide = ggplot2::guide_legend(order = 1))
  }
  #颜色
  pvaluenum <- max(plotdata$Pvalue)
  plot <- plot + ggplot2::scale_colour_continuous(low="red",high="blue",limits=c(0,pvaluenum),
                                                  guide = ggplot2::guide_colorbar(order = 3))#图例顺序
  ##轴标题
    plot <- plot + ggplot2::theme_bw() 
    plot <- plot + 
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = text_size,family = font_family),#8
        legend.title=ggplot2::element_text(size=legend_title_size,family = font_family),
        legend.text=ggplot2::element_text(size=legend_text_size,family = font_family),
        plot.title = ggplot2::element_text(hjust = 0.5,size = plot_title_size,family = font_family)
      )
    #if(!("cor_status" %in% colnames(plotdata))){
      # plot <- plot +  ggplot2::theme(
      #   axis.title = ggplot2::element_text(size = 12),#12
      #   axis.text = ggplot2::element_text(size = 12),
      #   # legend.title=ggplot2::element_text(size=12),
      #   # legend.text=ggplot2::element_text(size=12),
      #   # axis.text.y = ggplot2::element_text(size=12,angle=0),
      #   plot.title = ggplot2::element_text(hjust = 0.5,size = plot_title_size))
    if("cor_status" %in% colnames(plotdata)){#8
      plot <- plot + 
        ggplot2::labs(x="Differential Correlation Classes")+
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle=45,hjust = 1),
        )
    }

    if(axis_title_size == 0){
      plot <- plot + ggplot2::theme(axis.title = ggplot2::element_blank())
    }else{
      plot <- plot + ggplot2::theme(axis.title = ggplot2::element_text(size = axis_title_size, face="bold", family = font_family))
    }
    

  return(plot)
}

#' Create example network data for demonstration
#'
#' This function creates example network data for testing and demonstration purposes.
#'
#' @return A list containing example network data
#' @export
#'
#' @examples
#' example_data <- create_example_network_data()
#' print(example_data)
create_example_network_data <- function() {
  # Create example nodes
  nodes <- data.frame(
    node = c("A", "B", "C", "D", "E"),
    Class = c("Type1", "Type2", "Type1", "Type3", "Type2"),
    feature_name = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE"),
    omics_name = c("Transcriptomics", "Transcriptomics", "Proteomics", "Metabolomics", "Transcriptomics"),
    case_norm_mean = c(0.8, 1.2, 0.9, 1.1, 0.7),
    control_norm_mean = c(0.7, 1.0, 0.8, 1.0, 0.6),
    FC = c(1.14, 1.2, 1.13, 1.1, 1.17),
    Log2FC = c(0.19, 0.26, 0.18, 0.14, 0.23),
    membership = c("subnet_1", "subnet_1", "subnet_2", "subnet_2", "subnet_1")
  )
  
  # Create example edges
  edges <- data.frame(
    from = c("A", "B", "C", "D", "A"),
    to = c("B", "C", "D", "E", "E"),
    cor = c(0.8, 0.7, -0.6, 0.9, 0.85),
    cor_status = c("positive", "positive", "negative", "positive", "positive"),
    multiplex_status = c("cor", "cor", "cor", "cor", "cor")
  )
  
  # Create example enrichment data
  enrichment_data <- data.frame(
    Pathway = c("Pathway1", "Pathway2", "Pathway3"),
    Pvalue = c(0.001, 0.01, 0.05),
    RichFactor = c(0.5, 0.3, 0.2),
    Count = c(10, 8, 5),
    Types = c("KEGG", "GO", "KEGG"),
    cor_status = c("Enhanced in case", "Only in control", "Enhanced in control")
  )
  
  return(list(
    nodes = nodes,
    edges = edges,
    enrichment_data = enrichment_data
  ))
}

#' Create igraph object from nodes and edges
#'
#' This function creates an igraph object from node and edge data frames.
#'
#' @param nodes Data frame containing node information
#' @param edges Data frame containing edge information
#'
#' @return An igraph object
#' @export
#'
#' @examples
#' example_data <- create_example_network_data()
#' igraph_obj <- run_igraph(nodes = example_data$nodes, edges = example_data$edges)
#' print(igraph_obj)
run_igraph <- function(nodes=NULL,edges=NULL){
  net1 <- igraph::graph_from_data_frame(d=as.data.frame(edges),vertices=as.data.frame(nodes),directed = F)
  ##边线颜色\透明度按相关性
  # E(net1)$color <- ifelse(E(net1)$cor > 0,"red", "green")
  igraph::E(net1)$color <- ifelse(is.na(igraph::E(net1)$cor), "other",#"gray",
                                  ifelse(igraph::E(net1)$cor > 0,"positive","negative"))# "red", "green"))
  igraph::E(net1)$edge_width<-ifelse(is.na(igraph::E(net1)$cor), 0.5, igraph::E(net1)$cor)
  igraph::E(net1)$edge_width<-abs(igraph::E(net1)$edge_width)#不取正数的话负相关会透明
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

#' Create ggraph plot
#'
#' This function creates a network visualization using ggraph.
#'
#' @param igraph An igraph object
#' @param colormap Color mapping for nodes
#' @param edgecolormap Color mapping for edges
#' @param color_type Type of node coloring
#' @param plot_layout Layout for the plot
#' @param edge_color_type Type of edge coloring
#' @param edgelinetypemap Line type mapping for edges
#' @param plot_title_size Size of plot title
#' @param font_family Font family for text
#' @param legend_title_size Size of legend title
#' @param legend_text_size Size of legend text
#' @param alpha Transparency level
#' @param show_edge_legend Whether to show edge legend
#' @param show_node_legend Whether to show node legend
#' @param case_control Case-control label
#'
#' @return A ggraph plot
#' @export
#'
#' @examples
#' example_data <- create_example_network_data()
#' igraph_obj <- run_igraph(nodes = example_data$nodes, edges = example_data$edges)
#' plot <- run_ggraph_plot(igraph = igraph_obj,
#' edgecolormap=setNames(c("red", "blue"), c("positive", "negative")),colormap=c("red","blue","orange"))
#' print(plot)
run_ggraph_plot <- function(igraph=NULL, colormap=NULL,edgecolormap=NULL, color_type = "Class",plot_layout="fr",
                            edge_color_type="color",edgelinetypemap=NULL,plot_title_size=1,
                            font_family="Arial",legend_title_size=8,legend_text_size=8,
                            alpha = 1,show_edge_legend=FALSE,show_node_legend=FALSE,case_control="(case:control)") {
  plot <- ggraph::ggraph(igraph, layout = plot_layout)
  ##线形
  if(all(!is.null(edgelinetypemap))){
    plot <- plot + ggraph::geom_edge_link(
      ggplot2::aes(edge_color =!!rlang::sym(edge_color_type), edge_width = edge_width, linetype = multiplex_status),
      alpha = alpha,
      show.legend = c(edge_color = show_edge_legend, edge_width = FALSE, linetype = TRUE)#, size = FALSE)
    ) +
      ggraph::scale_edge_linetype_manual(values = edgelinetypemap, name = "Type",
                                         guide = ggplot2::guide_legend(order=1)) #图例顺序
    ##############################
  }else{
    plot <- plot + ggraph::geom_edge_link(
      ggplot2::aes(
        edge_color =!!rlang::sym(edge_color_type), edge_width = edge_width),
      alpha = alpha,
      show.legend = c(edge_color = show_edge_legend, edge_width = FALSE, linetype = FALSE)#, size = FALSE
      #show.legend = FALSE
    )
  }
  #边线颜色调整
  if(edge_color_type=="color"){
    plot <- plot + ggraph::scale_edge_color_manual(values = edgecolormap, name = "Correlation",
                                                   guide = ggplot2::guide_legend(order=2))#"Line colour")
  }else{
    plot <- plot + ggraph::scale_edge_color_manual(values = edgecolormap, name = "Correlation status",
                                                   guide = ggplot2::guide_legend(order=2))#"Line colour")
  }
  ######################################
  
  
  #ggraph::scale_edge_alpha_identity() +  # 使用自定义的透明度
  # ggraph::scale_edge_linetype_identity() +  # 使用自定义的线型
  
  #节点颜色大小
  omics_num<-length(unique(igraph::V(igraph)$omics_name))
  if(omics_num>1){
    nodeshape<-c(21,24,22,23,25)#圆形，上三角，方形，菱形，下三角
    # 获取唯一的 omics_name 并按字母排序
    omics_levels <- sort(unique(igraph::V(igraph)$omics_name))
    # 确保 omics_num 不超过 nodeshape 的长度
    omics_num <- min(length(omics_levels), length(nodeshape))
    # 创建一个映射表：omics_name -> shape
    shape_values <- setNames(nodeshape[1:omics_num], omics_levels)
    plot <- plot + ggiraph::geom_point_interactive(#ggraph::geom_node_point
      ggplot2::aes(x=x,y=y,size = size,tooltip = name,data_id = name,#交互
                   fill = !!rlang::sym(color_type),shape=omics_name),#color
      # shape = 21,               # 带填充和边框的点
      colour = "black",
      show.legend = c(size = FALSE, fill = show_node_legend,shape=show_node_legend)#color
    ) +
      ggplot2::scale_shape_manual(name="Omics",values =shape_values) 
  }else{
    plot <- plot + ggiraph::geom_point_interactive(#ggraph::geom_node_point
      ggplot2::aes(x=x,y=y,tooltip = name,data_id = name,#交互
                   size = size, fill = !!rlang::sym(color_type)),#color
      shape = 21,               # 带填充和边框的点
      colour = "black",
      show.legend = c(size = FALSE, fill = show_node_legend)#color
    ) 
  }
  
  #节点颜色调整
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
    #ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = plot_title_size)) +
    ggplot2::theme(
      legend.title=ggplot2::element_text(size=legend_title_size,family = font_family),
      legend.text=ggplot2::element_text(size=legend_text_size,family = font_family),
      plot.title = ggplot2::element_text(hjust = 0.5,size = plot_title_size,family = font_family)
    )+
    ggplot2::coord_equal()

  return(plot)
}

#' Add centrality measures to igraph object and adjust node size
#'
#' This function calculates specified centrality measures and adjusts node sizes accordingly.
#'
#' @param igraph An igraph object
#' @param add_Centrality Character vector specifying centrality measures to add.
#'                      Options: "betweenness", "degree", "eigenvector"
#'
#' @return An igraph object with updated node sizes based on centrality measures
#' @export
#'
#' @examples
#' # Create example network
#' library(igraph)
#' g <- graph_from_data_frame(
#'   data.frame(from = c("A", "B", "C", "D"), 
#'              to = c("B", "C", "D", "A")),
#'   directed = FALSE,
#'   vertices = data.frame(name = c("A", "B", "C", "D"),betweenness=c(0.3,0.4,0.6,0.7))
#' )
#' 
#' # Add betweenness centrality
#' g_with_centrality <- run_add_Centrality(g, "betweenness")
#' 
#' # Check node sizes (now based on betweenness centrality)
#' V(g_with_centrality)$size
run_add_Centrality<-function(igraph=NULL,add_Centrality="betweenness"){
  if(add_Centrality[1]=="betweenness"){
    node_center<-igraph::V(igraph)$betweenness
  }else if(add_Centrality[1]=="degree"){
    node_center<-igraph::V(igraph)$degree
  }else if(add_Centrality[1]=="eigenvector"){
    node_center<-igraph::V(igraph)$eigenvector
  }
  node_size <-node_center# scales::rescale(node_center, to = c(1, 3))  # 调整范围以适应你的需求
  igraph::V(igraph)$size <- node_size
  return(igraph)
}

#' Create centrality scatterplot
#'
#' This function creates a scatterplot comparing different centrality measures
#' for nodes in a network.
#'
#' @param igraph An igraph object
#' @param add_Centrality Character vector specifying centrality measures to plot
#' @param text_size Size of text elements (default: 8)
#' @param legend_title_size Size of legend title (default: 8)
#' @param legend_text_size Size of legend text (default: 8)
#' @param font_family Font family for text (default: "Arial")
#' @param nodes Data frame containing node information
#' @param node_name_type Type of node labels. Options: "id", "name" (default: "id")
#'
#' @return A list containing:
#'   - Centrality_plot: ggplot object of the scatterplot
#'   - Centrality_data: Data frame with centrality values
#' @export
#'
#' @examples
#' # Create example network
#' library(igraph)
#' g <- graph_from_data_frame(
#'   data.frame(from = c("A", "B", "C", "D"), 
#'              to = c("B", "C", "D", "A")),
#'   directed = FALSE,
#'   vertices = data.frame(
#'     name = c("A", "B", "C", "D"),
#'     feature_name = c("GeneA", "GeneB", "GeneC", "GeneD")
#'   )
#' )
#' 
#' # Calculate centrality measures
#' V(g)$betweenness <- betweenness(g)
#' V(g)$degree <- degree(g)
#' V(g)$eigenvector <- eigen_centrality(g)$vector
#' 
#' # Create centrality scatterplot
#' nodes_df <- as_data_frame(g, what = "vertices")
#' centrality_result <- run_plot_Centrality(
#'   igraph = g,
#'   add_Centrality = c("betweenness", "degree"),
#'   nodes = nodes_df,
#'   node_name_type = "name"
#' )
#' 
#' # Display the plot
#' print(centrality_result$Centrality_plot)
#' 
#' # View the data
#' head(centrality_result$Centrality_data)
run_plot_Centrality<-function(igraph=NULL,add_Centrality="betweenness",
                              text_size=8,#轴刻度文本大小
                              legend_title_size=8,#图例标题大小
                              legend_text_size=8,#图例文字大小
                              font_family="Arial",
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
      cols = c("betweenness", "degree", "eigenvector"),  # 或使用 all_of(...) if variable is stored
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
    # ggplot2::theme(
    #   axis.text = ggplot2::element_text(size = 15, face="bold"),
    #   plot.title = ggplot2::element_text(hjust = 0.5))+
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = text_size,family = font_family),#8
      legend.title=ggplot2::element_text(size=legend_title_size,family = font_family),
      legend.text=ggplot2::element_text(size=legend_text_size,family = font_family)
    )+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  return(list(Centrality_plot=Centrality_plot,Centrality_data=Centrality_data))
}

#' Add enrichment annotations to network plot
#'
#' This function adds enrichment pathway annotations to network plots
#' based on the most significant pathway for each subnetwork.
#'
#' @param Network A network object with subnetworks and enrichment data
#' @param plot_layout Layout matrix for the network plot
#'
#' @return A data frame with enrichment annotations and coordinates for plotting
#' @export
#'
#' @examples
#' # This function requires a specific network structure
#' # 
#' 
#' # Add annotations (this would work with a real Network object)
#' # annotations <- run_add_annotation(Network = mock_network, plot_layout = mock_layout)
run_add_annotation<-function(Network=NULL,plot_layout=NULL){
  subnetwork_name1<-names(Network@SubNetwork@subnetworks)
  annotationslist <-lapply(subnetwork_name1, function(subnetwork_name1) {
    annotations_filter<-Network@Subnet_Enrichment@network[[subnetwork_name1]][["annotations_filter"]]
    if(!is.null(annotations_filter) && nrow(annotations_filter)>0){
      annotations_filter<- annotations_filter |>
        dplyr::arrange(Pvalue)
      annotations_filter$Pathway[1]
    }else{
      NA
    }
  })
  names(annotationslist)<-subnetwork_name1
  annotationdata<-tibble::enframe(annotationslist, name = "membership", value = "enrichment") |>
    tidyr::unnest(enrichment)
  # 计算每个 membership 类别的中心坐标（用于标注），排除 "other"
  group_centroids <- Network@SubNetwork@overall_cluster_network[["nodes"]] |>#"kk"
    dplyr::mutate(x=plot_layout[,1],y=plot_layout[,2]) |>
    dplyr::filter(membership != "subnet_other") |>   # 只为非 "other" 组标注
    dplyr::group_by(membership) |>
    dplyr::summarise(x = median(x), y = median(y)) |>
    dplyr::left_join(annotationdata, by = "membership")
  return(group_centroids)
}

#' Recursively replace column names containing case/control with custom names
#'
#' This function recursively traverses a nested list structure and replaces
#' column names containing "case" or "control" with custom names.
#'
#' @param obj A list or data frame object
#' @param casename Custom name to replace "case"
#' @param controlname Custom name to replace "control"
#'
#' @return The object with updated column names
#' @export
#'
#' @examples
#' # Create a nested list with data frames
#' nested_list <- list(
#'   df1 = data.frame(
#'     case_value = c(1, 2, 3),
#'     control_value = c(4, 5, 6),
#'     other_col = c(7, 8, 9)
#'   ),
#'   df2 = data.frame(
#'     case_mean = c(1.1, 2.2),
#'     control_mean = c(3.3, 4.4)
#'   ),
#'   simple_vector = c("a", "b", "c")
#' )
#' 
#' # Replace case/control with custom names
#' updated_list <- replace_colnames_recursive(
#'   nested_list, 
#'   casename = "treatment", 
#'   controlname = "placebo"
#' )
#' 
#' # Check updated column names
#' names(updated_list$df1)
#' names(updated_list$df2)
replace_colnames_recursive <- function(obj, casename, controlname) {
  if (is.list(obj)) {
    # 如果是列表，递归处理每个元素
    lapply(obj, function(x) replace_colnames_recursive(x, casename, controlname))
  } else if (is.data.frame(obj)) {
    # 如果是数据框，替换列名
    colnames(obj) <- gsub("case", casename, colnames(obj))
    colnames(obj) <- gsub("control", controlname, colnames(obj))
    obj
  } else {
    # 其他类型对象保持不变
    obj
  }
}


# Example usage:
# plot_graph(your_igraph, your_color_mapping, color_type = "Class")
#中介度筛选
# run_net_filter<-function(igraph=NULL,cut=0.2,net_filter="betweenness"){
#   #1. 计算中介中心性
#   if(net_filter=="betweenness"){
#     node_center <- igraph::betweenness(
#       graph = igraph,
#       v = igraph::V(igraph),
#       directed = FALSE
#     )
#   # }else if(net_filter=="closeness"){
#   #   node_center <- igraph::closeness(
#   #     graph = igraph,
#   #     vids= igraph::V(igraph)
#   #   )
#   }else if(net_filter=="degree"){
#     node_center <- igraph::degree(
#       graph = igraph,
#       v= igraph::V(igraph)
#     )
#   }else if(net_filter=="eigenvector"){
#     node_center <- igraph::eigen_centrality(
#       graph = igraph
#       , directed = FALSE)$vector
#     
#   }
# 
#   
#   # 2. 将节点按中介度降序排序
#   sorted_nodes <- names(sort(node_center, decreasing = TRUE))
#   
#   # 3. 取前 20% 的节点,最少保留5个
#   top_n <- ceiling(length(sorted_nodes) * cut)
#   if(top_n<5){
#     top_n=5
#   }
#   top_nodes <- sorted_nodes[1:top_n]
#   
#   # 4. 提取只包含这些节点的子图
#   subgraph_top_betweenness <- igraph::induced_subgraph(graph = igraph, vids = top_nodes)
#   return(subgraph_top_betweenness)
# }
###demo:直接做相关性的结果 + 稳定性抽样的结果
# R_threshold=0.8
# color_mapping<-Network@colormapping
# edge_status <- c("positive", "negative", "other")
# edgecolormapping <- setNames(c("red", "green", "gray"), edge_status)
# bootnet_list<-Network@StableNetwork@bootnet_result@bootnet_list
# stable_num<-4#length(Network@StableNetwork@bootnet_result@bootnet_list)
# igraphlist<-lapply(c(1:stable_num),function(stable_number){
#   nodes<-bootnet_list[[stable_number]]$nodes
#   edges<-bootnet_list[[stable_number]]$edges |>
#     dplyr::filter(abs(cor)>R_threshold)
#   run_igraph(nodes,edges)
# })
# bootnetsample=Network@StableNetwork@bootnet_result@bootnet_object$sampleTable
# sample_nodes<-bootnet_list[[1]]$nodes
# sample_edges=data.frame(from=bootnetsample$node1,to=bootnetsample$node2,cor=bootnetsample$value) |>
#   dplyr::filter(from %in% sample_nodes$node) |>
#   dplyr::filter(to %in% sample_nodes$node) |>
#   dplyr::filter(abs(cor)>R_threshold)
#  sample_igraph<-run_igraph(nodes=sample_nodes,edges=sample_edges)
#  LAYOUT<-qgraph::averageLayout(igraphlist,sample_igraph)
# colormap <- color_mapping[unique(igraph::V(sample_igraph)$Class)]
# edgecolormap <- edgecolormapping[unique(igraph::E(sample_igraph)$color)]
#  sampleplot <- run_ggraph_plot(igraph=sample_igraph,colormap=colormap,edgecolormap, color_type = "Class",plot_layout=LAYOUT
#                                ,show_edge_legend = TRUE,show_node_legend = TRUE)
# plotlist<-list()
# for (boot_num in 1:stable_num) {
#   igraph<-igraphlist[[boot_num]]
#   colormap <- color_mapping[unique(igraph::V(igraph)$Class)]
#   edgecolormap <- edgecolormapping[unique(igraph::E(igraph)$color)]
#   plotlist[[boot_num]] <- run_ggraph_plot(igraph=igraph,colormap=colormap,edgecolormap, color_type = "Class",plot_layout=LAYOUT) +
#     ggplot2::ggtitle(paste0("boot_",boot_num)) +
#     ggplot2::theme(
#       plot.background = element_rect(color = "black"))
#   # panel.border = element_rect(fill=NA,color="black", size=0.5))
# }
# npar<-ceiling(sqrt(length(plotlist)))
# plot<-patchwork::wrap_plots(plotlist, ncol = 2, nrow = 2)
# ggsave("bootnet.png", plot = plot, width = 8, height = 8, dpi = 300,bg = "white")
# ggsave("directcor.png", plot = sampleplot, width =16, height = 8, dpi = 300,bg = "white")
