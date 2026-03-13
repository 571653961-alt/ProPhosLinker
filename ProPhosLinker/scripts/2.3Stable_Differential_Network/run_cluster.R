#' 子网络类
#'
#' 该类用于存储子网络分析的结果，包括模块度、颜色映射、整体聚类网络和子网络信息。
#'
#' @slot group_name 字符型，分析的组名称
#' @slot modularity 数值型，网络的模块度指数
#' @slot colormapping 颜色映射方案，用于子网络可视化
#' @slot overall_cluster_network 列表，包含整体聚类网络的节点和边信息
#' @slot subnetworks 列表，包含所有识别出的子网络信息
#'
#' @keywords internal
#dplyr igraph bootnet qgraph RColorBrewer
setClass("SubNetwork",
         slots = c(
           group_name= "character",
           modularity="numeric",
           colormapping = "ANY",
          # Conditional_network_layout="ANY",
           overall_cluster_network = "ANY",
           subnetworks = "ANY"
         )
)
#' 执行子网络聚类分析
#'
#' 该函数使用igraph包进行子网络聚类分析，支持快速贪婪算法和Louvain算法。
#' 当最大社区规模超过指定阈值时，会自动调整分辨率参数。
#'
#' @param nodes 数据框，包含节点信息
#' @param edges 数据框，包含边信息
#' @param clustersize 数值型，社区的最大规模阈值，默认30
#'
#' @return 返回聚类结果对象
#'
#' @examples
#' \dontrun{
#' # 运行子网络聚类分析
#' cluster_result <- run_subnet_cluster(
#' nodes = node_data,
#' edges = edge_data,
#' clustersize = 30
#' )
#' }
#'
#' @export
run_subnet_cluster<-function(nodes=NULL,edges=NULL,clustersize=30){
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
      resolution <- resolution + 0.1  # 你可以根据需要调整步长
    }
  }
  #########节点数多的设置为subnet1
  # Get membership vector
  # membership_vec <- igraph::membership(cfg_t1)
  # community_sizes<-table(membership_vec)
  # # Count sizes of communities
  # raw_membership <- names(community_sizes)
  # # Sort community IDs by size (descending) and get the order of community IDs
  # new_membership <- names(sort(community_sizes, decreasing = TRUE))
  # # Create a mapping from old community ID to new label (starting from 1)
  # map=setNames(new_membership,raw_membership)
  # # Map original membership to new labels
  # new_membership <- as.character(map[as.character(membership_vec)])
  # 
  # cfg_t1$membership <- new_membership
  #########
  return(cfg_t1)
}

#' 添加聚类信息到网络
#'
#' 该函数将聚类结果添加到网络数据中，计算模块度，识别子网络，
#' 并为每个子网络分配颜色映射。
#'
#' @param cfg_t1 聚类结果对象
#' @param nodes 数据框，包含节点信息
#' @param edges 数据框，包含边信息
#' @param group_name 字符型，分析组名称
#' @param diffmessage 字符型，差异分析类型标识，可选值："NULL"、"diff"、"multi"
#'
#' @return 返回SubNetwork对象，包含完整的子网络分析结果
#'
#' @examples
#' \dontrun{
#' # 添加聚类信息到网络
#' subnetwork_result <- run_add_cluster(
#' cfg_t1 = cluster_result,
#' nodes = node_data,
#' edges = edge_data,
#' group_name = "Treatment"
#' )
#' }
#'
#' @export
run_add_cluster<-function(cfg_t1=NULL,nodes=NULL,edges=NULL,group_name=NULL,
                          diffmessage="NULL"#,Conditional_network_layout=NULL#bootnet_list=NULL
                          #module_select="NULL"#非必须参数
){
  #######################group_name
  if(length(group_name)>1){
    precor_group_name<-paste(group_name, collapse = "-vs-")
  }else{
    precor_group_name<-gsub(":","-vs-", group_name)#diff
  }
  #######################modularity:模块度（Q值）衡量网络中社区内部连接的密集程度与随机情况下的差异
  #通常介于**-0.5到1之间，但实际中多为0.3~0.7**（显著社区结构），接近0表示划分与随机无异，负值则表明划分不合理。
  # 高Q值（如>0.3）：社区结构较明显，划分合理。
  # 低Q值（如<0.3）：可能网络本身社区结构模糊，或划分算法未找到最优解。
  # 负Q值：当前划分比随机分配更差，需检查算法或数据
  modularity <- igraph::modularity(cfg_t1) # 网络的模块指数:0.6982418
  ########################membership
  cor_igraph<-run_igraph(nodes,edges)
  t_modules_1 <- sort(table(igraph::membership(cfg_t1)),decr=T)#table(igraph::membership(cfg_t1))#
  t_modules<-table(igraph::membership(cfg_t1))
  # 创建数据
  # old_membership = as.character(igraph::membership(cfg_t1))
  # membership_map<-setNames(names(t_modules),names(t_modules_1))
  # membership_arrange<- as.character(membership_map[old_membership])
  membership_df <- data.frame(
    old_membership =as.character(igraph::membership(cfg_t1)),#membership_arrange,
    node = names(igraph::membership(cfg_t1))
  )
  node_names <- names(igraph::V(cor_igraph))
  # 匹配 membership 值
  old_membership_values <- membership_df$old_membership[match(node_names, membership_df$node)]
  # 将社区结果添加到图的节点属性中(社区命名是按节点数从大到小)
  igraph::V(cor_igraph)$membership <- old_membership_values
  
  ###拆分子模块:结果是一个list (社区命名是按节点数从大到小)
  subgraphs1 <- lapply(names(t_modules_1), function(i) {
    modelcut(igraph::subgraph(cor_igraph, cfg_t1$membership == i))####测试一下删除小于5连通组件
  })
  #统计有多少非na的相关性
  sum_cor <- lapply(seq(1,length(subgraphs1)), function(i) {
    if(diffmessage %in% c("diff","multi")){
      cor_case <- igraph::E(subgraphs1[[i]])$cor_case
      cor_control <- igraph::E(subgraphs1[[i]])$cor_control
      # 统计至少有一个不是 NA 的数量
      sum(!is.na(cor_case) | !is.na(cor_control))
    }else{
      sum(!is.na(igraph::E(subgraphs1[[i]])$cor))
    }
   
  })
  sum_cor<-unlist(sum_cor)
  #统计每个模块有多少节点
  sum_node <- lapply(seq(1,length(subgraphs1)), function(i) {
    sum(length(names(igraph::V(subgraphs1[[i]]))))
  })
  sum_node<-unlist(sum_node)
  ########################membershipfilter:子模块至少包含1个cor且至少有五个节点,否则membership为other
  # 找到至少包含1个cor且至少有五个节点的模块的位置,注意这是已经按代谢物数量排序后的
  non_zero_positions <- which(sum_cor != 0 & sum_node>4)
  if(length(non_zero_positions)==0){
    return(NULL)
    message("No subnetwork containing at least 5 nodes, the network clustering analysis is terminated.")
  }else{
  # if(module_select=="top1"){
  #   non_zero_positions<-non_zero_positions[1]
  # }
  # if(module_select=="top3"){
  #   ## 查看包含代谢物最多的3个模块,且至少有一条相关关系
  #   # 获取前三个至少包含1个cor的模块的位置,注意这是已经按代谢物数量排序后的
  #   selected_positions <- non_zero_positions[1:3]
  #   # 根据选中的位置提取子集
  #   subgraphs1 <- subgraphs1[selected_positions]
  #   #获取对应的原始模块位置(还是在大图中)
  #   t_modules_top3_1 <- t_modules_1[selected_positions]
  #   cor_igraph_top3<-cor_igraph
  #   membership_df<-membership_df |>
  #     mutate(membership_top3=ifelse(membership %in% names(t_modules_top3_1),membership,"other"))
  #   membership_values_top3 <- membership_df$membership_top3[match(node_names, membership_df$node)]
  #   V(cor_igraph_top3)$membership <- membership_values_top3
  #   return(list(color_mapping_cluster,cor_igraph,cor_igraph_top3,subgraphs1))
  # }else if(module_select=="all"){
  ## 查看至少有一条相关关系的模块
  # 根据选中的位置提取子集
  subgraphs1 <- subgraphs1[non_zero_positions]
  t_modules_top <- t_modules_1[non_zero_positions]
  cor_igraph_top<-cor_igraph
  membership_df<-membership_df |>
    dplyr::mutate(membership_top=ifelse(old_membership %in% names(t_modules_top),old_membership,"other"))
  # membership_values_top <- membership_df$membership_top[match(node_names, membership_df$node)]
  # V(cor_igraph_top)$membership <- membership_values_top
  ######################################overall_cluster_network
  bootnet_cluster_edges<-edges
  bootnet_cluster_nodes<-nodes
  ##cluster更名对照表
  membership_mapping<-data.frame(membership_top=as.character(names(t_modules_top)),membership_top_new =as.character(c(1:length(non_zero_positions))))
  membership_df <- membership_df |>
    dplyr::left_join(membership_mapping,by=c("membership_top"="membership_top")) |>
    dplyr::mutate(membership_top_new = dplyr::coalesce(membership_top_new, "other"))
  
  membership_values_top <- membership_df$membership_top_new[match(bootnet_cluster_nodes$node, membership_df$node)]
  bootnet_cluster_nodes$membership<-paste0("subnet_",membership_values_top)
  #######################################subnetworks
  #重命名
  names(subgraphs1)<-paste0("subnet_",membership_mapping$membership_top_new)
  subnet_bootnet_list <- lapply(subgraphs1, function(subgraphs) {
    sub_edge<-igrah_to_edgedata(cor_ppi_igraph_f=subgraphs,diffmessage=diffmessage)
    sub_node<-data.frame(node = union(sub_edge$from,sub_edge$to)) |>
      dplyr::left_join(bootnet_cluster_nodes,by=c("node"="node"))
    x<-sub_node$node
    # Function to filter nodes and edges based on x vector
    # filter_subset <- function(subset, x) {
    #   subset$nodes <- subset$nodes[subset$nodes$node %in% x, ]
    #   subset$edges <- subset$edges[subset$edges$from %in% x & subset$edges$to %in% x, ]
    #   return(subset)
    # }
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
    # Apply the function using lapply
    #filtered_bootnet_list <- lapply(bootnet_list, filter_subset, x = x)
   #,bootnet_list=filtered_bootnet_list)
  })
  ########################color_mapping_cluster
  unique_classes<-as.factor(paste0("subnet_",unique(membership_values_top)))
  unique_classes <- unique_classes[!(unique_classes %in% "subnet_other")]
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  if(length(unique_classes)>length(cluster_Palette)){
    cluster_Palette <-grDevices::colorRampPalette(cluster_Palette)(length(unique_classes))#避免颜色不够
  }
  color_mapping_cluster <- setNames(c(cluster_Palette[1:length(unique_classes)],"grey"), c(as.character(unique_classes),"subnet_other"))
  ##################################
  SubNetwork<-new("SubNetwork",
                  group_name=precor_group_name,
                  modularity=modularity,
                  colormapping = color_mapping_cluster,
                  #Conditional_network_layout=Conditional_network_layout,
                  overall_cluster_network = list(
                    nodes=bootnet_cluster_nodes,
                    edges=bootnet_cluster_edges
                  ),
                  subnetworks = subnet_bootnet_list
                  
  )
  return(SubNetwork)
 }
  
} 

#' 从igraph对象提取边数据
#'
#' 该辅助函数从igraph对象中提取边数据，支持不同类型的网络分析结果。
#'
#' @param cor_ppi_igraph_f igraph对象，包含网络信息
#' @param diffmessage 字符型，差异分析类型标识，可选值："NULL"、"diff"、"multi"
#'
#' @return 返回数据框，包含提取的边信息
#'
#' @keywords internal
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
      p_adjust= igraph::E(cor_ppi_igraph_f)$p_adjust
    ) 
  }
 
  return(edge_data)
}

#' 识别并过滤连通组件
#'
#' 该辅助函数识别网络中的连通组件，并过滤掉规模过小的组件。
#'
#' @param net igraph对象，输入网络
#'
#' @return 返回过滤后的igraph对象，仅保留大规模连通组件
#'
#' @keywords internal
modelcut<-function(net){
  comp <- igraph::components(net)
  # 查看连通组件的大小
  if(max(comp$csize)>5){
    # 移除孤立的节点
    net <- igraph::delete_vertices(net, igraph::degree(net) == 0)
    # 设置阈值，去除小于阈值的连通组件
    threshold <- 5  # 例如，去除小于n个节点的连通组件
    # 过滤连通组件
    large_comp_indices <- which(comp$csize >= threshold)
    # 获取大型连通组件的节点
    large_comp_nodes <- unlist(comp$membership[comp$membership %in% large_comp_indices])
    # 生成新的网络
    g_filtered <- igraph::induced_subgraph(net, names(large_comp_nodes))
    return(g_filtered)
  }else{
    return(net)
  }
}


