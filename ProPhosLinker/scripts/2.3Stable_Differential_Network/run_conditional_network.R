setClass("Conditional_network", slots = c(
  group_name="character",
#  Conditional_network_layout="ANY",
  network_case = "ANY",
  network_control = "ANY"
))

run_conditional_network<-function(count_table=NULL,samplelist=NULL,compare_group=NULL,Diff_anno=NULL,node_list=NULL,
                                  cor_method="spearman",nBoots=50,nCores=NULL,stability_threshold=0.2,bootnet_R_threshold=0){
  #compare_group:Experimental group:control group
  comparison<-gsub(":","-vs-", compare_group)
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
 #提取每组的定量表
  conditional_count_table<-run_conditional_data(count_table=count_table,samplelist=samplelist,compare_group=compare_group,node_list=node_list)
  count_table_case<-conditional_count_table$count_table_case
  count_table_control<-conditional_count_table$count_table_control
  ##稳定性:为了统一布局，保留所有节点
  network_case<-run_corStability(count_table = count_table_case, group_name = cgroup[1], annotation_table = Diff_anno,
  cor_method = cor_method, nBoots = nBoots, nCores = nCores, bootnet_R_threshold=bootnet_R_threshold,
  stability_threshold = stability_threshold,
  node_list = node_list,uniform_layout=TRUE)
  network_control<-run_corStability(count_table = count_table_control, group_name = cgroup[2], annotation_table = Diff_anno,
                                 cor_method = cor_method, nBoots = nBoots, nCores = nCores, bootnet_R_threshold=bootnet_R_threshold,
                                 stability_threshold = stability_threshold,
                                 node_list = node_list,uniform_layout=TRUE)
  
  #只保留至少在一个组里有边的节点
  if(all(is.null(network_control)) || all(is.null(network_case)) ){
  return(NULL)
  }
  control_edge<-network_control@bootnet_result_filter@bootnet_edge
  case_edge<-network_case@bootnet_result_filter@bootnet_edge  
  cor_data_all1<-run_compare_edge(case_edge) |> #run_diff_network.R
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::rename(cor_case=cor) |>
    dplyr::rename(CIrange_case=CIrange) |>
    dplyr::rename(p_adjust_case=p_adjust)
  cor_data_all2 <-run_compare_edge(control_edge) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::rename(cor_control=cor) |>
    dplyr::rename(CIrange_control=CIrange) |>
    dplyr::rename(p_adjust_control=p_adjust)
  cor_data_all <- dplyr::full_join(cor_data_all1, cor_data_all2, by = c("from", "to", "from_to"))
  cor_data_all<-cor_data_all |>
    dplyr::distinct(from_to, .keep_all = TRUE) 
  #更新
  network_control@bootnet_result_filter@bootnet_node<-network_control@bootnet_result_filter@bootnet_node |>
  dplyr::filter(node %in% union(cor_data_all$from,cor_data_all$to))
  network_case@bootnet_result_filter@bootnet_node<-network_case@bootnet_result_filter@bootnet_node |>
    dplyr::filter(node %in% union(cor_data_all$from,cor_data_all$to))
  
  ############取实验组和对照组边的并集，例如实验组中没有但是对照组中有的边在实验组相关性标记为NA
  edgeboth<-cor_data_all|>
    dplyr::select(from,to,from_to)
  controledge <-run_compare_edge(network_control@bootnet_result_filter@bootnet_edge) |>
    dplyr::select(-from,-to)
  caseedge <-run_compare_edge(network_case@bootnet_result_filter@bootnet_edge) |>
    dplyr::select(-from,-to)
  #edgeboth<-data.frame(from_to=union(controledge$from_to,caseedge$from_to)) 
  network_control@bootnet_result_filter@bootnet_edge<-edgeboth |>
    dplyr::left_join(controledge,by=c("from_to"="from_to")) |>
   # dplyr::select(-from,-to) |>
   # tidyr::separate(col = from_to, into = c("from", "to"), sep = "_",remove = FALSE) |>
    dplyr::select(-from_to)
  network_case@bootnet_result_filter@bootnet_edge<-edgeboth |>
    dplyr::left_join(caseedge,by=c("from_to"="from_to")) |>
    #dplyr::select(-from,-to) |>
    #tidyr::separate(col = from_to, into = c("from", "to"), sep = "_",remove = FALSE) |>
    dplyr::select(-from_to)
  ##########统一布局
  # network_case_igraph<- run_igraph(nodes=network_case@bootnet_result_filter@bootnet_node,
  #                                  edges=network_case@bootnet_result_filter@bootnet_edge)
  # network_control_igraph<- run_igraph(nodes=network_control@bootnet_result_filter@bootnet_node,
  #                                     edges=network_control@bootnet_result_filter@bootnet_edge)
  # Conditional_network_layout<-qgraph::averageLayout(network_case_igraph,network_control_igraph)
  ##################################
  Conditional_network_result <- new("Conditional_network",
                                    group_name=comparison,
                                   # Conditional_network_layout=Conditional_network_layout,
                                    network_case=network_case,
                                    network_control=network_control
  )
 return(Conditional_network_result)
  
}

run_conditional_data<-function(count_table=NULL,samplelist=NULL,compare_group=NULL,node_list=NULL){
  ###count_table
  count_table<-as.data.frame(count_table)
  if("feature_ID" %in% colnames(count_table)){
    core_table <- count_table |>
      dplyr::distinct(!!dplyr::sym("feature_ID"), .keep_all = TRUE) |>
      tidyr::drop_na()
  }else{
    stop("No feature_ID column in the count_table, no conditional analysis will be performed.")
  }
  
  #compare_group:Experimental group:control group
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
  
  ##samplelist
  if(all(c("group","sample") %in% colnames(samplelist))){
    compare_table <- tibble::as_tibble_col(cgroup) |>
      dplyr::group_by(.data$value) |> 
      dplyr::mutate(
        sample_list = list(
          samplelist[samplelist$group %in% unlist(stringr::str_split(.data$value, "\\+")), ]$sample
        )
      ) |>
      dplyr::mutate(sample_count = purrr::map_int(.data$sample_list, length)) 
  }else{
    stop("No group column or sample column in the samplelist, no conditional analysis will be performed.")
  }  
  #提取指定特征的定量表
  if(any(!is.null(node_list))){
    node_table<-count_table |>
      dplyr::filter(feature_ID %in% node_list)  
  }
  #提取case和control的定量表
  count_table_case<-node_table |>
    dplyr::select(feature_ID,all_of(compare_table$sample_list[[1]]))
  count_table_control<-node_table |>
    dplyr::select(feature_ID,all_of(compare_table$sample_list[[2]]))
  return(list(count_table_case=count_table_case,count_table_control=count_table_control))
}
