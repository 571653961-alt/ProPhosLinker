#dplyr parallel bootnet
# 定义 StableNetwork 类
 
# 数据准备和验证
# 
# 网络估计和bootstrap抽样
# 
# 稳定性评估（基于置信区间）
# 
# 多重过滤（R值、p值、稳定性）
# 
# 结果整理和对象创建

setClass("StableNetwork",
         slots = c(
           group_name = "character",
           bootnet_result = "ANY",
           bootnet_result_filter = "ANY"
         )
)

# 定义 bootnet_result 和 bootnet_result_filter 的结构
setClass("BootnetResult",
         slots = c(
           bootnet_object = "ANY",  # 假设 bootnet 对象是任意类型
           bootnet_summary = "data.frame",
           bootnet_stable = "data.frame",
           # bootnet_edge = "data.frame",
           # bootnet_node = "data.frame",
           bootnet_list = "list"
         )
)

setClass("BootnetResultFilter",
         slots = c(
           bootnet_summary = "data.frame",
           bootnet_stable = "data.frame",
           bootnet_edge = "data.frame",
           bootnet_node = "data.frame"#,
           # bootnet_list = "list"
         )
)


#相关性网络的边的稳定性检验
run_corStability<-function(count_table=NULL,group_name="data",annotation_table=NULL,
                     cor_method="spearman", nBoots=500, nCores=NULL,bootnet_R_threshold=0,
                     stability_threshold=0.2,p_filter_table=NULL,bootnet_p_threshold=0.05,
                     node_list=NULL,uniform_layout=FALSE){#'arg' should be one of “cor”, “cov”, “cor_auto”, “npn”, “spearman”# R_threshold=NULL,
  # 1. 数据验证：检查数据行数是否足够
  if(nrow(count_table)<2){
    message("Insufficient data rows for stability analysis. Computation terminated.")
    return(NULL)
  }
  # 2. 处理组名
  #######################group_name
  if(length(group_name)>1){
    precor_group_name<-paste(group_name, collapse = "-vs-")
  }else{
    precor_group_name<-group_name
  }
  ##########################
  # 3. 预处理：如果没有提供p值过滤表，先运行预处理网络分析
if(all(is.null(p_filter_table))){
  PreCor <- run_prenetwork(count_table = count_table, group_name = group_name,
                           cor_method = cor_method, R_threshold = 0, p_threshold = bootnet_p_threshold)
  p_filter_table=PreCor@precor$edges |> dplyr::select(from,to,padj)
}

  #######################BootnetResult
  # 4. 准备数据矩阵
    rownames(count_table)<-count_table$feature_ID
    count_matrix<-count_table |>
    dplyr::select(-feature_ID)
    count_matrix<-as.matrix(t(count_matrix))
    # 5. 设置并行计算核心数
    if(is.null(nCores)){
      nCores <- parallel::detectCores() - 1
    }else{
      nCores <- nCores#parallel::detectCores() - 1
    }
    # 6. 估计网络
  Network<-suppressWarnings({bootnet::estimateNetwork(count_matrix,default = "cor",corMethod = cor_method,nonPositiveDefinite = "continue")  })
  #  bootnetlist[[i]] <- bootnet(Network, statistics=c("strength","closeness"), nBoots = nBoots, nCores = 8, type = "case")
  # 7. 进行bootstrap分析
  bootnetresult <- tryCatch({
    suppressMessages({bootnet::bootnet(Network, nBoots = nBoots, nCores = nCores)})    #边抽样
    #更多次数的抽样示例
    #bootnet(Network, nBoots=1000,nCores=8)
  }, error = function(e) {
    warning(paste("Error occurred for bootnet:", e$message))
    NULL
  })
  # 8. 处理bootstrap结果
  bootnet_sta<-bootnetresult$bootTable |>
    dplyr::filter(type=="edge") |>
    dplyr::mutate(num = as.numeric(gsub("boot ", "", name))) |>
    dplyr::arrange(num) |>
    dplyr::select(-num)
  bootnet_sum<- summary(bootnetresult) |>
    dplyr::filter(type=="edge")
  bootnet_sum<-as.data.frame(bootnet_sum)
  bootnet_sum$CIrange<-bootnet_sum$CIupper-bootnet_sum$CIlower # 计算置信区间范围
  
  #R filter
  # 9. 准备bootstrap列表数据
  bootnet_sta_R<- bootnet_sta 
  #  dplyr::filter(abs(value)>R_threshold) ###保留大于R阈值的边的抽样结果
  # bootnet_sum_e0<- bootnet_sum |>
  #   dplyr::select(node1,node2,mean)
  # colnames(bootnet_sum_e0)<-c("from","to","cor")
  # bootnet_sum_n0<-data.frame(node = union(bootnet_sum_e0$from,bootnet_sum_e0$to)) |>
  #   dplyr::left_join(annotation_table,by=c("node"="feature_ID"))
  
  bootnet_list_e0<- bootnet_sta_R |>
    dplyr::select(node1,node2,value)
  colnames(bootnet_list_e0)<-c("from","to","cor")
  bootnet_list_e0<-split(bootnet_list_e0,  factor(bootnet_sta_R$name, levels = unique(bootnet_sta_R$name)))#bootnet_sta_R$name)

  # 10. 创建bootstrap列表（包含节点和边）
  bootnet_list0 <- lapply(bootnet_list_e0, function(edges) {#为统一布局，节点统一为为筛选R之前的所有节点
    # Create nodes dataframe
    nodes <-data.frame(node =union(bootnet_sum$node1, bootnet_sum$node2)) |>#data.frame(node = union(edges$from, edges$to)) |>
      dplyr::left_join(annotation_table, by = c("node" = "feature_ID"))
    # Return a list containing edges and nodes
    list(nodes = nodes,edges = edges)
  })

  #####布局
  # igraphlist<-lapply(c(1:length(bootnet_list0)),function(stable_number){
  #   nodes<-bootnet_list0[[stable_number]]$nodes
  #   edges<-bootnet_list0[[stable_number]]$edges 
  #   run_igraph(nodes,edges)
  # })
  # #####
  # LAYOUT<-qgraph::averageLayout(igraphlist)
  #######################BootnetResultFilter
  #CI阈值,小于该阈值的被视为稳定的边
  # 11. 稳定性过滤：基于置信区间范围
  CI_cut<-stability_threshold
  # CI_cut<-quantile(bootnet_sum$CIrange, 0.25)
  bootnet_sum$CI_Status <- ifelse(bootnet_sum$CIrange < CI_cut, "stable", "Unstable")
  #Filter    
  # 12. 应用过滤条件
  bootnet_sum_f<-bootnet_sum |>
    dplyr::filter(abs(mean)>bootnet_R_threshold) |>
    dplyr::filter(CI_Status =="stable")
  #是否只整理指定特征的抽样表
  # 13. 如果指定了节点列表，进一步过滤
  if(any(!is.null(node_list))){
    bootnet_sum_f<-bootnet_sum_f |>
      dplyr::filter(node1 %in% node_list) |>
      dplyr::filter(node2 %in% node_list)
  }
  # 14. 检查是否有稳定的边
  if(nrow(bootnet_sum_f)==0){
    message("No stable edges, the operation is terminated.")
    return(NULL)
  }

  #p值
  # 15. 合并p值信息
  bootnet_sum_f<- bootnet_sum_f |>
    dplyr::mutate(from=node1,to=node2)
  bootnet_sum_f <-run_compare_edge(bootnet_sum_f) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::select(-from,-to)
  p_edges<-run_compare_edge(p_filter_table) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::select(-from,-to)
  
  bootnet_sum_f <- dplyr::inner_join(bootnet_sum_f, p_edges, by = c("from_to"))
  if(nrow(bootnet_sum_f)==0){
    message("No significant edges, the operation is terminated.")
    return(NULL)
  }
  #######################
  #ppicor_bootnet_s1<-read_delim(file.path(processdir,"cor.txt"))
  # 16. 准备最终的边和节点数据
  bootnet_sum_e<- bootnet_sum_f |>#"CIrange"    "CI_Status"
    dplyr::select(node1,node2,mean,CIrange,padj)
  colnames(bootnet_sum_e)<-c("from","to","cor","CIrange","p_adjust")
  # 17. 处理节点数据（统一布局或过滤布局）
  if(uniform_layout){
    #为统一布局，节点统一为筛选之前的所有节点
    bootnet_sum_n<-data.frame(node =union(bootnet_sum$node1, bootnet_sum$node2))  |>
      dplyr::left_join(annotation_table,by=c("node"="feature_ID"))
    if(any(!is.null(node_list))){
      bootnet_sum_n1<-data.frame(node=node_list) |>
        dplyr::left_join(bootnet_sum_n,by=c("node"="node"))
      bootnet_sum_n=bootnet_sum_n1
    }
  }else{
    bootnet_sum_n<-data.frame(node = union(bootnet_sum_e$from,bootnet_sum_e$to)) |>
      dplyr::left_join(annotation_table,by=c("node"="feature_ID")) 
  }
  

  # 18. 过滤稳定性数据
  bootnet_sta_f<-  bootnet_sta |>
    dplyr::filter(id %in% bootnet_sum_f$id) |>
  dplyr::filter(abs(value)>bootnet_R_threshold) ###保留大于R阈值的边的抽样结果
  # bootnet_list_e<- bootnet_sta_f |>
  #   dplyr::select(node1,node2,value)
  # colnames(bootnet_list_e)<-c("from","to","cor")
  # bootnet_list_e<-split(bootnet_list_e, bootnet_sta_f$name)
  # bootnet_list <- lapply(bootnet_list_e, function(edges) {
  #   # Create nodes dataframe
  #   nodes <- data.frame(node = union(edges$from, edges$to)) |>
  #     dplyr::left_join(annotation_table, by = c("node" = "feature_ID"))
  #   # Return a list containing edges and nodes
  #   list(nodes = nodes,edges = edges)
  # })
  # 19. 创建结果对象
  bootnet_result <- new("BootnetResult",
                        bootnet_object = bootnetresult,  # 示例中使用空列表
                        bootnet_summary = bootnet_sum,
                        bootnet_stable = bootnet_sta,
                        # bootnet_edge = bootnet_sum_e0,
                        # bootnet_node = bootnet_sum_n0,
                        bootnet_list= bootnet_list0
  )
  bootnet_result_filter <- new("BootnetResultFilter",
                               bootnet_summary = bootnet_sum_f,
                               bootnet_stable = bootnet_sta_f,
                               bootnet_edge = bootnet_sum_e,
                               bootnet_node = bootnet_sum_n#,
                               # bootnet_list = bootnet_list
  )
  # 创建一个 StableNetwork 对象
  # 20. 创建最终的稳定性网络对象
  stable_network <- new("StableNetwork",
                        group_name = precor_group_name,
                        bootnet_result = bootnet_result,
                        bootnet_result_filter = bootnet_result_filter
  )
  invisible(stable_network)
}


