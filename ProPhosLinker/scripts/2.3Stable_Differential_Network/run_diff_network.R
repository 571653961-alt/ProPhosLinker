setClass("Differential_network", slots = c(
  group_name="character",
 # Conditional_network_layout="ANY",
  diff_nodes = "data.frame",
  diff_edges = "data.frame"
))

run_diff_network<-function(Conditional_network=NULL,Conditional_multiplexnetwork=NULL,
                           edge_FC_threshold=1.2,edge_p_threshold=0.05,compare_group=NULL){
  if(any(!is.null(Conditional_multiplexnetwork))){
    case_node<-Conditional_multiplexnetwork@network_case@nodes
    case_edge<-Conditional_multiplexnetwork@network_case@edges  
    control_node<-Conditional_multiplexnetwork@network_control@nodes
    control_edge<-Conditional_multiplexnetwork@network_control@edges 
   # Conditional_network_layout<-NULL
  }else{
    case_node<-Conditional_network@network_case@bootnet_result_filter@bootnet_node
    case_edge<-Conditional_network@network_case@bootnet_result_filter@bootnet_edge  
    control_node<-Conditional_network@network_control@bootnet_result_filter@bootnet_node
    control_edge<-Conditional_network@network_control@bootnet_result_filter@bootnet_edge
   # Conditional_network_layout<-Conditional_network@Conditional_network_layout
  }

  comparison<-gsub(":","-vs-", compare_group)
  split_result <- strsplit(comparison, "-vs-")[[1]]
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
  
  case_bootnet_stable<-Conditional_network@network_case@bootnet_result_filter@bootnet_stable |>
    dplyr::select(node1,node2,value)
  colnames(case_bootnet_stable)<-c("from","to","cor")
  
  control_bootnet_stable<-Conditional_network@network_control@bootnet_result_filter@bootnet_stable |>
    dplyr::select(node1,node2,value)
  colnames(control_bootnet_stable)<-c("from","to","cor")

  ###edges diff
  cor_data_all1<-run_compare_edge(case_edge) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::rename(cor_case=cor) |>
    dplyr::rename(CIrange_case=CIrange) |>
    dplyr::rename(p_adjust_case=p_adjust)
  if("multiplex_status" %in% colnames(cor_data_all1)){
    cor_data_all1 <- cor_data_all1 |>
      dplyr::mutate(
        # 从 multiplex_status 列提取 relationship 值
        relationship = ifelse(
          grepl("_cor$", multiplex_status),  # 检查是否以 "_cor" 结尾
          sub("_cor$", "", multiplex_status),  # 如果是，去掉 "_cor" 后缀
          multiplex_status  # 如果不是，保持原值
        )
      ) |>
      dplyr::select(-multiplex_status)  # 删除 multiplex_status 列
  }
  
  cor_data_all2 <-run_compare_edge(control_edge) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::rename(cor_control=cor) |>
    dplyr::rename(CIrange_control=CIrange) |>
    dplyr::rename(p_adjust_control=p_adjust)
  if("multiplex_status" %in% colnames(cor_data_all2)){
    cor_data_all2 <- cor_data_all2 |>
      dplyr::mutate(
        # 从 multiplex_status 列提取 relationship 值
        relationship = ifelse(
          grepl("_cor$", multiplex_status),  # 检查是否以 "_cor" 结尾
          sub("_cor$", "", multiplex_status),  # 如果是，去掉 "_cor" 后缀
          multiplex_status  # 如果不是，保持原值
        )
      ) |>
      dplyr::select(-multiplex_status)  # 删除 multiplex_status 列
  }
  
  if("score" %in% colnames(cor_data_all1) && "score" %in% colnames(cor_data_all2)){
    cor_data_all <- dplyr::full_join(cor_data_all1, cor_data_all2, by = c("from", "to","score", "from_to","relationship")) |>
      dplyr::mutate(multiplex_status=dplyr::case_when(
        (!is.na(cor_case) | !is.na(cor_control)) & !is.na(score) ~ paste(relationship, "cor", sep = "_"),
        (!is.na(cor_case) | !is.na(cor_control)) & is.na(score) ~ "cor",
        (is.na(cor_case) & is.na(cor_control)) & !is.na(score) ~ relationship
      )) |>
      dplyr::select(-relationship) 

    }else{
      cor_data_all <- dplyr::full_join(cor_data_all1, cor_data_all2, by = c("from", "to", "from_to"))
    }

  cor_data_all<-cor_data_all |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    as.data.frame() |>
    dplyr::mutate(
      case_control=cor_case*cor_control,
      cor_status = dplyr::case_when(
        (!is.na(case_control)) & (case_control<0) ~ "Conflict relation",
        (!is.na(case_control)) & (case_control>0) ~ "Both group",
      (is.na(cor_case)) & (!is.na(cor_control)) ~ Only_in_control,#"Only in control",
      (!is.na(cor_case)) & (is.na(cor_control)) ~ Only_in_case,#"Only in case",
      (is.na(cor_case)) & (is.na(cor_control)) ~ "Neither group"
    )) |>
    dplyr::select(-case_control)

  #如果存在两组共有的边，计算差异
  cor_data_both<-cor_data_all |>
    dplyr::filter(cor_status=="Both group")
  if(nrow(cor_data_both)>0){
    cor_data_boot1 <-run_compare_edge(case_bootnet_stable) |>
      dplyr::filter(from_to  %in% cor_data_both$from_to) 
    cor_data_boot2 <-run_compare_edge(control_bootnet_stable) |>
      dplyr::filter(from_to  %in% cor_data_both$from_to)
    #diff
    results<-lapply(cor_data_both$from_to,function(x){
      casedata<-cor_data_boot1 |>
        dplyr::filter(from_to==x) |>
        dplyr::pull(cor) |>
        abs()
      condata<-cor_data_boot2 |>
        dplyr::filter(from_to==x) |>
        dplyr::pull(cor) |>
        abs()
      ca_mean <- mean(casedata)
      co_mean <- mean(condata)
      fd <- ca_mean/co_mean 
      p <- t.test(log2(casedata),log2(condata))
      k <- c(fd,p$p.value)
      names(k) <- c("cor_FC","cor_p_value") 
      return(k)
    })
    results_df <- do.call(rbind, results)
    results_df<-as.data.frame(results_df) |>   
      dplyr::mutate(cor_status= dplyr::case_when(
        cor_p_value < edge_p_threshold & cor_FC > edge_FC_threshold ~ Enhanced_in_case,#"Enhanced in case",#"Up",
        cor_p_value < edge_p_threshold & cor_FC < 1 / edge_FC_threshold ~Enhanced_in_control,#"Enhanced in control",# "Down",
        TRUE ~ "Non-significant"
      ))
    results_df$from_to<-cor_data_both$from_to
    #整合差异结果
    cor_data_diff <- cor_data_all |>
      dplyr::left_join(results_df, by = "from_to")
    
    # Update cor_status in cor_data_diff where there is a match
    cor_data_diff$cor_status <- ifelse(!is.na(cor_data_diff$cor_status.y), 
                                       cor_data_diff$cor_status.y, 
                                       cor_data_diff$cor_status.x)
  
    # Remove the temporary columns used for joining
    cor_data_diff <- cor_data_diff |>
      dplyr::select(-cor_status.x, -cor_status.y)
    
  }else{
    cor_data_diff<- cor_data_all |>
      dplyr::mutate(cor_status=cor_status,
                    cor_FC = NA,
                    cor_p_value = NA) 
  }
  cor_data_diff<- cor_data_diff |>
    dplyr::mutate(cor = cor_FC)

  ##nodes diff
  merge_node <- rbind(case_node,control_node) |>
    dplyr::distinct(node, .keep_all = TRUE) #|>
  #  dplyr::filter(node %in% union(cor_data_diff$from,cor_data_diff$to)) #只保留至少在一个组里有边的节点
  Differential_network_result <- new("Differential_network",
                                     group_name=comparison,
                                   #  Conditional_network_layout=Conditional_network_layout,
                                     diff_nodes=merge_node,
                                     diff_edges=cor_data_diff
  )
  
  return(Differential_network_result)
}




#为了比对边，进行预处理
# 函数：处理网络边数据，标准化边的表示方式a
run_compare_edge<-function(edges){
  # 1. 标准化边的方向：确保每条边的两个节点按固定顺序排列
  # 使用pmin和pmax创建排序后的节点对，消除边的方向性
  edges <-transform(edges, 
                    sorted_from = pmin(from, to),  # 取两个节点中较小的作为sorted_from
                    sorted_to = pmax(from, to))    # 取两个节点中较大的作为sorted_to
  # 2. 创建唯一的边标识符：使用排序后的节点对创建无向边标识
  # 格式为"较大节点_较小节点"，确保同一条无向边有相同的标识符
  edges$from_to<-paste0(edges$sorted_to,"_",edges$sorted_from)
  # 3. 清理临时列：移除排序过程中创建的临时列
  edges<-edges |>
    dplyr::select(-sorted_from, -sorted_to) 
  # 4. 返回处理后的边数据
  return(edges)
}


