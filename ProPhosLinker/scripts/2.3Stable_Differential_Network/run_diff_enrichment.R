
run_diff_enrichment<-function(Differential_network=NULL,species = NULL,nCores=NULL,
                                    omics_name=NULL,enrichment_p_threshold=0.05,
                                    compare_group="data",database_path=NULL,glist_path=NULL){
  
split_edge_data <- split(Differential_network@diff_edges, Differential_network@diff_edges$cor_status)
subnetworks_nodes <- lapply(names(split_edge_data), function(x) {
  data<-split_edge_data[[x]]
  Differential_network@diff_nodes |>
    dplyr::filter(node %in% union(data$from,data$to))
})
names(subnetworks_nodes)<-names(split_edge_data)
Enrichment <- run_enrichment(annotation_table_select = subnetworks_nodes, species = species,nCores=nCores,omics_name=omics_name,
                             group_name=compare_group,enrichment_p_threshold=enrichment_p_threshold,database_path=database_path,
                             glist_path=glist_path)

if(all(is.null(Enrichment))){
  return(NULL)
}
  network=Enrichment@network
if(all(is.null(network))){
  return(NULL)
}

  # Step 1: 给每个数据集添加 cor_status 列
  dfs_by_type <-lapply(names(network), function(type_name) {
    type_data <- network[[type_name]]
    modified_type <- lapply(type_data, function(df) {
      if (any(!is.na(df))) {
        df$cor_status <- type_name
      }
      df
    })
    modified_type  # 返回整个修改后的 type（如 only in case）
  }) |>
    setNames(names(network))  # 恢复外层名称（如 only in case）
  
  # Step 2: 按数据集名合并不同 type 的数据
  all_df_names <- unique(unlist(lapply(dfs_by_type, names)))
  
  combined <- setNames(lapply(all_df_names, function(df_name) {
    do.call(rbind, lapply(dfs_by_type, `[[`, df_name))
  }), all_df_names)

  Enrichment@network<-combined
  return(Enrichment)
}


