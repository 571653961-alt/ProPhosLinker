#igraph dplyr Hmisc
# Define S4 Classes
setClass("PreCor", slots = c(
  group_name = "character",
  # samplelist = "data.frame",
  # raw_table = "data.frame",
  precor = "list",
  filter_num="numeric",
  filter_table = "data.frame"
))

run_prenetwork<-function(count_table=NULL,group_name=NULL,cor_method="spearman",R_threshold=0.8,p_threshold=0.05,filter_num=1000){
  if(length(group_name)>1){
    precor_group_name<-paste(group_name, collapse = "-vs-")
  }else{
    precor_group_name<-group_name
  }
  ############precor
  cor_table<-count_table
  rownames(cor_table)<-cor_table$feature_ID
  cor_table<-cor_table |>
    dplyr::select(-feature_ID)
  cor_table<-as.matrix(t(cor_table))
  t_its_cor <- Hmisc::rcorr(cor_table,type=cor_method)#
  CorrDF <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      from = rownames(cormat)[col(cormat)[ut]],
      to = rownames(cormat)[row(cormat)[ut]],
      cor =(cormat)[ut],
      p = pmat[ut]
    )
  }
  t_its_cor_df <- CorrDF(t_its_cor$r,t_its_cor$P)
  t_its_cor_df$padj <- p.adjust(t_its_cor_df$p, method = "none") #method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  edges <- t_its_cor_df |>
    dplyr::filter(abs(cor) > R_threshold & abs(cor) < 0.999 & padj < p_threshold)
  nodes <- data.frame(node = union(edges$from,edges$to))
  net <- igraph::graph_from_data_frame(d=edges,vertices=nodes,directed = F)
  node_degrees <- igraph::degree(net, mode = "all")
  # Step 4: Sort nodes by degree and select top 300
  num_node<-length(names(sort(node_degrees, decreasing = TRUE)))
  if(num_node>filter_num){
    top_nodes <- names(sort(node_degrees, decreasing = TRUE))[1:filter_num]
  }else{
    top_nodes <- names(sort(node_degrees, decreasing = TRUE))
  }
  filter_table<-count_table |>
  dplyr::filter(feature_ID %in% top_nodes)
  ################################
  PreCordata<-new("PreCor", 
      group_name = precor_group_name,
      # samplelist = samplelist,
      # raw_table = raw_table,
      precor = list(
        nodes = nodes,
        raw_edges=t_its_cor_df,
        edges = edges,
        node_degrees=node_degrees
      ),
      filter_num=filter_num,
      filter_table = filter_table

  )
  return(PreCordata)#unique_classes,colors,tits_nodesizes,
}
