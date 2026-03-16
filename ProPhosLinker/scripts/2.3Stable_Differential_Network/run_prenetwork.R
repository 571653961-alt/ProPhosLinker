#' Construct a preliminary correlation network from a feature count table and extract top-connected features.
#'
#' This function computes pairwise correlations (Spearman or Pearson) among features (e.g., genes, taxa)
#' across samples in a given group, constructs an undirected correlation network using igraph,
#' filters edges based on correlation strength and p-value thresholds, and then selects the top
#' features by node degree (connectivity) up to a specified limit (`filter_num`). The result is
#' encapsulated in a `PreCor` S4 object for downstream analysis.
#'
#' @param count_table A data frame with rows as features (must include a column named `feature_ID`)
#'                    and columns as samples. Numeric values represent feature abundances or counts.
#' @param group_name A character vector specifying the name(s) of the biological/experimental group(s).
#'                   If multiple names are provided, they will be collapsed into a single label using "-vs-".
#' @param cor_method Correlation method passed to `Hmisc::rcorr()`. Default is `"spearman"`.
#'                   Other options include `"pearson"`.
#' @param R_threshold Minimum absolute correlation coefficient to retain an edge. Default is `0.8`.
#' @param p_threshold Maximum unadjusted p-value to retain an edge. Default is `0.05`.
#' @param filter_num Maximum number of top-connected nodes (features) to retain after network construction.
#'                   Default is `1000`.
#'
#' @return An object of S4 class `PreCor`, containing:
#'   \describe{
#'     \item{group_name}{The formatted group comparison label.}
#'     \item{precor}{A list with: `nodes`, `raw_edges` (all pairwise correlations), `edges` (filtered edges),
#'                   and `node_degrees` (connectivity scores).}
#'     \item{filter_num}{The requested maximum number of top nodes.}
#'     \item{filter_table}{Subset of `count_table` containing only the top-connected features.}
#'   }
#'
#' @importFrom dplyr filter select
#' @importFrom stats p.adjust
#' @importFrom igraph graph_from_data_frame degree
#' @import Hmisc
#' @import igraph
#' @import dplyr
#'
#' @examples
#' # Example count table
#' count_table <- data.frame(
#'   feature_ID = c("GeneA", "GeneB", "GeneC"),
#'   Sample1 = c(10, 20, 30),
#'   Sample2 = c(15, 25, 35),
#'   Sample3 = c(12, 22, 32)
#' )
#' result <- run_prenetwork(count_table, group_name = "Treatment", R_threshold = 0.7, filter_num = 2)
#' print(result@group_name)
#' head(result@filter_table)
#'
#' @export

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