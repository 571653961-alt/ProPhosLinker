#' Run pre-network analysis for microbiome/transcriptomics data
#'
#' This function calculates correlations between features from an abundance table,
#' filters edges based on correlation and significance thresholds, constructs a network,
#' and selects top features based on node degree. Results are returned as a PreCor S4 object.
#'
#' @param count_table data.frame, containing a 'feature_ID' column and multiple sample columns
#'                    (abundance values). Rows represent features, columns represent samples.
#' @param group_name character, group label for result identification. If length > 1,
#'                    joined with "-vs-".
#' @param cor_method character, correlation method passed to Hmisc::rcorr. Options:
#'                    "pearson" or "spearman" (default).
#' @param R_threshold numeric, absolute correlation coefficient threshold. Edges with
#'                     |cor| > R_threshold are retained. Default 0.8.
#' @param p_threshold numeric, adjusted P-value threshold. Edges with padj < p_threshold
#'                     are retained. Default 0.05.
#' @param filter_num numeric, maximum number of features to retain based on node degree.
#'                    Default 1000.
#'
#' @return A PreCor S4 object with the following slots:
#'   \item{group_name}{character, group name}
#'   \item{precor}{list, containing nodes (data.frame), edges (data.frame),
#'                 node_degrees (named numeric vector)}
#'   \item{filter_num}{numeric, actual node count limit used for filtering}
#'   \item{filter_table}{data.frame, filtered abundance table with high-degree features
#'                       (includes feature_ID column)}
#'
#' @details 
#' Analysis steps:
#' 1. Transpose count_table to sample × feature matrix and compute correlation matrix
#'    using Hmisc::rcorr.
#' 2. Convert correlation matrix to edge list and adjust P-values using p.adjust
#'    (method = "none").
#' 3. Filter edges by R_threshold and p_threshold, then construct an igraph network.
#' 4. Calculate node degree, select top filter_num features by degree, and return
#'    the filtered abundance table.
#'
#' @importFrom dplyr select filter
#' @importFrom Hmisc rcorr
#' @importFrom igraph graph_from_data_frame degree
#'
#' @examples
#' \dontrun{
#' # Example analysis with count_data containing feature_ID column
#' res <- run_prenetwork(
#'   count_table = count_data,
#'   group_name = "TreatmentA",
#'   cor_method = "spearman",
#'   R_threshold = 0.7,
#'   p_threshold = 0.01,
#'   filter_num = 500
#' )
#' # Check number of edges in network
#' nrow(res@precor$edges)
#' # View filtered feature table
#' head(res@filter_table)
#' }
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
        edges = edges,
        node_degrees=node_degrees
      ),
      filter_num=filter_num,
      filter_table = filter_table

  )
  return(PreCordata)#unique_classes,colors,tits_nodesizes,
}
