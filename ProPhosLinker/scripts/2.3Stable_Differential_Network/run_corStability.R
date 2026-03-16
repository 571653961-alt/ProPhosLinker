#' S4 class to store a stable network derived from bootstrapped correlation analysis.
#'
#' This class holds the group name and two main components:
#' the full bootnet result (`bootnet_result`) and its filtered version (`bootnet_result_filter`)
#' based on stability, correlation strength, and significance thresholds.
#'
#' @slot group_name Character string identifying the sample group or condition.
#' @slot bootnet_result Object of class \code{BootnetResult} containing full bootnet output.
#' @slot bootnet_result_filter Object of class \code{BootnetResultFilter} with filtered edges/nodes.
#' @name StableNetwork-class
#' @rdname StableNetwork-class
#' @exportClass StableNetwork

setClass("StableNetwork",
         slots = c(
           group_name = "character",
           bootnet_result = "ANY",
           bootnet_result_filter = "ANY"
         )
)

#' S4 class to encapsulate raw bootnet analysis results.
#'
#' Contains the original \code{bootnet} object, summary statistics, stable edge estimates,
#' and a list of network realizations (edges per bootstrap sample).
#'
#' @slot bootnet_object The raw \code{bootnet} object returned by \code{bootnet::bootnet()}.
#' @slot bootnet_summary Data frame from \code{summary(bootnetresult)} for edges.
#' @slot bootnet_stable Long-format data frame of bootstrap edge estimates (\code{bootTable}).
#' @slot bootnet_list List of networks (edges + nodes) for each bootstrap iteration.
#' @name BootnetResult-class
#' @rdname BootnetResult-class
#' @exportClass BootnetResult
setClass("BootnetResult",
         slots = c(
           bootnet_object = "ANY",
           bootnet_summary = "data.frame",
           bootnet_stable = "data.frame",
           bootnet_list = "list"
         )
)

#' S4 class to store filtered bootnet results after applying stability and significance criteria.
#'
#' Represents the final usable network: only edges that are stable, sufficiently strong,
#' and statistically significant are retained, along with associated node metadata.
#'
#' @slot bootnet_summary Filtered edge summary (stable, |cor| > threshold, significant).
#' @slot bootnet_stable Filtered bootstrap edge estimates matching the summary.
#' @slot bootnet_edge Final edge table with columns: from, to, cor, CIrange, p_adjust.
#' @slot bootnet_node Final node table with annotation merged from external source.
#' @name BootnetResultFilter-class
#' @rdname BootnetResultFilter-class
#' @exportClass BootnetResultFilter
setClass("BootnetResultFilter",
         slots = c(
           bootnet_summary = "data.frame",
           bootnet_stable = "data.frame",
           bootnet_edge = "data.frame",
           bootnet_node = "data.frame"
         )
)


#' Estimate a stable correlation-based network using bootstrapping and filtering.
#'
#' This function constructs a correlation network from a feature-by-sample count table,
#' assesses edge stability via non-parametric bootstrapping using \code{bootnet},
#' and applies multiple filters: correlation magnitude, confidence interval width (stability),
#' and adjusted p-value from a precomputed significance test. It returns an S4 object
#' of class \code{StableNetwork} containing both raw and filtered results.
#'
#' @param count_table A data frame with rows as features (must include column \code{feature_ID})
#'                    and columns as samples. Numeric values represent abundances.
#' @param group_name Character vector specifying the group label(s). If length > 1, collapsed with "-vs-".
#' @param annotation_table Optional data frame to annotate nodes (merged by \code{feature_ID}).
#' @param cor_method Correlation method passed to \code{bootnet::estimateNetwork()} (e.g., "spearman").
#' @param nBoots Number of bootstrap iterations (default: 500).
#' @param nCores Number of CPU cores for parallel bootstrapping. If NULL, uses all but one core.
#' @param bootnet_R_threshold Minimum absolute correlation to retain an edge (default: 0).
#' @param stability_threshold Maximum allowed CI range for an edge to be considered "stable" (default: 0.2).
#' @param p_filter_table Optional data frame with columns \code{from}, \code{to}, \code{padj}
#'                       used to filter edges by significance. If NULL, computed internally via \code{run_prenetwork}.
#' @param bootnet_p_threshold P-value threshold used in internal pre-network if \code{p_filter_table} is NULL.
#' @param node_list Optional character vector: restrict final network to these nodes.
#' @param uniform_layout Logical: if TRUE, include all nodes from the full bootnet summary (not just those in edges)
#'                       when building the node table; useful for consistent layout across groups.
#'
#' @return An invisible object of class \code{StableNetwork}, or \code{NULL} if:
#'         - input has <2 rows,
#'         - no stable edges pass filters,
#'         - or no significant edges remain after p-value filtering.
#'
#' @importFrom dplyr select filter mutate arrange distinct left_join inner_join
#' @importFrom stats setNames
#' @importFrom parallel detectCores
#' @importFrom utils head
#' @import bootnet
#' @import dplyr
#'
#' @examples
#' # Example usage requires count_table and annotation_table
#' # result <- run_corStability(count_table, group_name = "Tumor", cor_method = "spearman")
#'
#' @export

run_corStability<-function(count_table=NULL,group_name="data",annotation_table=NULL,
                     cor_method="spearman", nBoots=500, nCores=NULL,bootnet_R_threshold=0,
                     stability_threshold=0.2,p_filter_table=NULL,bootnet_p_threshold=0.05,
                     node_list=NULL,uniform_layout=FALSE){#'arg' should be one of “cor”, “cov”, “cor_auto”, “npn”, “spearman”# R_threshold=NULL,
  if(nrow(count_table)<2){
    message("Insufficient data rows for stability analysis. Computation terminated.")
    return(NULL)
  }
  
  if(length(group_name)>1){
    precor_group_name<-paste(group_name, collapse = "-vs-")
  }else{
    precor_group_name<-group_name
  }
 
if(all(is.null(p_filter_table))){
  PreCor <- run_prenetwork(count_table = count_table, group_name = group_name,
                           cor_method = cor_method, R_threshold = 0, p_threshold = bootnet_p_threshold)
  p_filter_table=PreCor@precor$edges |> dplyr::select(from,to,padj)
}
    rownames(count_table)<-count_table$feature_ID
    count_matrix<-count_table |>
    dplyr::select(-feature_ID)
    count_matrix<-as.matrix(t(count_matrix))
    if(is.null(nCores)){
      nCores <- parallel::detectCores() - 1
    }else{
      nCores <- nCores#parallel::detectCores() - 1
    }
  Network<-suppressWarnings({bootnet::estimateNetwork(count_matrix,default = "cor",corMethod = cor_method,nonPositiveDefinite = "continue")  })
  bootnetresult <- tryCatch({
    suppressMessages({bootnet::bootnet(Network, nBoots = nBoots, nCores = nCores)})
  }, error = function(e) {
    warning(paste("Error occurred for bootnet:", e$message))
    NULL
  })

  bootnet_sta<-bootnetresult$bootTable |>
    dplyr::filter(type=="edge") |>
    dplyr::mutate(num = as.numeric(gsub("boot ", "", name))) |>
    dplyr::arrange(num) |>
    dplyr::select(-num)
  bootnet_sum<- summary(bootnetresult) |>
    dplyr::filter(type=="edge")
  bootnet_sum<-as.data.frame(bootnet_sum)
  bootnet_sum$CIrange<-bootnet_sum$CIupper-bootnet_sum$CIlower
  
  
  bootnet_sta_R<- bootnet_sta 
  
  bootnet_list_e0<- bootnet_sta_R |>
    dplyr::select(node1,node2,value)
  colnames(bootnet_list_e0)<-c("from","to","cor")
  bootnet_list_e0<-split(bootnet_list_e0,  factor(bootnet_sta_R$name, levels = unique(bootnet_sta_R$name)))

  bootnet_list0 <- lapply(bootnet_list_e0, function(edges) {
    # Create nodes dataframe
    nodes <-data.frame(node =union(bootnet_sum$node1, bootnet_sum$node2)) |>
      dplyr::left_join(annotation_table, by = c("node" = "feature_ID"))
    # Return a list containing edges and nodes
    list(nodes = nodes,edges = edges)
  })

  
  CI_cut<-stability_threshold
  bootnet_sum$CI_Status <- ifelse(bootnet_sum$CIrange < CI_cut, "stable", "Unstable")
  #Filter    
  bootnet_sum_f<-bootnet_sum |>
    dplyr::filter(abs(mean)>bootnet_R_threshold) |>
    dplyr::filter(CI_Status =="stable")

  if(any(!is.null(node_list))){
    bootnet_sum_f<-bootnet_sum_f |>
      dplyr::filter(node1 %in% node_list) |>
      dplyr::filter(node2 %in% node_list)
  }

  if(nrow(bootnet_sum_f)==0){
    message("No stable edges, the operation is terminated.")
    return(NULL)
  }

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
  
  bootnet_sum_e<- bootnet_sum_f |>#"CIrange"    "CI_Status"
    dplyr::select(node1,node2,mean,CIrange,padj)
  colnames(bootnet_sum_e)<-c("from","to","cor","CIrange","p_adjust")
  
  if(uniform_layout){
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
  

  bootnet_sta_f<-  bootnet_sta |>
    dplyr::filter(id %in% bootnet_sum_f$id) |>
  dplyr::filter(abs(value)>bootnet_R_threshold)
  
  bootnet_result <- new("BootnetResult",
                        bootnet_object = bootnetresult,
                        bootnet_summary = bootnet_sum,
                        bootnet_stable = bootnet_sta,
                        bootnet_list= bootnet_list0
  )
  bootnet_result_filter <- new("BootnetResultFilter",
                               bootnet_summary = bootnet_sum_f,
                               bootnet_stable = bootnet_sta_f,
                               bootnet_edge = bootnet_sum_e,
                               bootnet_node = bootnet_sum_n
  )
  
  stable_network <- new("StableNetwork",
                        group_name = precor_group_name,
                        bootnet_result = bootnet_result,
                        bootnet_result_filter = bootnet_result_filter
  )
  invisible(stable_network)
}