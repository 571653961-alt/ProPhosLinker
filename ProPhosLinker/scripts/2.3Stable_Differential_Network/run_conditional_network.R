#' S4 class to store two condition-specific stable networks for comparative analysis.
#'
#' This class holds the comparison label (e.g., "Tumor-vs-Normal") and two network objects:
#' one for the case/experimental group and one for the control group, each of class \code{StableNetwork}.
#'
#' @slot group_name Character string describing the group comparison (formatted as "A-vs-B").
#' @slot network_case Object of class \code{StableNetwork} for the first group (case/experimental).
#' @slot network_control Object of class \code{StableNetwork} for the second group (control).
#' @name Conditional_network-class
#' @rdname Conditional_network-class
#' @exportClass Conditional_network

setClass("Conditional_network", slots = c(
  group_name="character",
  network_case = "ANY",
  network_control = "ANY"
))

#' Construct condition-specific stable correlation networks for two groups.
#'
#' This function splits the input count data by experimental and control groups (defined in \code{samplelist}),
#' then builds a bootstrapped stable correlation network for each group using \code{run_corStability}.
#' It ensures both networks share the same node set (intersection based on edges) and aligns their edge tables
#' to enable downstream differential network analysis.
#'
#' @param count_table A data frame with rows as features (must include column \code{feature_ID}) and columns as samples.
#' @param samplelist A data frame with columns \code{group} and \code{sample}, mapping each sample to its group.
#' @param compare_group A character string in the format "Case:Control" specifying which groups to compare.
#' @param Diff_anno Optional annotation data frame (e.g., differential expression results) used to annotate nodes.
#' @param node_list Optional character vector of feature IDs to restrict analysis to.
#' @param cor_method Correlation method (e.g., "spearman") passed to \code{run_corStability}.
#' @param nBoots Number of bootstrap iterations for stability assessment (default: 50).
#' @param nCores Number of CPU cores for parallel computation; if NULL, auto-detected.
#' @param stability_threshold Maximum CI range for an edge to be considered stable (default: 0.2).
#' @param bootnet_R_threshold Minimum absolute correlation to retain an edge (default: 0).
#'
#' @return An object of class \code{Conditional_network} containing two aligned \code{StableNetwork} objects,
#'         or \code{NULL} if either group fails to produce a valid network.
#'
#' @importFrom dplyr select filter distinct rename left_join
#' @importFrom tidyr drop_na
#' @importFrom tibble as_tibble_col
#' @importFrom stringr str_split
#' @importFrom purrr map_int
#' @import dplyr
#'
#' @examples
#' # Assuming proper count_table, samplelist, and annotation are available
#' # cond_net <- run_conditional_network(count_table, samplelist, compare_group = "T:N")
#'
#' @export

run_conditional_network<-function(count_table=NULL,samplelist=NULL,compare_group=NULL,Diff_anno=NULL,node_list=NULL,
                                  cor_method="spearman",nBoots=50,nCores=NULL,stability_threshold=0.2,bootnet_R_threshold=0){
  #compare_group:Experimental group:control group
  comparison<-gsub(":","-vs-", compare_group)
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
  conditional_count_table<-run_conditional_data(count_table=count_table,samplelist=samplelist,compare_group=compare_group,node_list=node_list)
  count_table_case<-conditional_count_table$count_table_case
  count_table_control<-conditional_count_table$count_table_control
  network_case<-run_corStability(count_table = count_table_case, group_name = cgroup[1], annotation_table = Diff_anno,
  cor_method = cor_method, nBoots = nBoots, nCores = nCores, bootnet_R_threshold=bootnet_R_threshold,
  stability_threshold = stability_threshold,
  node_list = node_list,uniform_layout=TRUE)
  network_control<-run_corStability(count_table = count_table_control, group_name = cgroup[2], annotation_table = Diff_anno,
                                 cor_method = cor_method, nBoots = nBoots, nCores = nCores, bootnet_R_threshold=bootnet_R_threshold,
                                 stability_threshold = stability_threshold,
                                 node_list = node_list,uniform_layout=TRUE)
  
  if(all(is.null(network_control)) || all(is.null(network_case)) ){
  return(NULL)
  }
  control_edge<-network_control@bootnet_result_filter@bootnet_edge
  case_edge<-network_case@bootnet_result_filter@bootnet_edge  
  cor_data_all1<-run_compare_edge(case_edge) |> 
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

  network_control@bootnet_result_filter@bootnet_node<-network_control@bootnet_result_filter@bootnet_node |>
  dplyr::filter(node %in% union(cor_data_all$from,cor_data_all$to))
  network_case@bootnet_result_filter@bootnet_node<-network_case@bootnet_result_filter@bootnet_node |>
    dplyr::filter(node %in% union(cor_data_all$from,cor_data_all$to))
  
  edgeboth<-cor_data_all|>
    dplyr::select(from,to,from_to)
  controledge <-run_compare_edge(network_control@bootnet_result_filter@bootnet_edge) |>
    dplyr::select(-from,-to)
  caseedge <-run_compare_edge(network_case@bootnet_result_filter@bootnet_edge) |>
    dplyr::select(-from,-to)
  network_control@bootnet_result_filter@bootnet_edge<-edgeboth |>
    dplyr::left_join(controledge,by=c("from_to"="from_to")) |>
    dplyr::select(-from_to)
  network_case@bootnet_result_filter@bootnet_edge<-edgeboth |>
    dplyr::left_join(caseedge,by=c("from_to"="from_to")) |>
    dplyr::select(-from_to)
  Conditional_network_result <- new("Conditional_network",
                                    group_name=comparison,
                                    network_case=network_case,
                                    network_control=network_control
  )
 return(Conditional_network_result)
  
}


#' Split count data into case and control subsets based on sample grouping.
#'
#' This helper function validates the input data, extracts samples belonging to each group
#' specified in \code{compare_group} (format: "Group1:Group2"), and returns separate count tables
#' for case and control. It supports filtering by a predefined \code{node_list}.
#'
#' @param count_table A data frame with a \code{feature_ID} column and sample columns.
#' @param samplelist A data frame with \code{group} and \code{sample} columns defining sample-group mapping.
#' @param compare_group A string of the form "Case:Control" indicating the two groups to extract.
#' @param node_list Optional character vector of feature IDs to subset the data.
#'
#' @return A list with two elements: \code{count_table_case} and \code{count_table_control},
#'         each a data frame with \code{feature_ID} and group-specific sample columns.
#'
#' @importFrom dplyr distinct filter select all_of sym
#' @importFrom tidyr drop_na
#' @importFrom tibble as_tibble_col
#' @importFrom stringr str_split
#' @importFrom purrr map_int
#' @import dplyr
#'
#' @keywords internal

run_conditional_data<-function(count_table=NULL,samplelist=NULL,compare_group=NULL,node_list=NULL){
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
  
  # samplelist
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

  if(any(!is.null(node_list))){
    node_table<-count_table |>
      dplyr::filter(feature_ID %in% node_list)  
  }

  count_table_case<-node_table |>
    dplyr::select(feature_ID,all_of(compare_table$sample_list[[1]]))
  count_table_control<-node_table |>
    dplyr::select(feature_ID,all_of(compare_table$sample_list[[2]]))
  return(list(count_table_case=count_table_case,count_table_control=count_table_control))
}
