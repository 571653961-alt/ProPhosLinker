#' Conditional Network S4 Class
#'
#' An S4 class that stores paired network results for two experimental conditions,
#' enabling direct comparison of network structures between case and control groups.
#'
#' @slot group_name character, comparison label describing the two conditions
#'                   (e.g., "Treatment-vs-Control").
#' @slot network_case ANY, network result object (typically BootnetResult) for
#'                    the experimental/case group.
#' @slot network_control ANY, network result object (typically BootnetResult) for
#'                      the control group.
#'
#' @name Conditional_network-class
#' @rdname Conditional_network-class
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a conditional network object
#' cond_net <- new("Conditional_network",
#'   group_name = "Treatment-vs-Control",
#'   network_case = case_network_result,
#'   network_control = control_network_result
#' )
#' 
#' # Access components
#' cond_net@group_name
#' case_net <- cond_net@network_case
#' control_net <- cond_net@network_control
#' }
setClass("Conditional_network", slots = c(
  group_name="character",
#  Conditional_network_layout="ANY",
  network_case = "ANY",
  network_control = "ANY"
))

#' Build conditional networks comparing two experimental conditions
#'
#' This function constructs and compares correlation networks between experimental
#' and control groups. It performs bootstrap stability analysis on both networks,
#' then aligns edges to enable direct comparison of correlation patterns across
#' conditions.
#'
#' @param count_table data.frame, containing a 'feature_ID' column and abundance
#'                    values for each sample as additional columns.
#' @param samplelist data.frame, sample metadata containing 'sample' and 'group'
#'                   columns for group assignment.
#' @param compare_group character, comparison specification in format
#'                      "Experimental_group:Control_group". Groups can be combined
#'                      using "+" (e.g., "Treatment1+Treatment2:Control").
#' @param Diff_anno data.frame, annotation table for features (e.g., differential
#'                  expression results) used for node coloring and labeling.
#' @param node_list character vector, optional. Subset of features to include in
#'                  network analysis. If NULL, all features are used.
#' @param cor_method character, correlation method. Options: "pearson" or "spearman".
#'                   Default "spearman".
#' @param nBoots integer, number of bootstrap resamples for stability analysis.
#'               Default 50.
#' @param nCores integer, number of CPU cores for parallel processing. If NULL,
#'               detects available cores automatically.
#' @param stability_threshold numeric, threshold for edge stability (0-1). Edges
#'                             with stability below this threshold are filtered out.
#'                             Default 0.2.
#' @param bootnet_R_threshold numeric, correlation threshold for bootnet network
#'                             construction. Default 0.
#'
#' @return A Conditional_network S4 object with the following slots:
#'   \item{group_name}{character, comparison label (e.g., "Treatment-vs-Control")}
#'   \item{network_case}{BootnetResult object, stable network for experimental group}
#'   \item{network_control}{BootnetResult object, stable network for control group}
#'
#' @details
#' Analysis steps:
#' 1. Parse compare_group to identify experimental and control groups.
#' 2. Split count_table into case and control subsets using run_conditional_data.
#' 3. Build stable networks for each group using run_corStability.
#' 4. Extract and compare edges from both networks.
#' 5. Align edge sets to ensure both networks contain identical node pairs.
#' 6. Filter nodes to only those present in aligned edges.
#' 7. Return combined conditional network object.
#'
#' @importFrom dplyr full_join distinct rename select filter left_join
#' @importFrom tidyr separate
#'
#' @examples
#' \dontrun{
#' # Build conditional networks comparing Treatment vs Control
#' cond_network <- run_conditional_network(
#'   count_table = abundance_data,
#'   samplelist = sample_metadata,
#'   compare_group = "Treatment:Control",
#'   Diff_anno = differential_results,
#'   cor_method = "spearman",
#'   nBoots = 100,
#'   stability_threshold = 0.3
#' )
#' 
#' # Access networks for each condition
#' case_network <- cond_network@network_case
#' control_network <- cond_network@network_control
#' }
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
 
  ##################################
  Conditional_network_result <- new("Conditional_network",
                                    group_name=comparison,
                                   # Conditional_network_layout=Conditional_network_layout,
                                    network_case=network_case,
                                    network_control=network_control
  )
 return(Conditional_network_result)
  
}

#' Prepare data for conditional network analysis
#'
#' Internal function that splits the abundance table into case and control subsets
#' based on sample grouping and optional feature filtering.
#'
#' @param count_table data.frame, containing a 'feature_ID' column and abundance
#'                    values for each sample as additional columns.
#' @param samplelist data.frame, sample metadata containing 'sample' and 'group'
#'                   columns for group assignment.
#' @param compare_group character, comparison specification in format
#'                      "Experimental_group:Control_group". Groups can be combined
#'                      using "+".
#' @param node_list character vector, optional. Subset of features to include.
#'                  If NULL or empty, all features are retained.
#'
#' @return A list containing two data.frames:
#'   \item{count_table_case}{data.frame, abundance table for experimental group}
#'   \item{count_table_control}{data.frame, abundance table for control group}
#'
#' @details
#' Processing steps:
#' 1. Validate count_table contains 'feature_ID' column.
#' 2. Remove duplicate feature IDs and drop rows with NA values.
#' 3. Parse compare_group to identify experimental and control groups.
#' 4. Extract sample lists for each group from samplelist.
#' 5. If node_list is provided, filter features to specified subset.
#' 6. Create separate abundance tables for case and control groups containing
#'    only samples from their respective groups.
#'
#' @importFrom dplyr distinct filter select
#' @importFrom tidyr drop_na
#' @importFrom tibble as_tibble_col
#' @importFrom purrr map_int
#' @importFrom stringr str_split
#'
#' @keywords internal
#'
#' @noRd
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
