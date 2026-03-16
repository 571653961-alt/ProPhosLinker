#' Perform functional enrichment analysis on subnetworks stratified by differential edge status.
#'
#' This function takes a `Differential_network` object (containing edges labeled by categories such as 
#' "Enhanced in case", "Only in control", etc.) and splits the network into subnetworks based on the 
#' `cor_status` of edges. It then performs functional enrichment analysis (e.g., GO, KEGG) separately 
#' for the nodes in each subnetwork using an external `run_enrichment` function. Finally, it reattaches 
#' the `cor_status` label to each enrichment result and aggregates them into a unified structure within 
#' the returned enrichment object.
#'
#' @param Differential_network An object of class `Differential_network` with slots `diff_edges` 
#'        (must contain column `cor_status`) and `diff_nodes`.
#' @param species Character string specifying the organism (e.g., "hsa" for human, "mmu" for mouse).
#' @param nCores Numeric: number of CPU cores to use for parallel enrichment computation.
#' @param omics_name Character: name of the omics layer (e.g., "transcriptome", "metabolome") for labeling.
#' @param enrichment_p_threshold Numeric: p-value cutoff for significant enrichment terms (default: 0.05).
#' @param compare_group Character: group comparison label used in enrichment metadata (default: `"data"`).
#' @param database_path Optional path to custom enrichment database files.
#' @param glist_path Optional path to gene list or background definition file.
#'
#' @return An enriched object (typically S4) returned by `run_enrichment`, with its `@network` slot 
#'         updated to include a `cor_status` column indicating which differential subnetwork each term belongs to.
#'         Returns `NULL` if no enrichment results are obtained.
#'
#' @importFrom dplyr filter
#' @import dplyr
#'
#' @examples
#' # Assuming `diff_net` is a valid Differential_network object
#' # enrich_result <- run_diff_enrichment(Differential_network = diff_net, species = "hsa")
#'
#' @export
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


  dfs_by_type <-lapply(names(network), function(type_name) {
    type_data <- network[[type_name]]
    modified_type <- lapply(type_data, function(df) {
      if (any(!is.na(df))) {
        df$cor_status <- type_name
      }
      df
    })
    modified_type
  }) |>
    setNames(names(network))
  

  all_df_names <- unique(unlist(lapply(dfs_by_type, names)))
  
  combined <- setNames(lapply(all_df_names, function(df_name) {
    do.call(rbind, lapply(dfs_by_type, `[[`, df_name))
  }), all_df_names)

  Enrichment@network<-combined
  return(Enrichment)
}