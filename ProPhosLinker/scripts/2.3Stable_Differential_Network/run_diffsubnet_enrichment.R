#' Differential Subnetwork Enrichment Analysis
#'
#' This function performs enrichment analysis specifically for differential subnetworks
#' identified from network analysis. It processes each subnetwork separately, splits
#' edges based on correlation status, and runs enrichment analysis on the nodes
#' associated with each edge type.
#'
#' @param Differential_subnetwork A DiffSubnetwork object containing subnetworks with
#'        nodes and edges, where edges have cor_status indicating their differential
#'        behavior (e.g., "Enhanced in case", "Enhanced in control", "Conflict relation")
#' @param species Character string specifying the species for KEGG enrichment
#'        (e.g., "hsa" for Homo sapiens)
#' @param nCores Integer specifying number of CPU cores for parallel processing.
#'        If NULL, uses detectCores() - 1
#' @param omics_name Character vector specifying omics types to include in enrichment
#' @param enrichment_p_threshold Numeric p-value threshold for significance filtering
#'        (default: 0.05)
#' @param compare_group Character string specifying comparison groups in format
#'        "case:control" (e.g., "T:N")
#' @param database_path Character string path to local KEGG database files
#' @param glist_path Character string path to gene list files for species
#'
#' @return A DiffSubnet_Enrichment object containing enrichment results for each
#'         subnetwork, with results split by correlation status. Returns NULL if
#'         no subnetworks exist or enrichment fails.
#'
#' @details
#' The function processes each subnetwork by:
#'   1. Splitting edges based on cor_status (e.g., different differential patterns)
#'   2. Extracting nodes associated with each edge type
#'   3. Running enrichment analysis on each node set using run_enrichment()
#'   4. Combining results across subnetworks with proper naming
#'
#' The resulting object contains nested lists where:
#'   - Top level: Subnetwork names
#'   - Second level: Correlation status types
#'   - Third level: Standard enrichment result data frames
#'
#' @seealso
#'   \code{\link{run_enrichment}} for the underlying enrichment analysis
#'   \code{\link{DiffSubnet_Enrichment-class}} for the result class structure
#'
#' @examples
#' \dontrun{
#' # Run enrichment on differential subnetworks
#' enrichment_results <- run_diffsubnet_enrichment(
#'   Differential_subnetwork = diff_subnet_object,
#'   species = "hsa",
#'   nCores = 4,
#'   omics_name = c("Proteomics", "Phosphoproteomics"),
#'   enrichment_p_threshold = 0.05,
#'   compare_group = "T:N"
#' )
#' 
#' # Access results for a specific subnetwork and correlation status
#' subnet1_enhanced_case <- enrichment_results@network$subnet_1$`Enhanced in case`
#' }
#'
#' @export


setClass("DiffSubnet_Enrichment", slots = c(
  group_name="character",
  network = "ANY"
))
run_diffsubnet_enrichment<-function(Differential_subnetwork=NULL,species = NULL,nCores=NULL,omics_name=NULL,enrichment_p_threshold=0.05,
                                    compare_group="data",database_path=NULL,glist_path=NULL){
  
  if(length(compare_group)>1){
    group_name<-paste(compare_group, collapse = "-vs-")
  }else{
    group_name<-gsub(":","-vs-", compare_group)#diff
  }
  subnet_name<-names(Differential_subnetwork@subnetworks)
  if(length(subnet_name)>0){
    Enrichment_results <- list()  
    for (i in 1:length(subnet_name)) {
 #    lapply(c(1:length(subnet_name)), function(i){
      subnetname<-subnet_name[i]
      subnet_edges <-Differential_subnetwork@subnetworks[[subnetname]]$edges
      subnet_nodes <-Differential_subnetwork@subnetworks[[subnetname]]$nodes
      split_edge_data <- split(subnet_edges, subnet_edges$cor_status)
      subnetworks_nodes <- lapply(names(split_edge_data), function(x) {
        data<-split_edge_data[[x]]
        subnet_nodes |>
          dplyr::filter(node %in% union(data$from,data$to))
      })
      names(subnetworks_nodes)<-names(split_edge_data)
      Enrichment_results[[i]] <- run_enrichment(annotation_table_select = subnetworks_nodes, species = species,nCores=nCores,omics_name=omics_name,
                                                group_name=subnetname,enrichment_p_threshold=enrichment_p_threshold,database_path=database_path,
                                                glist_path=glist_path)
   }
    #  })
   if (all(vapply(Enrichment_results, function(x) is.null(x) || is.null(x@network), TRUE))) {
      return(NULL)
    } 
    names(Enrichment_results)=subnet_name
   # Enrichment_subnetwork_list <- purrr::map(Enrichment_results, ~ .x@network)
    
     Enrichment_subnetwork_list <-purrr::map( Enrichment_results,function(x){
       if(is.null(x)){
         return(NULL)
       }else{
         return(x@network)
       }
     }) 
    #   purrr::map(~ .x@network) |>
    #   purrr::imap(~ purrr::set_names(.x, paste0(.y, "_", names(.x))) ) |>
    #   purrr::flatten()

      result <- lapply(names(Enrichment_subnetwork_list), function(net) {
        subnet=Enrichment_subnetwork_list[[net]]
        if(all(is.null(subnet))){
          return(NULL)
        }
       
        dfs_by_type <-lapply(names(subnet), function(type_name) {
          type_data <- subnet[[type_name]]
          modified_type <- lapply(type_data, function(df) {
            if (any(!is.na(df))) {
              df$cor_status <- type_name
            }
            df
          })
          modified_type  
        }) |>
          setNames(names(subnet)) 
        
       
        all_df_names <- unique(unlist(lapply(dfs_by_type, names)))
        
        combined <- setNames(lapply(all_df_names, function(df_name) {
          do.call(rbind, lapply(dfs_by_type, `[[`, df_name))
        }), all_df_names)
        combined
      }) |>
        setNames(names(Enrichment_subnetwork_list))  

    Enrichment_subnetwork_result <- new("DiffSubnet_Enrichment",
                                      group_name=group_name,
                                      network =result
    )
  }else{
    Enrichment_subnetwork_result<-NULL
    message("No subnetwork result.")
  }
  return(Enrichment_subnetwork_result)
}

