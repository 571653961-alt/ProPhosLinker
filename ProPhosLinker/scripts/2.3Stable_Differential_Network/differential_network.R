#' Stable_DifferentialNetwork S4 Class
#'
#' An S4 class that stores complete results from differential network analysis,
#' including data preprocessing, conditional networks, differential networks,
#' subnetworks, and functional enrichment results.
#'
#' @slot group_name character, comparison label (e.g., "Treatment-vs-Control")
#' @slot colormapping ANY, named color vector for node class annotations
#' @slot PreData ANY, PreData object containing preprocessed count table and sample metadata
#' @slot node_list ANY, character vector of feature IDs included in network analysis
#' @slot Diff_anno ANY, differential analysis results merged with feature annotations
#' @slot Conditional_network ANY, Conditional_network object with case and control networks
#' @slot Differential_network ANY, Differential_network object with edge status classifications
#' @slot Enrichment ANY, functional enrichment results for the differential network
#' @slot Differential_subnetwork ANY, SubNetwork object with clustered subnetworks
#' @slot DiffSubnet_Enrichment ANY, functional enrichment results for individual subnetworks
#'
#' @name Stable_DifferentialNetwork-class
#' @rdname Stable_DifferentialNetwork-class
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a Stable_DifferentialNetwork object
#' diff_network <- new("Stable_DifferentialNetwork",
#'   group_name = "Treatment-vs-Control",
#'   colormapping = color_map,
#'   PreData = preprocessed_data,
#'   node_list = selected_features,
#'   Diff_anno = differential_annotations,
#'   Conditional_network = cond_network,
#'   Differential_network = diff_net,
#'   Enrichment = enrichment_results,
#'   Differential_subnetwork = subnetwork_results,
#'   DiffSubnet_Enrichment = subnetwork_enrichment
#' )
#' 
#' # Access components
#' diff_network@group_name
#' diff_network@Differential_network
#' }
setClass("Stable_DifferentialNetwork", slots = c(
  group_name="character",
  colormapping = "ANY",
  PreData = "ANY",
  node_list = "ANY",
  Diff_anno = "ANY",
  Conditional_network="ANY",
  Differential_network="ANY",
  Enrichment="ANY",
  Differential_subnetwork="ANY",
  DiffSubnet_Enrichment="ANY"
))

setMethod("show", "Stable_DifferentialNetwork", function(object) {
  cat("=== Stable Differential Network Results ===\n")
  cat("Use network_show(Network=differential_network1,plot_target='DifferentialNetwork',plot_type='control_stable_test',stable_num=4,R_threshold = 0.5) 
      to plot the network constructed from the first four bootstraps (only edges with an absolute value of R greater than 0.5 are kept)\n")
})

#' Perform comprehensive differential network analysis
#'
#' This function orchestrates a complete differential network analysis pipeline,
#' integrating data preprocessing, differential expression analysis, conditional
#' network construction, differential network identification, subnetwork clustering,
#' and functional enrichment analysis.
#'
#' @param count_table data.frame, feature abundance table with 'feature_ID' column
#' @param quantitative_table data.frame, quantitative table for differential analysis.
#'                           If NULL, count_table is used. Default NULL.
#' @param samplelist data.frame, sample metadata with 'sample' and 'group' columns
#' @param compare_group character, comparison specification in format
#'                      "Experimental_group:Control_group"
#' @param annotation_table data.frame, feature annotations with at least 'feature_ID'
#'                         and 'Class' columns
#' @param diff_table data.frame, pre-computed differential analysis results.
#'                   If NULL, run_diff is called automatically. Default NULL.
#' @param node_list character vector, pre-specified feature list for network construction.
#'                  If NULL, features with differential status are used. Default NULL.
#' @param FC_threshold numeric, fold change threshold for differential analysis.
#'                     Default 1.2.
#' @param p_threshold numeric, p-value threshold for differential analysis. Default 0.05.
#' @param p_value_type character, which p-value to use: "p_value" or "q_value".
#'                     Default "q_value".
#' @param filter_num integer, maximum number of features to include in network.
#'                   If node_list exceeds this, features are filtered by omics type
#'                   or top-ranked. Default 1000.
#' @param nBoots integer, number of bootstrap resamples for stability analysis.
#'               Default 50.
#' @param nCores integer, number of CPU cores for parallel processing. If NULL,
#'               detects available cores automatically.
#' @param bootnet_R_threshold numeric, correlation threshold for bootnet construction.
#'                            Default 0.
#' @param stability_threshold numeric, edge stability threshold (0-1). Edges with
#'                             stability below this are filtered. Default 0.2.
#' @param cor_method character, correlation method: "pearson" or "spearman".
#'                   Default "spearman".
#' @param edge_FC_threshold numeric, fold change threshold for classifying differential
#'                          edges (case vs control correlation differences).
#'                          Default 1.2.
#' @param edge_p_threshold numeric, p-value threshold for differential edges.
#'                         Default 0.05.
#' @param enrichment_p_threshold numeric, p-value threshold for functional enrichment.
#'                                Default 0.05.
#' @param clustersize integer, maximum size for network clusters/subnetworks.
#'                    Default 25.
#' @param run_enrich logical, whether to perform functional enrichment analysis on
#'                   the differential network. Default TRUE.
#' @param species character, species name for enrichment analysis (e.g., "hsa" for human,
#'                "mmu" for mouse). Required if run_enrich is TRUE.
#' @param omics_name character, name of omics type for enrichment analysis
#'                   (e.g., "UNIPROT", "ENTREZID"). Default NULL.
#' @param run_diffsubnet_enrich logical, whether to perform enrichment analysis on
#'                              individual subnetworks. Default TRUE.
#' @param colormapping ANY, pre-defined color mapping for node classes.
#'                     If NULL, generated automatically. Default NULL.
#' @param PreData ANY, pre-existing PreData object. If NULL, created from inputs.
#'                Default NULL.
#' @param Diff_anno ANY, pre-existing differential annotations. If NULL, created.
#'                  Default NULL.
#' @param Conditional_network ANY, pre-existing Conditional_network object. Default NULL.
#' @param Differential_network ANY, pre-existing Differential_network object. Default NULL.
#' @param Enrichment ANY, pre-existing enrichment results. Default NULL.
#' @param Differential_subnetwork ANY, pre-existing SubNetwork object. Default NULL.
#' @param DiffSubnet_Enrichment ANY, pre-existing subnetwork enrichment results. Default NULL.
#' @param database_path character, path to enrichment database files. Default NULL.
#' @param glist_path character, path to gene list files for enrichment. Default NULL.
#'
#' @return A Stable_DifferentialNetwork S4 object containing all analysis results:
#'   \item{group_name}{Comparison label}
#'   \item{colormapping}{Color mapping for node classes}
#'   \item{PreData}{Preprocessed data object}
#'   \item{node_list}{List of features included in analysis}
#'   \item{Diff_anno}{Differential results with annotations}
#'   \item{Conditional_network}{Conditional network results (case and control)}
#'   \item{Differential_network}{Differential network with edge classifications}
#'   \item{Enrichment}{Functional enrichment results}
#'   \item{Differential_subnetwork}{Subnetwork clustering results}
#'   \item{DiffSubnet_Enrichment}{Subnetwork enrichment results (if run)}
#'
#' @details
#' Complete analysis pipeline:
#'
#' **1. Data Preprocessing (run_predata)**
#'    - Validates count_table and samplelist
#'    - Filters samples by group
#'    - Removes duplicates
#'    - Ensures proper data types
#'
#' **2. Differential Analysis (run_diff)**
#'    - Performs t-tests between groups
#'    - Calculates fold changes and p-values
#'    - Classifies features as up/down-regulated
#'
#' **3. Node Selection**
#'    - Identifies differentially expressed features
#'    - Filters to top N features if exceeding filter_num
#'    - Optionally balances by omics type
#'
#' **4. Conditional Network Construction (run_conditional_network)**
#'    - Builds separate correlation networks for case and control
#'    - Performs bootstrap stability analysis
#'    - Filters unstable edges
#'
#' **5. Differential Network (run_diff_network)**
#'    - Compares case and control edge correlations
#'    - Classifies edges as:
#'      * Only in case/control
#'      * Enhanced in case/control
#'      * Conflict relation
#'      * Non-significant
#'
#' **6. Subnetwork Clustering (run_subnet_cluster, run_add_cluster)**
#'    - Applies community detection algorithms
#'    - Ensures clusters ≤ clustersize
#'    - Extracts individual subnetworks
#'
#' **8. Subnetwork Enrichment (Optional)**
#'    - Performs enrichment on individual subnetworks
#'
#' @importFrom dplyr filter pull left_join
#' @importFrom methods new
#'
#' @examples
#' \dontrun{
#' # Complete differential network analysis
#' diff_network <- differential_network(
#'   count_table = abundance_data,
#'   samplelist = sample_metadata,
#'   compare_group = "Treatment:Control",
#'   annotation_table = feature_annotations,
#'   FC_threshold = 1.5,
#'   p_threshold = 0.05,
#'   nBoots = 100,
#'   stability_threshold = 0.3,
#'   species = "hsa",
#'   omics_name = "UNIPROT"
#' )
#' 
#' # Access results
#' diff_network@Differential_network       # Differential network
#' diff_network@Differential_subnetwork    # Subnetwork clusters
#' diff_network@Enrichment                 # Enrichment results
#' 
#' # Visualize using network_show
#' net_plot <- network_show(
#'   Network = diff_network,
#'   plot_type = "differential_network",
#'   show_node_name = TRUE,
#'   node_colortype = "Log2FC"
#' )
#' }
#'
#' @export

differential_network <- function(count_table = NULL,quantitative_table = NULL, samplelist = NULL, compare_group = NULL, annotation_table = NULL, diff_table= NULL,
                              node_list = NULL,FC_threshold=1.2, p_threshold = 0.05,p_value_type="q_value",filter_num=1000, nBoots = 50, nCores = NULL,bootnet_R_threshold=0,stability_threshold = 0.2, 
                              cor_method = "spearman", edge_FC_threshold=1.2,edge_p_threshold=0.05,enrichment_p_threshold=0.05,clustersize=25, 
                              run_enrich = TRUE,species = NULL,omics_name=NULL,run_diffsubnet_enrich= TRUE,
                              colormapping = NULL, PreData=NULL,Diff_anno = NULL,Conditional_network = NULL, Differential_network = NULL, Enrichment = NULL,
                              Differential_subnetwork=NULL,DiffSubnet_Enrichment=NULL,database_path=NULL,glist_path=NULL) {

  comparison<-gsub(":","-vs-", compare_group)
  groupname <-unlist(strsplit(x = compare_group, split = ":"))
  # Color Mapping
  if (is.null(colormapping)) {
    tryCatch({
      colormapping <- run_color(annotation_table)
    }, error = function(e) {
      message("Error in Color Mapping: ", e$message)
    })
  }
  ####Data checking
  if (is.null(PreData)) {
    tryCatch({
      PreData <- run_predata(count_table = count_table, samplelist = samplelist, group_name = groupname)
    }, error = function(e) {
      message("Error in Data checking: ", e$message)
      return(new("Stable_DifferentialNetwork",
                 group_name=comparison,
                 colormapping = colormapping,
                 PreData = PreData,
                 node_list = node_list,
                 Diff_anno = Diff_anno,
                 Conditional_network=Conditional_network,
                 Differential_network=Differential_network,
                 Enrichment=Enrichment,
                 Differential_subnetwork=Differential_subnetwork,
                 DiffSubnet_Enrichment=DiffSubnet_Enrichment
      ))
    })
  }
  count_table<-PreData@count_table
  samplelist<-PreData@samplelist

  # Differential analysis
  if (is.null(diff_table)) {
    if(is.null(quantitative_table)){
      message("Without quantitative_table, count_table will be used for differential analysis, please make sure that count_table does not have negative values!")
      quantitative_table<-count_table
    }
    tryCatch({
      diff_table <- run_diff(quantitative_table=quantitative_table,samplelist=samplelist,compare_group=compare_group,
               p_threshold=p_threshold,FC_threshold=FC_threshold,p_value_type=p_value_type)
    }, error = function(e) {
      message("Error in Differential analysis: ", e$message)
      return(new("Stable_DifferentialNetwork",
                 group_name=comparison,
                 colormapping = colormapping,
                 PreData = PreData,
                 node_list = node_list,
                 Diff_anno = Diff_anno,
                 Conditional_network=Conditional_network,
                 Differential_network=Differential_network,
                 Enrichment=Enrichment,
                 Differential_subnetwork=Differential_subnetwork,
                 DiffSubnet_Enrichment=DiffSubnet_Enrichment
      ))
    })
  }
  #mean_table
  if (!is.null(diff_table)) {
  tryCatch({
    if("Intensity(2)" %in% colnames(diff_table) && "Intensity(1)" %in% colnames(diff_table)){
      diff_table<-run_process_intensity(diff_table)
    }else{
      mean_table<-run_mean(count_table=count_table,samplelist=samplelist,compare_group=compare_group)
      diff_table<-diff_table |>
        dplyr::left_join(mean_table,by=c("feature_ID"="feature_ID"))}

  }, error = function(e){
    message("Error in mean_table: ", e$message)
    return(new("Stable_DifferentialNetwork",
               group_name=comparison,
               colormapping = colormapping,
               PreData = PreData,
               node_list = node_list,
               Diff_anno = Diff_anno,
               Conditional_network=Conditional_network,
               Differential_network=Differential_network,
               Enrichment=Enrichment,
               Differential_subnetwork=Differential_subnetwork,
               DiffSubnet_Enrichment=DiffSubnet_Enrichment
    ))
  })
  }
  # node_list
  if (is.null(node_list) && !is.null(diff_table)) {
    tryCatch({
    node_list<-diff_table |>
      dplyr::filter(!(State %in% c("Non-significant","Non"))) |>
      dplyr::pull(feature_ID)
    }, error = function(e) {
      message("Error in node_list: ", e$message)
      return(new("Stable_DifferentialNetwork",
                 group_name=comparison,
                 colormapping = colormapping,
                 PreData = PreData,
                 node_list = node_list,
                 Diff_anno = Diff_anno,
                 Conditional_network=Conditional_network,
                 Differential_network=Differential_network,
                 Enrichment=Enrichment,
                 Differential_subnetwork=Differential_subnetwork,
                 DiffSubnet_Enrichment=DiffSubnet_Enrichment
      ))
    })
  }
  node_list<-unique(node_list)
  if(length(node_list)==0){
    stop("There are no nodes for building the network. 
         Please enter node_list or adjust the difference analysis threshold (FC_threshold, p_threshold) to obtain sufficient difference features!")
  }
  if(length(node_list)>filter_num){
    message(paste0("The number of node_list is ", length(node_list), ",only the top ",filter_num, " features are analysed!"))
    if("omics_name" %in% colnames(annotation_table)){
      omics_num<-length(unique(annotation_table$omics_name))
      node_list<-annotation_table |>
        dplyr::filter(feature_ID %in% node_list) |>
        dplyr::group_by(omics_name) |>
        dplyr::mutate(row_num = dplyr::row_number()) |>
        dplyr::filter(row_num <= floor(filter_num/omics_num)) |>
        dplyr::select(-row_num) |>
        dplyr::pull(feature_ID)
      
    }else{
      node_list<-node_list[1:filter_num]
    }
  }
  # Integration of Differential table and annotation table
  if (is.null(Diff_anno) && !is.null(diff_table)) {
    tryCatch({
      Diff_anno <- diff_table |>
       # dplyr::select(feature_ID,FC,p_value,State) |>
        dplyr::left_join(annotation_table, by=c("feature_ID"="feature_ID"))
        
    }, error = function(e) {
      message("Error in Integration analysis: ", e$message)
      return(new("Stable_DifferentialNetwork",
                 group_name=comparison,
                 colormapping = colormapping,
                 PreData = PreData,
                 node_list = node_list,
                 Diff_anno = Diff_anno,
                 Conditional_network=Conditional_network,
                 Differential_network=Differential_network,
                 Enrichment=Enrichment,
                 Differential_subnetwork=Differential_subnetwork,
                 DiffSubnet_Enrichment=DiffSubnet_Enrichment
      ))
    })
  }
  # Conditional network
  if (is.null(Conditional_network) && !is.null(Diff_anno)  && !is.null(node_list)) {
    tryCatch({
      Conditional_network<-run_conditional_network(count_table = count_table, samplelist = samplelist,compare_group=compare_group,
                                                   Diff_anno= Diff_anno, node_list = node_list,cor_method = cor_method, nBoots = nBoots,
                                                   nCores = nCores, stability_threshold = stability_threshold,bootnet_R_threshold=bootnet_R_threshold
                                                   )
      
    }, error = function(e) {
      message("Error in Conditional network analysis: ", e$message)
      return(new("Stable_DifferentialNetwork",
                 group_name=comparison,
                 colormapping = colormapping,
                 PreData = PreData,
                 node_list = node_list,
                 Diff_anno = Diff_anno,
                 Conditional_network=Conditional_network,
                 Differential_network=Differential_network,
                 Enrichment=Enrichment,
                 Differential_subnetwork=Differential_subnetwork,
                 DiffSubnet_Enrichment=DiffSubnet_Enrichment
      ))
    })
  }
  
  # Differential network
  if (is.null(Differential_network) && !is.null(Conditional_network)) {
    tryCatch({
      Differential_network<-run_diff_network(Conditional_network=Conditional_network,
                                             edge_FC_threshold=edge_FC_threshold,edge_p_threshold=edge_p_threshold,
                                             compare_group=compare_group)

    }, error = function(e) {
      message("Error in Differential network analysis: ", e$message)
      return(new("Stable_DifferentialNetwork",
                 group_name=comparison,
                 colormapping = colormapping,
                 PreData = PreData,
                 node_list = node_list,
                 Diff_anno = Diff_anno,
                 Conditional_network=Conditional_network,
                 Differential_network=Differential_network,
                 Enrichment=Enrichment,
                 Differential_subnetwork=Differential_subnetwork,
                 DiffSubnet_Enrichment=DiffSubnet_Enrichment
      ))
    })
  }
  # #Functional Enrichment Analysis
  # if (run_enrich && !is.null(Differential_network)) {
  #   tryCatch({
  #     Enrichment <- run_diff_enrichment(Differential_network=Differential_network,species = species,
  #                     nCores=nCores,omics_name=omics_name,enrichment_p_threshold=enrichment_p_threshold,
  #                     compare_group=compare_group,database_path=database_path,glist_path=glist_path)
  #   }, error = function(e) {
  #     message("Error in Functional Enrichment Analysis: ", e$message)
  #     return(new("Stable_DifferentialNetwork",
  #                group_name=comparison,
  #                colormapping = colormapping,
  #                PreData = PreData,
  #                node_list = node_list,
  #                Diff_anno = Diff_anno,
  #                Conditional_network=Conditional_network,
  #                Differential_network=Differential_network,
  #                Enrichment=Enrichment,
  #                Differential_subnetwork=Differential_subnetwork,
  #                DiffSubnet_Enrichment=DiffSubnet_Enrichment
  #     ))
  #   })
  # }
  
  # Subnetwork Identification 
  if (is.null(Differential_subnetwork) && !is.null(Differential_network)) {
    tryCatch({
      nodes <- Differential_network@diff_nodes
      edges <- Differential_network@diff_edges
     # Conditional_network_layout<- Differential_network@Conditional_network_layout
      clusters <- run_subnet_cluster(nodes=nodes, edges=edges,clustersize=clustersize)
      Differential_subnetwork<-run_add_cluster(cfg_t1 = clusters, nodes = nodes, 
                                               edges = edges, group_name = comparison,
                                               diffmessage="diff")#,Conditional_network_layout=Conditional_network_layout)#,
                                              # module_select="top1")#, bootnet_list = bootnet_list)
      
    }, error = function(e) {
      message("Error in Subnetwork Identification analysis: ", e$message)
      return(new("Stable_DifferentialNetwork",
                 group_name=comparison,
                 colormapping = colormapping,
                 PreData = PreData,
                 node_list = node_list,
                 Diff_anno = Diff_anno,
                 Conditional_network=Conditional_network,
                 Differential_network=Differential_network,
                 Enrichment=Enrichment,
                 Differential_subnetwork=Differential_subnetwork,
                 DiffSubnet_Enrichment=DiffSubnet_Enrichment
      ))
    })
  }
  # # message("Subnetwork Done")
  # # Subnetwork Enrichment Analysis 
  # if (run_diffsubnet_enrich && !is.null(Differential_subnetwork)) {
  #   tryCatch({
  #     DiffSubnet_Enrichment<-run_diffsubnet_enrichment(Differential_subnetwork=Differential_subnetwork,species = species,nCores=nCores,omics_name=omics_name,enrichment_p_threshold=enrichment_p_threshold,
  #                                         compare_group=compare_group,database_path=database_path,glist_path=glist_path)
  #     
  #   }, error = function(e) {
  #     message("Error in Subnetwork Functional Enrichment Analysis: ", e$message)
  #     return(new("Stable_DifferentialNetwork",
  #                group_name=comparison,
  #                colormapping = colormapping,
  #                PreData = PreData,
  #                node_list = node_list,
  #                Diff_anno = Diff_anno,
  #                Conditional_network=Conditional_network,
  #                Differential_network=Differential_network,
  #                Enrichment=Enrichment,
  #                Differential_subnetwork=Differential_subnetwork,
  #                DiffSubnet_Enrichment=DiffSubnet_Enrichment
  #     ))
  #   })
  # }
  # 

  
  # Prepare output
  Differential_network <- new("Stable_DifferentialNetwork",
                           group_name=comparison,
                           colormapping = colormapping,
                           PreData = PreData,
                           node_list = node_list,
                           Diff_anno = Diff_anno,
                           Conditional_network=Conditional_network,
                           Differential_network=Differential_network,
                           Enrichment=Enrichment,
                           Differential_subnetwork=Differential_subnetwork,
                           DiffSubnet_Enrichment=DiffSubnet_Enrichment
  )
  
  return(Differential_network)
}
