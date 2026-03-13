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
# 定义 show 方法
setMethod("show", "Stable_DifferentialNetwork", function(object) {
  cat("=== Stable Differential Network Results ===\n")
  cat("Use network_show(Network=differential_network1,plot_target='DifferentialNetwork',plot_type='control_stable_test',stable_num=4,R_threshold = 0.5) 
      to plot the network constructed from the first four bootstraps (only edges with an absolute value of R greater than 0.5 are kept)\n")
})
differential_network <- function(count_table = NULL,quantitative_table = NULL, samplelist = NULL, compare_group = NULL, annotation_table = NULL, diff_table= NULL,
                                 node_list = NULL,FC_threshold=1.2, p_threshold = 0.05,p_value_type="q_value",filter_num=1000, nBoots = 50, nCores = NULL,bootnet_R_threshold=0,stability_threshold = 0.2, 
                                 cor_method = "spearman", edge_FC_threshold=1.2,edge_p_threshold=0.05,enrichment_p_threshold=0.05,clustersize=25, 
                                 run_enrich = TRUE,species = NULL,omics_name=NULL,run_diffsubnet_enrich= TRUE,run_mediation = FALSE,
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
        diff_table<-run_process_intensity(diff_table)#集群蛋白的特殊处理
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
    if("omics_name" %in% colnames(annotation_table)){#多组学
      omics_num<-length(unique(annotation_table$omics_name))
      node_list<-annotation_table |>
        dplyr::filter(feature_ID %in% node_list) |>
        ## 按 omics_name 分组，并从每个组中选择前 100 个 feature_ID，如果组中少于 100 个则选择所有
        dplyr::group_by(omics_name) |>
        dplyr::mutate(row_num = dplyr::row_number()) |>
        dplyr::filter(row_num <= floor(filter_num/omics_num)) |>####100
        dplyr::select(-row_num) |>
        dplyr::pull(feature_ID)
      
    }else{#单组学
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
  #Functional Enrichment Analysis
  if (run_enrich && !is.null(Differential_network)) {
    tryCatch({
      Enrichment <- run_diff_enrichment(Differential_network=Differential_network,species = species,
                                        nCores=nCores,omics_name=omics_name,enrichment_p_threshold=enrichment_p_threshold,
                                        compare_group=compare_group,database_path=database_path,glist_path=glist_path)
    }, error = function(e) {
      message("Error in Functional Enrichment Analysis: ", e$message)
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
  # Functional Enrichment Analysis
  # if (run_enrich && !is.null(Differential_network)) {
  #   tryCatch({
  #     split_edge_data <- split(Differential_network@diff_edges, Differential_network@diff_edges$cor_status)
  #     subnetworks_nodes <- lapply(names(split_edge_data), function(x) {
  #       data<-split_edge_data[[x]]
  #       Differential_network@diff_nodes |>
  #         dplyr::filter(node %in% union(data$from,data$to))
  #     })
  #     names(subnetworks_nodes)<-names(split_edge_data)
  #     Enrichment <- run_enrichment(annotation_table_select = subnetworks_nodes, species = species,nCores=nCores,omics_name=omics_name,
  #                                  group_name=compare_group,enrichment_p_threshold=enrichment_p_threshold)
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
  #message("Enrichment Done")
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
  # message("Subnetwork Done")
  # Subnetwork Enrichment Analysis 
  if (run_diffsubnet_enrich && !is.null(Differential_subnetwork)) {
    tryCatch({
      DiffSubnet_Enrichment<-run_diffsubnet_enrichment(Differential_subnetwork=Differential_subnetwork,species = species,nCores=nCores,omics_name=omics_name,enrichment_p_threshold=enrichment_p_threshold,
                                                       compare_group=compare_group,database_path=database_path,glist_path=glist_path)
      
    }, error = function(e) {
      message("Error in Subnetwork Functional Enrichment Analysis: ", e$message)
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
  # Subnetwork Enrichment Analysis :only top1 subnet
  if (run_enrich && !is.null(Differential_subnetwork)) {
    tryCatch({
      TOPsubnet<-names(Differential_subnetwork@subnetworks)
      if(length(TOPsubnet)>0){
        TOPsubnet_edges <-Differential_subnetwork@subnetworks[[TOPsubnet]]$edges
        TOPsubnet_nodes <-Differential_subnetwork@subnetworks[[TOPsubnet]]$nodes
        split_edge_data <- split(TOPsubnet_edges, TOPsubnet_edges$cor_status)
        subnetworks_nodes <- lapply(names(split_edge_data), function(x) {
          data<-split_edge_data[[x]]
          TOPsubnet_nodes |>
            dplyr::filter(node %in% union(data$from,data$to))
        })
        names(subnetworks_nodes)<-names(split_edge_data)
        Enrichment_subnetwork <- run_enrichment(annotation_table_select = subnetworks_nodes, species = species,nCores=nCores,omics_name=omics_name,
                                                group_name=paste0(compare_group,"_",TOPsubnet),enrichment_p_threshold=enrichment_p_threshold)
      }else{
        Enrichment_subnetwork<-NULL
        message("No subnetwork result.")
      }
      
    }, error = function(e) {
      message("Error in Subnetwork Functional Enrichment Analysis: ", e$message)
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
                 Enrichment_subnetwork=Enrichment_subnetwork
      ))
    })
  }
  #######
  if (run_mediation && !is.null(Differential_subnetwork)) {
    tryCatch({
      TOPsubnet<-names(Differential_subnetwork@subnetworks)
      if(length(TOPsubnet)>0){
        TOPsubnet_nodes <-Differential_subnetwork@subnetworks[[TOPsubnet]]$nodes$node
        Mediation_subnetwork <- run_mediation(count_table = count_table, node_list=TOPsubnet_nodes,
                                              group_name=paste0(compare_group,"_",TOPsubnet),
                                              mediation_R2_threshold=mediation_R2_threshold,
                                              mediation_p_threshold=mediation_p_threshold,nCores=nCores)
      }else{
        Mediation_subnetwork<-NULL
        message("No subnetwork result.")
      }
      
    }, error = function(e) {
      message("Error in Subnetwork Mediation Analysis: ", e$message)
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
                 Enrichment_subnetwork=Enrichment_subnetwork,
                 Mediation_subnetwork=Mediation_subnetwork
      ))
    })
  }
  
  
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


