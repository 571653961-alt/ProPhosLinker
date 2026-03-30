#' Save and visualize network analysis pipeline results
#'
#' This comprehensive pipeline function saves network analysis results including
#' stability tests, overall networks, clustering outputs, and subnetworks.
#' It generates visualizations, exports data files, and performs functional
#' enrichment analysis for identified modules.
#'
#' @param Network S4 object, network result object (Conditional_network or
#'                other network class containing analysis results)
#' @param stable_num integer, number of stable edges to retain in stability
#'                   visualization. Default 4.
#' @param richfactor_threshold numeric, threshold for enrichment factor filtering.
#'                             Default 0.
#' @param plot_title_size numeric, font size for plot titles. Default 12.
#' @param axis_title_size numeric, font size for axis titles. Default 8.
#' @param text_size numeric, general text size for plots. Default 8.
#' @param legend_title_size numeric, font size for legend titles. Default 8.
#' @param legend_text_size numeric, font size for legend text. Default 8.
#' @param font_family character, font family for text in plots. Default "sans".
#' @param image_margin_size numeric, margin size around images. Default 0.3.
#' @param node_name_size numeric, size of node labels. Default 1.5.
#' @param R_threshold numeric, correlation threshold for edge filtering.
#'                    Default 0.3.
#' @param ModuleSize_show integer, minimum module size to display. Default 0.
#' @param top_module_num integer, number of top modules to show. Default 20.
#' @param omics1_name character, name of first omics type (e.g., 'Pro' for
#'                    proteomics). Default 'Pro'.
#' @param omics2_name character, name of second omics type (e.g., 'Phos' for
#'                    phosphoproteomics). Default 'Phos'.
#' @param phos_pro data.frame, mapping between phosphosite and protein identifiers.
#' @param max_subnet_num integer, maximum number of subnetworks to process.
#'                       Default 8.
#' @param enrich_fromType character, database identifier type for enrichment
#'                        analysis (e.g., 'UNIPROT'). Default 'UNIPROT'.
#' @param edge_color_pos character, color for positive correlation edges.
#'                       Default "#9b6a65".
#' @param edge_color_neg character, color for negative correlation edges.
#'                       Default "#5d8992".
#' @param Enhanced_in_N character, color for edges enhanced in control group.
#'                      Default "#5d8992".
#' @param Enhanced_in_T character, color for edges enhanced in treatment group.
#'                      Default "#9b6a65".
#' @param Only_in_N character, color for edges present only in control group.
#'                  Default "#0c2b32".
#' @param Only_in_T character, color for edges present only in treatment group.
#'                  Default "#381512".
#' @param Conflict_relation character, color for conflicting edge relationships.
#'                          Default '#808080'.
#' @param color_gradient_low character, low color for continuous gradients.
#'                           Default "#175663".
#' @param color_gradient_high character, high color for continuous gradients.
#'                            Default "#90362d".
#' @param fill_gradientn_color character vector, colors for gradient fill
#'                              in subnetwork plots. Default c("#175663", "#dce6e9", "#90362d").
#' @param outdir character, output directory path for saving results. Default "./".
#'
#' @return No return value. Results are saved to files in the output directory:
#'   \item{2.3.0StabilityTest/}{Stability test visualizations and node/edge tables}
#'   \item{2.3.1Overall_Network/}{Overall differential network visualizations,
#'         node/edge tables, and functional enrichment results}
#'   \item{2.3.2Network_Clustering/}{Clustered network visualizations and
#'         community detection results}
#'   \item{2.3.3Sub_Network/}{Individual subnetwork visualizations,
#'         node/edge tables, and functional enrichment for each subnetwork}
#'
#' @details
#' Pipeline execution steps:
#'
#' 1. **Stability Testing** (2.3.0StabilityTest):
#'    - Generate stability plots for both case and control groups
#'    - Export node and edge tables for stable networks
#'    - Save stability test visualizations as PNG files
#'
#' 2. **Overall Network** (2.3.1Overall_Network):
#'    - Create differential network visualization
#'    - Export node and edge data with layout coordinates
#'    - Perform functional enrichment on omics features:
#'      * GO (Gene Ontology) enrichment
#'      * KEGG pathway enrichment
#'    - Save enrichment results and plots
#'
#' 3. **Network Clustering** (2.3.2Network_Clustering):
#'    - Generate clustered network visualization
#'    - Export cluster node/edge tables
#'    - Perform community detection analysis
#'
#' 4. **Subnetwork Analysis** (2.3.3Sub_Network):
#'    - Process top N subnetworks (controlled by max_subnet_num)
#'    - For each subnetwork:
#'      * Create differential subnetwork visualization
#'      * Export node and edge tables with layouts
#'      * Generate subnetwork-specific plots with edge color coding
#'      * Perform functional enrichment analysis
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggsave
#' @importFrom readr write_delim
#' @importFrom tryCatch withCallingHandlers
#'
#' @examples
#' \dontrun{
#' # Run full pipeline with custom parameters
#' pipline_save(
#'   Network = conditional_network,
#'   stable_num = 5,
#'   plot_title_size = 14,
#'   node_name_size = 2,
#'   max_subnet_num = 10,
#'   omics1_name = "Protein",
#'   omics2_name = "Phosphosite",
#'   phos_pro = phosphosite_mapping,
#'   outdir = "./network_results"
#' )
#' 
#' # Minimal execution with defaults
#' pipline_save(
#'   Network = network_object,
#'   phos_pro = mapping_table,
#'   outdir = "./results"
#' )
#' }
#'
#' @export

pipline_save<-function(Network=NULL,
                       stable_num=4,
                       richfactor_threshold=0,
                       plot_title_size=12,
                       axis_title_size=8,
                       text_size=8,
                       legend_title_size=8,
                       legend_text_size=8,
                       font_family="sans",
                       image_margin_size=0.3,
                       node_name_size=1.5,
                       R_threshold=0.3,
                       ModuleSize_show = 0,
                       top_module_num =20,
                       omics1_name = 'Pro',
                       omics2_name = 'Phos',
                       phos_pro = phos_pro,
                       max_subnet_num = 8,
                       enrich_fromType = 'UNIPROT',
                       edge_color_pos = "#9b6a65",
                       edge_color_neg = "#5d8992", 
                       Enhanced_in_N = "#5d8992", 
                       Enhanced_in_T = "#9b6a65",
                       Only_in_N = "#0c2b32",
                       Only_in_T = "#381512",
                       Conflict_relation = '#808080',
                       color_gradient_low = "#175663",
                       color_gradient_high = "#90362d",
                       fill_gradientn_color = c("#175663", "#dce6e9", "#90362d"),
                       outdir="./"){ 

  group_name=Network@group_name                   #T-vs-N
  split_result <- strsplit(group_name, "-vs-")[[1]]
  Networkclass=class(Network)[1]
  if(length(split_result)>1){
    casename <- split_result[1] 
    controlname <- split_result[2] 
  }else{
    casename <- "case" 
    controlname <- "control"
  }
  

  differential_network1=Network
  # 0. Stable_DifferentialNetwork\StabilityTest\T
 
  outdir1<-file.path(outdir,"2.3.0StabilityTest")
  tryCatch({
  bootnetplot1<-network_show(Network=differential_network1,plot_type="case_stable_test",image_margin_size=image_margin_size,
                             plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                             stable_num=stable_num,R_threshold = R_threshold)
  outdir1<-file.path(outdir,"2.3.0StabilityTest",casename)
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  boot_names <- names(bootnetplot1@data)
  lapply(boot_names,function(net_name){
    boot_data <- bootnetplot1@data[[net_name]]
    edges<-as.data.frame(boot_data$edges)
    nodes<-as.data.frame(boot_data$nodes) |>
      dplyr::mutate(x=boot_data$plot_layout[,1],y=boot_data$plot_layout[,2])
    net_name <- gsub(" ", "", net_name)
    readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",net_name,"_",casename,".tsv")),delim="\t")
    readr::write_delim(edges,file.path(outdir1,paste0("edges_",net_name,"_",casename,".tsv")),delim="\t")
  })
  save_path <- file.path(outdir, "2.3.0StabilityTest", paste0("Bootnet_", casename, ".png"))
  ggplot2::ggsave(
    filename = save_path, 
    plot = bootnetplot1@plot,
    width = 10,    
    height = 8,    
    dpi = 300,     
    bg = "white"   
  )

  }, error = function(e) {
    message("Error in case_stable_test: ", e$message)
  })
  
  # 0. Stable_DifferentialNetwork\StabilityTest\T
  tryCatch({
  bootnetplot2<-network_show(Network=differential_network1,plot_type="control_stable_test",image_margin_size=image_margin_size,
                             plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                             stable_num=stable_num,R_threshold = R_threshold)
  outdir1<-file.path(outdir,"2.3.0StabilityTest",controlname)
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  boot_names <- names(bootnetplot2@data)
  lapply(boot_names,function(net_name){
    boot_data <- bootnetplot2@data[[net_name]]
    edges<-as.data.frame(boot_data$edges)
    nodes<-as.data.frame(boot_data$nodes) |>
      dplyr::mutate(x=boot_data$plot_layout[,1],y=boot_data$plot_layout[,2])
    net_name <- gsub(" ", "", net_name)
    readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",net_name,"_",controlname,".tsv")),delim="\t")
    readr::write_delim(edges,file.path(outdir1,paste0("edges_",net_name,"_",controlname,".tsv")),delim="\t")
  })

  save_path <- file.path(outdir, "2.3.0StabilityTest", paste0("Bootnet_", controlname, ".png"))
  ggplot2::ggsave(
    filename = save_path, 
    plot = bootnetplot1@plot,
    width = 10,    
    height = 8,    
    dpi = 300,     
    bg = "white"   
  )

  }, error = function(e) {
    message("Error in control_stable_test: ", e$message)
  })

  
  # 1. Stable_DifferentialNetwork\OverallNetwork
  tryCatch({
  if(any(!is.null(differential_network1@Differential_network))) {
    diffnetplot<-network_show(Network=differential_network1,plot_type="differential_network",
                              image_margin_size=image_margin_size,
                              plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                              node_colortype="Log2FC",focus=c("all"),node_size=3,node_name_size=node_name_size,
                              show_edge_legend = TRUE,show_node_legend = TRUE,show_node_name = TRUE
    )
    outdir1<-file.path(outdir,"2.3.1Overall_Network")
    if (!dir.exists(outdir1)) {
      dir.create(outdir1, recursive = TRUE)
    }
    for (net_name1 in names(diffnetplot@data)) {
        LAYOUT<-as.data.frame(diffnetplot@data[[net_name1]]$plot_layout)
        nodes<-as.data.frame(diffnetplot@data[[net_name1]]$nodes) |>
          dplyr::mutate(x=LAYOUT[,1],y=LAYOUT[,2]) 
        edges<-as.data.frame(diffnetplot@data[[net_name1]]$edges)
        readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",net_name1,".tsv")),delim="\t")
        readr::write_delim(edges,file.path(outdir1,paste0("edges_",net_name1,".tsv")),delim="\t") 
    }
    outdir1_enrich <- file.path(outdir1,"Functional_Enrichment")
    if (!dir.exists(outdir1_enrich)) {
      dir.create(outdir1_enrich, recursive = TRUE)
    }
    pro_list <- nodes[nodes$Class == omics1_name,]$node
    phos_list <- nodes[nodes$Class == omics2_name,]$node
    phos_list <- phos_pro[phos_pro[[omics2_name]]  %in% phos_list,][[omics1_name]]
    if(length(pro_list)==0){
      pro_list <- c()
    }
    if(length(phos_list)==0){
      phos_list <- c()
    }
   omics_enrichment_list(pro_list,phos_list,outdir1_enrich,
                          omics1_name=omics1_name,omics2_name=omics2_name,
                          enrich_fromType = enrich_fromType,
                          pvalueCutoff = 0.05,GO_showCategory=6,KEGG_showCategory=15,
                          color_gradient_low = color_gradient_low,
                          color_gradient_high = color_gradient_high
    )
    
  }
  }, error = function(e) {
    message("Error in differential_network: ", e$message)
  })



  #2. Stable_DifferentialNetwork\NetworkClustering
  
  tryCatch({
    if(any(!is.null(differential_network1@Differential_subnetwork)) 
       # &&
       # any(!is.null(differential_network1@Differential_subnetwork@overall_cluster_network))
    ){
  diff_subnetplot<-network_show(Network=differential_network1,plot_type="diff_overall_cluster_network",
                                focus=c("all"),
                                image_margin_size=image_margin_size,
                                plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                                show_edge_legend = TRUE,show_node_legend = TRUE,show_node_name = TRUE
                                )
  outdir1<-file.path(outdir,"2.3.2Network_Clustering")
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  nodes<-as.data.frame(diff_subnetplot@data$net$nodes) |>
    dplyr::mutate(x=diff_subnetplot@data$net$plot_layout[,1],y=diff_subnetplot@data$net$plot_layout[,2])
  edges<-as.data.frame(diff_subnetplot@data$net$edges)
  readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",group_name,".tsv")),delim="\t")
  readr::write_delim(edges,file.path(outdir1,paste0("edges_",group_name,".tsv")),delim="\t")
  diff_net_community_detection_plot(file.path(outdir1,paste0("nodes_",group_name,".tsv")),
                                    file.path(outdir1,paste0("edges_",group_name,".tsv")),
                                    outdir1,
                                    omics1_name,omics2_name,
                                    ModuleSize_show,top_module_num)

  }
  }, error = function(e) {
    message("Error in diff_overall_cluster_network: ", e$message)
  })
  
  #3. Stable_DifferentialNetwork\SubNetwork
if (any(!is.null(differential_network1@Differential_subnetwork))){

  tryCatch({
  netname<-names(differential_network1@Differential_subnetwork@subnetworks)
  if(length(netname)> max_subnet_num){
  netname_show <- netname[1:max_subnet_num]}
  else{
    netname_show <- netname
  }
  filelist<- lapply(netname_show, function(i){
    diff_subnet_top_plot<-network_show(Network=differential_network1,
                                       plot_type="differential_subnetwork",image_margin_size=image_margin_size,
                                       node_colortype="Log2FC",focus=c("all"),subnetwork_name = c(i),
                                       show_edge_legend = TRUE,show_node_legend = TRUE,plot_title_size=plot_title_size,
                                       font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                                       show_node_name = TRUE,node_size=5,node_name_size=node_name_size
    )
    outdir1<-file.path(outdir,"2.3.3Sub_Network",i) 
    if (!dir.exists(outdir1)) {
      dir.create(outdir1, recursive = TRUE)
    }
    for (net_name1 in names(diff_subnet_top_plot@data)) {
    nodes<-as.data.frame(diff_subnet_top_plot@data[[net_name1]]$nodes) |>
      dplyr::mutate(x=diff_subnet_top_plot@data[[net_name1]]$plot_layout[,1],y=diff_subnet_top_plot@data[[net_name1]]$plot_layout[,2]) 
    edges<-as.data.frame(diff_subnet_top_plot@data[[net_name1]]$edges)
    readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",i,"_",net_name1,".tsv")),delim="\t")
    readr::write_delim(edges,file.path(outdir1,paste0("edges_",i,"_",net_name1,".tsv")),delim="\t") 
    }
    Differential_subnetwork_plot(outdir1,i,group_name,omics1_name,omics2_name,
                                 edge_color_pos,
                                 edge_color_neg, 
                                 Enhanced_in_N, 
                                 Enhanced_in_T,
                                 Only_in_N,
                                 Only_in_T,
                                 Conflict_relation,
                                 fill_gradientn_color)
    outdir1_enrich <- file.path(outdir1,"Functional_Enrichment")
    if (!dir.exists(outdir1_enrich)) {
      dir.create(outdir1_enrich, recursive = TRUE)
    }
    pro_list <- nodes[nodes$Class == omics1_name,]$node
    phos_list <- nodes[nodes$Class == omics2_name,]$node
    phos_list <- phos_pro[phos_pro[[omics2_name]]  %in% phos_list,][[omics1_name]]
    if(length(pro_list)==0){
      pro_list <- c()
    }
    if(length(phos_list)==0){
      phos_list <- c()
    }
    omics_enrichment_list(pro_list,phos_list,outdir1_enrich,
                          omics1_name=omics1_name,omics2_name=omics2_name,
                          enrich_fromType = enrich_fromType,
                          pvalueCutoff = 0.05,GO_showCategory=6,KEGG_showCategory=15,
                          color_gradient_low = color_gradient_low,
                          color_gradient_high = color_gradient_high
                          )
  })
  }, error = function(e) {
    message("Error in differential_subnetwork: ", e$message)
  })
}
}