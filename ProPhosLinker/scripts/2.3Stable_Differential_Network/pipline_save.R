' Save network analysis results: This function saves the plots and data generated from the network analysis pipeline.
#'
#' @param Network A network object (e.g., Stable_SubNetwork, StableNetwork, etc.)
#' @param stable_num Number of bootstrap networks to display (for stable_test)
#' @param richfactor_threshold Threshold for enrichment rich factor (default: 0)
#' @param plot_title_size Size of plot title (default: 12)
#' @param axis_title_size Size of axis titles (default: 8)
#' @param text_size Size of text elements (default: 8)
#' @param legend_title_size Size of legend title (default: 8)
#' @param legend_text_size Size of legend text (default: 8)
#' @param font_family Font family for text (default: "Arial")
#' @param image_margin_size Margin size around plot (default: 0.3)
#' @param node_name_size Size of node labels (default: 2)
#' @param outdir output directory (default: ./)
#'
#' @return NULL. This function does not return any value, it only saves PDF images and data files.
#' @export
#'
#' @examples
#' data("stable_subnetwork_result")
#' pipline_save(stable_subnetwork_result)
#'
#' data("differential_network_result")
#' pipline_save(differential_network_result)
#'
#' data("multiplex_network_result")
#' pipline_save(multiplex_network_result)
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
  
if(Networkclass=="Stable_SubNetwork"){
  
  stable_subnetwork1=Network
  
  tryCatch({
  
  bootnetplot1<-network_show(Network=stable_subnetwork1,plot_type="stable_test",stable_num=stable_num,image_margin_size=image_margin_size,
                             plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                             R_threshold = R_threshold,node_size=3)
  outdir1<-file.path(outdir,Networkclass,"StabilityTest")
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
    readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",net_name,"_",group_name,".txt")),delim="\t")
    readr::write_delim(edges,file.path(outdir1,paste0("edges_",net_name,"_",group_name,".txt")),delim="\t") 
  })
  pdf(file = file.path(outdir1, paste0("Bootnet_",group_name,".pdf")), 
      width = bootnetplot1@picture_width, height = bootnetplot1@picture_height)
  print(bootnetplot1@plot)
  dev.off()
  
  }, error = function(e) {
    message("Error in stable_test: ", e$message)
  })
  
  tryCatch({

  stableplot1<-network_show(Network=stable_subnetwork1,plot_type="overall_network",R_threshold=R_threshold,image_margin_size=image_margin_size,
                            plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                            show_node_legend=TRUE)
  outdir1<-file.path(outdir,Networkclass,"OverallNetwork")
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  nodes<-as.data.frame(stableplot1@data$net$nodes) |>
    dplyr::mutate(x=stableplot1@data$net$plot_layout[,1],y=stableplot1@data$net$plot_layout[,2]) 
  edges<-as.data.frame(stableplot1@data$net$edges)
  readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",group_name,".txt")),delim="\t")
  readr::write_delim(edges,file.path(outdir1,paste0("edges_",group_name,".txt")),delim="\t")
  pdf(file = file.path(outdir1, paste0("Stable_Correlation_Network_",group_name,".pdf")), 
      width = stableplot1@picture_width, height = stableplot1@picture_height)
  print(stableplot1@plot)
  dev.off()
  }, error = function(e) {
    message("Error in overall_network: ", e$message)
  })
  
  if(any(!is.null(stable_subnetwork1@Enrichment)) &&
     any(!is.null(stable_subnetwork1@Enrichment@network)) &&
     any(!is.null(stable_subnetwork1@Enrichment@network$network))){
    ann_filter <- stable_subnetwork1@Enrichment@network$network$annotations_filter
  }else{
    ann_filter <-NULL
  }
  tryCatch({
    if (!is.null(ann_filter) && nrow(ann_filter) > 0) {
      bubbleDiagram_net<-network_show(Network=stable_subnetwork1,plot_type="enrichment",focus=c("all"),
                                      richfactor_threshold = richfactor_threshold,plot_title_size=plot_title_size,
                                      axis_title_size=axis_title_size,
                                      text_size=text_size,
                                      legend_title_size=legend_title_size,
                                      legend_text_size=legend_text_size
                                      
      )
      outdir1<-file.path(outdir,Networkclass,"OverallNetwork")  
      if (!dir.exists(outdir1)) {
        dir.create(outdir1, recursive = TRUE)
      }
      readr::write_delim(bubbleDiagram_net@data$net, file.path(outdir1, paste0("BubbleDiagram_",group_name,".txt")), delim = "\t")
      pdf(file = file.path(outdir1, paste0("BubbleDiagram_",group_name,".pdf")), 
          width = bubbleDiagram_net@picture_width, height = bubbleDiagram_net@picture_height)
      print(bubbleDiagram_net@plot)
      dev.off()
    }
  }, error = function(e) {
    message("Error in enrichment: ", e$message)
  })
  tryCatch({
    if(any(!is.null(stable_subnetwork1@SubNetwork)) &&
       any(!is.null(stable_subnetwork1@SubNetwork@overall_cluster_network)) 
    ){
  cluster_anno_plot1<-network_show(Network=stable_subnetwork1,plot_type="overall_cluster_network",image_margin_size=image_margin_size,
                                   plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                                   R_threshold = R_threshold,add_enrichement=TRUE,
                                   show_node_legend =TRUE)
  outdir1<-file.path(outdir,Networkclass,"NetworkClustering")
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  nodes<-as.data.frame(cluster_anno_plot1@data$net$nodes) |>
    dplyr::mutate(x=cluster_anno_plot1@data$net$plot_layout[,1],y=cluster_anno_plot1@data$net$plot_layout[,2]) 
  edges<-as.data.frame(cluster_anno_plot1@data$net$edges)
  add_anno<-as.data.frame(cluster_anno_plot1@other)
  readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",group_name,".txt")),delim="\t")
  readr::write_delim(edges,file.path(outdir1,paste0("edges_",group_name,".txt")),delim="\t")
  readr::write_delim(add_anno,file.path(outdir1,paste0("annotation_",group_name,".txt")),delim="\t")
  pdf(file = file.path(outdir1, paste0("Network_Clustering_",group_name,".pdf")), 
      width = cluster_anno_plot1@picture_width, height = cluster_anno_plot1@picture_height)
  print(cluster_anno_plot1@plot)
  dev.off()
    }
  }, error = function(e) {
    message("Error in overall_cluster_network: ", e$message)
  })
  tryCatch({
  if (any(!is.null(stable_subnetwork1@SubNetwork)) &&
        any(!is.null(stable_subnetwork1@SubNetwork@subnetworks))
        ){
  for (i in names(stable_subnetwork1@SubNetwork@subnetworks)) {
    subnet_plot<-network_show(Network=stable_subnetwork1,plot_type="sub_network",
                              subnetwork_name=c(i), node_colortype="Class",image_margin_size=image_margin_size,
                              plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                              R_threshold = R_threshold,show_node_name = TRUE,show_edge_legend = TRUE,show_node_legend=TRUE,
                              centrality_scatterplot=FALSE)
    outdir1<-file.path(outdir,Networkclass,"SubNetWork",i)
    outdir2<-file.path(outdir1,"Topological_analysis")
    if (!dir.exists(outdir2)) {
      dir.create(outdir2, recursive = TRUE)
    }
    nodes<-as.data.frame(subnet_plot@data[[1]]$nodes) |>
      dplyr::mutate(x=subnet_plot@data[[1]]$plot_layout[,1],y=subnet_plot@data[[1]]$plot_layout[,2]) 
    edges<-as.data.frame(subnet_plot@data[[1]]$edges)
    readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",i,"_",group_name,".txt")),delim="\t")
    readr::write_delim(edges,file.path(outdir1,paste0("edges_",i,"_",group_name,".txt")),delim="\t") 
    
    pdf(file = file.path(outdir1, paste0("SubNetWork_",i,"_",group_name,".pdf")), 
        width = subnet_plot@picture_width, height = subnet_plot@picture_height)
    print(subnet_plot@plot)
    dev.off()
    
    subnet_betweennessplot1<-network_show(Network=stable_subnetwork1,
                                          plot_type="sub_network",subnetwork_name=c(i), node_colortype="Class",image_margin_size=image_margin_size,
                                          plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,
                                          text_size=text_size,legend_text_size=legend_text_size,
                                          R_threshold = R_threshold,show_node_name = TRUE,show_edge_legend = TRUE,show_node_legend=TRUE,
                                          add_Centrality="betweenness",centrality_scatterplot=FALSE
                                          
    )
    nodes<-as.data.frame(subnet_betweennessplot1@data[[1]]$nodes) |>
      dplyr::mutate(x=subnet_betweennessplot1@data[[1]]$plot_layout[,1],y=subnet_betweennessplot1@data[[1]]$plot_layout[,2]) 
    edges<-as.data.frame(subnet_betweennessplot1@data[[1]]$edges)
    readr::write_delim(nodes,file.path(outdir2,paste0("Betweenness_",i,"_",group_name,"_nodes.txt")),delim="\t")
    readr::write_delim(edges,file.path(outdir2,paste0("Betweenness_",i,"_",group_name,"_edges.txt")),delim="\t") 
    
    pdf(file = file.path(outdir2, paste0("Betweenness_",i,"_",group_name,".pdf")), 
        width = subnet_betweennessplot1@picture_width, height = subnet_betweennessplot1@picture_height)
    print(subnet_betweennessplot1@plot)
    dev.off()
    
    subnet_degreeplot1<-network_show(Network=stable_subnetwork1,
                                     plot_type="sub_network",subnetwork_name=c(i), node_colortype="Class",image_margin_size=image_margin_size,
                                     plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,
                                     text_size=text_size,legend_text_size=legend_text_size,
                                     R_threshold = R_threshold,show_node_name = TRUE,show_edge_legend = TRUE,show_node_legend=TRUE,
                                     add_Centrality="degree",centrality_scatterplot=FALSE
                                     )
    nodes<-as.data.frame(subnet_degreeplot1@data[[1]]$nodes) |>
      dplyr::mutate(x=subnet_degreeplot1@data[[1]]$plot_layout[,1],y=subnet_degreeplot1@data[[1]]$plot_layout[,2]) 
    edges<-as.data.frame(subnet_degreeplot1@data[[1]]$edges)
    readr::write_delim(nodes,file.path(outdir2,paste0("Degree_",i,"_",group_name,"_nodes.txt")),delim="\t")
    readr::write_delim(edges,file.path(outdir2,paste0("Degree_",i,"_",group_name,"_edges.txt")),delim="\t") 
    
    pdf(file = file.path(outdir2, paste0("Degree_",i,"_",group_name,".pdf")), 
        width = subnet_degreeplot1@picture_width, height = subnet_degreeplot1@picture_height)
    print(subnet_degreeplot1@plot)
    dev.off()
    
    subnet_eigenvectorplot1<-network_show(Network=stable_subnetwork1,plot_type="sub_network",
                                          subnetwork_name=c(i), node_colortype="Class",image_margin_size=image_margin_size,
                                          plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,
                                          text_size=text_size,legend_text_size=legend_text_size,
                                          R_threshold = R_threshold,show_node_name = TRUE,show_edge_legend = TRUE,show_node_legend=TRUE,
                                          add_Centrality="eigenvector",centrality_scatterplot=FALSE
                                          )
    nodes<-as.data.frame(subnet_eigenvectorplot1@data[[1]]$nodes) |>
      dplyr::mutate(x=subnet_eigenvectorplot1@data[[1]]$plot_layout[,1],y=subnet_eigenvectorplot1@data[[1]]$plot_layout[,2]) 
    edges<-as.data.frame(subnet_eigenvectorplot1@data[[1]]$edges)
    readr::write_delim(nodes,file.path(outdir2,paste0("Eigenvector_",i,"_",group_name,"_nodes.txt")),delim="\t")
    readr::write_delim(edges,file.path(outdir2,paste0("Eigenvector_",i,"_",group_name,"_edges.txt")),delim="\t") 
    
    pdf(file = file.path(outdir2, paste0("Eigenvector_",i,"_",group_name,".pdf")), 
        width = subnet_eigenvectorplot1@picture_width, height = subnet_eigenvectorplot1@picture_height)
    print(subnet_eigenvectorplot1@plot)
    dev.off()
    
    subnet_centrality<- network_show(Network=stable_subnetwork1,plot_type="sub_network",subnetwork_name=c(i), node_colortype="Class",
                 R_threshold = R_threshold,show_node_name = TRUE,show_edge_legend = TRUE,image_margin_size=image_margin_size,
                 plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,
                 text_size=text_size,legend_text_size=legend_text_size,
                 show_node_legend=TRUE,add_Centrality=c("betweenness","degree","eigenvector")
                 )
    centrality<-as.data.frame(subnet_centrality@other@data[[1]]) 
    readr::write_delim(centrality,file.path(outdir2,paste0("Centrality_",i,"_",group_name,".txt")),delim="\t") 
    pdf(file = file.path(outdir2, paste0("Centrality_",i,"_",group_name,".pdf")), 
        width = subnet_centrality@picture_width, height = subnet_centrality@picture_height)
    print(subnet_centrality@other@plot[[1]])
    dev.off()
  }
  }
  }, error = function(e) {
    message("Error in sub_network: ", e$message)
  })
  tryCatch({
if (any(!is.null(stable_subnetwork1@Subnet_Enrichment))){
  netname<-names(stable_subnetwork1@Subnet_Enrichment@network)
  filelist<-lapply(netname, function(i){
    if(any(!is.null(stable_subnetwork1@Subnet_Enrichment)) &&
       any(!is.null(stable_subnetwork1@Subnet_Enrichment@network))){
      ann_filter <- stable_subnetwork1@Subnet_Enrichment@network[[i]]$annotations_filter
    }else{
      ann_filter<-NULL
    }
    if (!is.null(ann_filter) && nrow(ann_filter) > 0) {
      bubbleDiagram_subnet<-network_show(Network=stable_subnetwork1,subnetwork_name = c(i),
                                         plot_type="subnetwork_enrichment",
                                         richfactor_threshold = richfactor_threshold,plot_title_size=plot_title_size,
                                         axis_title_size=axis_title_size,
                                         text_size=text_size,
                                         legend_title_size=legend_title_size,
                                         legend_text_size=legend_text_size
                                         )
      outdir1<-file.path(outdir,Networkclass,"SubNetWork",i)
      if (!dir.exists(outdir1)) {
        dir.create(outdir1, recursive = TRUE)
      }
      readr::write_delim(bubbleDiagram_subnet@data[[1]], file.path(outdir1, paste0("BubbleDiagram_",i,"_",group_name,".txt")), delim = "\t")
      pdf(file = file.path(outdir1, paste0("BubbleDiagram_",i,"_",group_name,".pdf")), 
          width = bubbleDiagram_subnet@picture_width, height = bubbleDiagram_subnet@picture_height)
      print(bubbleDiagram_subnet@plot)
      dev.off()
    }
  })
}
  }, error = function(e) {
    message("Error in subnetwork_enrichment: ", e$message)
  })
}else if(Networkclass=="Stable_DifferentialNetwork"){

  differential_network1=Network

  
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
        readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",net_name1,".txt")),delim="\t")
        readr::write_delim(edges,file.path(outdir1,paste0("edges_",net_name1,".txt")),delim="\t") 
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
    if(any(!is.null(differential_network1@Differential_subnetwork)) &&
       any(!is.null(differential_network1@Differential_subnetwork@overall_cluster_network))
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
  readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",group_name,".txt")),delim="\t")
  readr::write_delim(edges,file.path(outdir1,paste0("edges_",group_name,".txt")),delim="\t")

  diff_net_community_detection_plot(file.path(outdir1,paste0("nodes_",group_name,".txt")),
                                    file.path(outdir1,paste0("edges_",group_name,".txt")),
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
    readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",i,"_",net_name1,".txt")),delim="\t")
    readr::write_delim(edges,file.path(outdir1,paste0("edges_",i,"_",net_name1,".txt")),delim="\t") 
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
  
  tryCatch({
if (any(!is.null(differential_network1@DiffSubnet_Enrichment))){
netname<-names(differential_network1@DiffSubnet_Enrichment@network)

filelist<-lapply(netname, function(i){
  if(any(!is.null(differential_network1@DiffSubnet_Enrichment)) &&
     any(!is.null(differential_network1@DiffSubnet_Enrichment@network))){
    sub_ann_filter <- differential_network1@DiffSubnet_Enrichment@network[[i]]$annotations_filter
  }else{
    sub_ann_filter<-NULL
  }
  if (!is.null(sub_ann_filter) && nrow(sub_ann_filter) > 0) {
    bubbleDiagram_subnet<-network_show(Network=differential_network1,subnetwork_name = c(i),
                                       plot_type="differential_subnetwork_enrichment",
                                       richfactor_threshold = richfactor_threshold,plot_title_size=plot_title_size,
                                       axis_title_size=axis_title_size,
                                       text_size=text_size,
                                       legend_title_size=legend_title_size,
                                       legend_text_size=legend_text_size)
    outdir1<-file.path(outdir,"SubNetwork",i) 
    if (!dir.exists(outdir1)) {
      dir.create(outdir1, recursive = TRUE)
    }
    readr::write_delim(as.data.frame(bubbleDiagram_subnet@data[[i]]), file.path(outdir1, paste0("BubbleDiagram_",i,"_",group_name,".txt")), delim = "\t")
    pdf(file = file.path(outdir1, paste0("BubbleDiagram_",i,"_",group_name,".pdf")), 
        width = bubbleDiagram_subnet@picture_width, height = bubbleDiagram_subnet@picture_height)
    dev.off()
    
  }else{
    print(paste0("subnet: The enrichment result for ",i," is null."))
  }
})
}
  }, error = function(e) {
    message("Error in differential_subnetwork_enrichment: ", e$message)
  })
  
  
}else if(Networkclass=="Stable_MultiplexNetwork"){
 

  
  multiplex_network1=Network
  tryCatch({
  ppiplot<-network_show(Network=multiplex_network1,plot_type="interaction_network",image_margin_size=image_margin_size,
                        plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size)
  outdir1<-file.path(outdir,Networkclass,"OverallNetwork","InteractionNetwork") 
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  nodes<-as.data.frame(ppiplot@data$net$nodes) |>
    dplyr::mutate(x=ppiplot@data$net$plot_layout[,1],y=ppiplot@data$net$plot_layout[,2]) 
  edges<-as.data.frame(ppiplot@data$net$edges)
  readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",group_name,".txt")),delim="\t")
  readr::write_delim(edges,file.path(outdir1,paste0("edges_",group_name,".txt")),delim="\t") 
  pdf(file = file.path(outdir1, paste0("Interaction_Network_",group_name,".pdf")), 
      width = ppiplot@picture_width, height = ppiplot@picture_height)
  print(ppiplot@plot)
  dev.off()
  }, error = function(e) {
    message("Error in interaction_network: ", e$message)
  })
  tryCatch({

  bootnetplot1<-network_show(Network=multiplex_network1,plot_type="case_stable_test",image_margin_size=image_margin_size,
                             plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                             stable_num=4,R_threshold = R_threshold)
  outdir1<-file.path(outdir,Networkclass,"StabilityTest",casename)
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
    readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",net_name,"_",casename,".txt")),delim="\t")
    readr::write_delim(edges,file.path(outdir1,paste0("edges_",net_name,"_",casename,".txt")),delim="\t") 
  })
  pdf(file = file.path(outdir1, paste0("Bootnet_",casename,".pdf")), 
      width = bootnetplot1@picture_width, height = bootnetplot1@picture_height)
  print(bootnetplot1@plot)
  dev.off()
  }, error = function(e) {
    message("Error in case_stable_test: ", e$message)
  })
  tryCatch({

  bootnetplot2<-network_show(Network=multiplex_network1,plot_type="control_stable_test",image_margin_size=image_margin_size,
                             plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                             stable_num=4,R_threshold = R_threshold)
  outdir1<-file.path(outdir,Networkclass,"StabilityTest",controlname)
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
    readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",net_name,"_",controlname,".txt")),delim="\t")
    readr::write_delim(edges,file.path(outdir1,paste0("edges_",net_name,"_",controlname,".txt")),delim="\t") 
  })
  pdf(file = file.path(outdir1, paste0("Bootnet_",controlname,".pdf")), 
      width = bootnetplot2@picture_width, height = bootnetplot2@picture_height)
  print(bootnetplot2@plot)
  dev.off()
}, error = function(e) {
  message("Error in control_stable_test: ", e$message)
})
  
  tryCatch({

  stableplot1<-network_show(Network=multiplex_network1,plot_type="case_overall_network",image_margin_size=image_margin_size,
                            plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                            R_threshold = R_threshold)
  outdir1<-file.path(outdir,Networkclass,"OverallNetwork","CorrelationNetwork",casename) 
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  nodes<-as.data.frame(stableplot1@data$net$nodes) |>
    dplyr::mutate(x=stableplot1@data$net$plot_layout[,1],y=stableplot1@data$net$plot_layout[,2]) 
  edges<-as.data.frame(stableplot1@data$net$edges)
  readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",casename,".txt")),delim="\t")
  readr::write_delim(edges,file.path(outdir1,paste0("edges_",casename,".txt")),delim="\t") 
  pdf(file = file.path(outdir1, paste0("Correlation_Network_",casename,".pdf")), 
      width = stableplot1@picture_width, height = stableplot1@picture_height)
  print(stableplot1@plot)
  dev.off()
  }, error = function(e) {
    message("Error in case_overall_network: ", e$message)
  })
  tryCatch({
  stableplot2<-network_show(Network=multiplex_network1,plot_type="control_overall_network",image_margin_size=image_margin_size,
                            plot_title_size=plot_title_size,font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                            R_threshold = R_threshold)
  outdir1<-file.path(outdir,Networkclass,"OverallNetwork","CorrelationNetwork",controlname) 
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  nodes<-as.data.frame(stableplot2@data$net$nodes) |>
    dplyr::mutate(x=stableplot2@data$net$plot_layout[,1],y=stableplot2@data$net$plot_layout[,2]) 
  edges<-as.data.frame(stableplot2@data$net$edges)
  readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",controlname,".txt")),delim="\t")
  readr::write_delim(edges,file.path(outdir1,paste0("edges_",controlname,".txt")),delim="\t") 
  pdf(file = file.path(outdir1, paste0("Correlation_Network_",controlname,".pdf")), 
      width = stableplot2@picture_width, height = stableplot2@picture_height)
  print(stableplot2@plot)
  dev.off()
  }, error = function(e) {
    message("Error in control_overall_network: ", e$message)
  })
  
  tryCatch({
  if(any(!is.null(multiplex_network1@Differential_multiplexnetwork))){ 
    diffnetplot<-network_show(Network=multiplex_network1,plot_type="differential_network",
                              image_margin_size=image_margin_size,plot_title_size=plot_title_size,
                              font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                              node_colortype="Log2FC",focus=c("all"),node_size=3,node_name_size=node_name_size,
                              show_edge_legend = TRUE,show_node_legend = TRUE,show_node_name = TRUE
                              )
    outdir1<-file.path(outdir,Networkclass,"OverallNetwork")
    if (!dir.exists(outdir1)) {
      dir.create(outdir1, recursive = TRUE)
    }
    for (net_name1 in names(diffnetplot@data)) {
      LAYOUT<-as.data.frame(diffnetplot@data[[net_name1]]$plot_layout)
      nodes<-as.data.frame(diffnetplot@data[[net_name1]]$nodes) |>
        dplyr::mutate(x=LAYOUT[,1],y=LAYOUT[,2]) 
      edges<-as.data.frame(diffnetplot@data[[net_name1]]$edges)
      readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",net_name1,".txt")),delim="\t")
      readr::write_delim(edges,file.path(outdir1,paste0("edges_",net_name1,".txt")),delim="\t") 
    }
    
    pdf(file = file.path(outdir1, paste0("Differential_network_",group_name,".pdf")), 
        width = diffnetplot@picture_width, height = diffnetplot@picture_height)
    print(diffnetplot@plot)
    dev.off()
  }
  }, error = function(e) {
    message("Error in differential_network: ", e$message)
  })
  tryCatch({
    if(any(!is.null(multiplex_network1@Enrichment)) &&
       any(!is.null(multiplex_network1@Enrichment@network))){
      ann_filter <- multiplex_network1@Enrichment@network$annotations_filter
    }else{
      ann_filter <-NULL
    }
    if (!is.null(ann_filter) && nrow(ann_filter) > 0) {

      bubbleDiagram_net<-network_show(Network=multiplex_network1,plot_type="enrichment",focus=c("all"),
                                      richfactor_threshold = richfactor_threshold,plot_title_size=plot_title_size,
                                      axis_title_size=axis_title_size,
                                      text_size=text_size,
                                      legend_title_size=legend_title_size,
                                      legend_text_size=legend_text_size
                                      )
      outdir1<-file.path(outdir,Networkclass,"OverallNetwork")  
      if (!dir.exists(outdir1)) {
        dir.create(outdir1, recursive = TRUE)
      }
      readr::write_delim(bubbleDiagram_net@data$net, file.path(outdir1, paste0("BubbleDiagram_",group_name,".txt")), delim = "\t")
      pdf(file = file.path(outdir1, paste0("BubbleDiagram_",group_name,".pdf")), 
          width = bubbleDiagram_net@picture_width, height = bubbleDiagram_net@picture_height)
      print(bubbleDiagram_net@plot)
      dev.off()
    }
  }, error = function(e) {
    message("Error in enrichment: ", e$message)
  })
  tryCatch({
  if(any(!is.null(multiplex_network1@Differential_subnetwork)) &&
     any(!is.null(multiplex_network1@Differential_subnetwork@overall_cluster_network)) 
     ){
  diff_subnetplot<-network_show(Network=multiplex_network1,plot_type="diff_overall_cluster_network",plot_title_size=plot_title_size,image_margin_size=image_margin_size,
                                focus=c("all"),font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                                show_edge_legend = TRUE,show_node_legend = TRUE,show_node_name = TRUE
                                )
  outdir1<-file.path(outdir,Networkclass,"NetworkClustering")  
  if (!dir.exists(outdir1)) {
    dir.create(outdir1, recursive = TRUE)
  }
  nodes<-as.data.frame(diff_subnetplot@data$net$nodes) |>
    dplyr::mutate(x=diff_subnetplot@data$net$plot_layout[,1],y=diff_subnetplot@data$net$plot_layout[,2]) 
  edges<-as.data.frame(diff_subnetplot@data$net$edges)
  readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",group_name,".txt")),delim="\t")
  readr::write_delim(edges,file.path(outdir1,paste0("edges_",group_name,".txt")),delim="\t") 
  pdf(file = file.path(outdir1, paste0("Network_Clustering_",group_name,".pdf")), 
      width = diff_subnetplot@picture_width, height = diff_subnetplot@picture_height)
  print(diff_subnetplot@plot)
  dev.off()
  }
  }, error = function(e) {
    message("Error in diff_overall_cluster_network: ", e$message)
  })
  
  
  if (any(!is.null(multiplex_network1@Differential_subnetwork))){
    tryCatch({
    netname<-names(multiplex_network1@Differential_subnetwork@subnetworks)
    filelist<- lapply(netname, function(i){
      diff_subnet_top_plot<-network_show(Network=multiplex_network1,
                                         plot_type="differential_subnetwork",image_margin_size=image_margin_size,
                                         node_colortype="Log2FC",focus=c("all"),subnetwork_name = c(i),
                                         show_edge_legend = TRUE,show_node_legend = TRUE,plot_title_size=plot_title_size,
                                         font_family=font_family,legend_title_size=legend_title_size,legend_text_size=legend_text_size,
                                         show_node_name = TRUE,node_name_size=node_name_size
                                         )
      outdir1<-file.path(outdir,Networkclass,"SubNetwork",i) 
      if (!dir.exists(outdir1)) {
        dir.create(outdir1, recursive = TRUE)
      }
      for (net_name1 in names(diff_subnet_top_plot@data)) {
        nodes<-as.data.frame(diff_subnet_top_plot@data[[net_name1]]$nodes) |>
          dplyr::mutate(x=diff_subnet_top_plot@data[[net_name1]]$plot_layout[,1],y=diff_subnet_top_plot@data[[net_name1]]$plot_layout[,2]) 
        edges<-as.data.frame(diff_subnet_top_plot@data[[net_name1]]$edges)
        readr::write_delim(nodes,file.path(outdir1,paste0("nodes_",i,"_",net_name1,".txt")),delim="\t")
        readr::write_delim(edges,file.path(outdir1,paste0("edges_",i,"_",net_name1,".txt")),delim="\t") 
      }
      pdf(file = file.path(outdir1, paste0("Differential_subnetwork_",i,"_",group_name,".pdf")), 
          width = diff_subnet_top_plot@picture_width, height = diff_subnet_top_plot@picture_height)
      print(diff_subnet_top_plot@plot)
      dev.off()

    })
    }, error = function(e) {
      message("Error in differential_subnetwork: ", e$message)
    })
    
  }
  
  tryCatch({
  if (any(!is.null(multiplex_network1@DiffSubnet_Enrichment))){
  netname<-names(multiplex_network1@DiffSubnet_Enrichment@network)
  # enrich_subnets<- c()
  filelist<-lapply(netname, function(i){
    if(any(!is.null(multiplex_network1@DiffSubnet_Enrichment)) &&
       any(!is.null(multiplex_network1@DiffSubnet_Enrichment@network))){
      ann_filter <- multiplex_network1@DiffSubnet_Enrichment@network[[i]]$annotations_filter
    }else{
      ann_filter<-NULL
    }
    if (!is.null(ann_filter) && nrow(ann_filter) > 0) {
      bubbleDiagram_subnet<-network_show(Network=multiplex_network1,subnetwork_name = c(i),
                                         plot_type="differential_subnetwork_enrichment",
                                         richfactor_threshold = richfactor_threshold,plot_title_size=plot_title_size,
                                         axis_title_size=axis_title_size,
                                         text_size=text_size,
                                         legend_title_size=legend_title_size,
                                         legend_text_size=legend_text_size
                                        )
      outdir1<-file.path(outdir,Networkclass,"SubNetwork",i) 
      if (!dir.exists(outdir1)) {
        dir.create(outdir1, recursive = TRUE)
      }
      readr::write_delim(as.data.frame(bubbleDiagram_subnet@data[[i]]), file.path(outdir1, paste0("BubbleDiagram_",i,"_",group_name,".txt")), delim = "\t")
      pdf(file = file.path(outdir1, paste0("BubbleDiagram_",i,"_",group_name,".pdf")), 
          width = bubbleDiagram_subnet@picture_width, height = bubbleDiagram_subnet@picture_height)
      print(bubbleDiagram_subnet@plot)
      dev.off()
    }
  })
    }
  }, error = function(e) {
    message("Error in differential_subnetwork_enrichment: ", e$message)
  })
}else{
  stop(paste0("‘pipline_save‘ does not applies to the target object :",Networkclass))
}
}