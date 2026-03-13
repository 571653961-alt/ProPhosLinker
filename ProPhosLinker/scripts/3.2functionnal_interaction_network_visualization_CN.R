#' @title Visualize a phosphoprotein network
#' @description This function visualizes a phosphoprotein network using the igraph package.param name 

suppressWarnings(library(optparse))   # 加载 optparse 包


# args <- commandArgs(trailingOnly=TRUE)
# print(paste("参数1:", args[1]))
# print(paste("参数2:", args[2]))


option_list <- list(
  make_option( c( "--function_enrich_rscript_path"),type = "character",default = NULL,help = "富集代码路径",metavar = "FILE"),
  #input parameter
  make_option(c("-n", "--nodes"), type = "character", default = NULL, metavar = "FILE",
              help = "Input nodes file (CSV)"),
  make_option(c("-e", "--edges"), type = "character", default = NULL, metavar = "FILE",
              help = "Input edges file (CSV)"),
  make_option(c("-m", "--module_info"), type = "character", default = NULL, metavar = "FILE",
              help = "Module information file (CSV)"),
  make_option(c("-i", "--module_id"), type = "integer", default = 0, metavar = "INT",
              help = "Module ID to visualize [default: %default]"),
  #output parameter
  make_option(c("-o", "--outdir"), type = "character", default = getwd(), metavar = "DIR",
              help = "Output directory [default: %default]"),
  make_option(c("--enrich_fromType"),type = "character",default = "UNIPROT",
              help = "P值类型 [默认: %default], 选项: UNIPROT, SYMBOL",metavar = "STRING"),
  #network visualization parameter
  make_option(c("-d", "--max_phosphoSite_displayed"), type = "integer", default = 5, metavar = "INT",
              help = "Max phosphosites per protein [default: %default]"),
  make_option(c("-t", "--network_mode"), type = "character", default = "ALL", metavar = "METHOD",
              help = "Subnetwork filter: ALL, KS, LR, TF, DEP, DEPUP, DEPDOWN, SEEDNODE [default: %default]"),
  make_option(c("-s", "--SEEDNODEID"), type = "character", default = NULL, metavar = "STRING",
              help = "SEEDNODE ID [default: %default]"),
  make_option(c("-f", "--node_filtering"), type = "character", default = "betweenness", metavar = "STRING",
              help = "Node filter: degree, betweenness, closeness, harmonic, eigenvalue, pagerank, alpha, hub, authority.[default: %default]"),
  make_option(c("-l", "--network_layout"), type = "character", default = "kk", metavar = "STRING",
              help = "Layout method: fr, kk, dh, stress, tree, gem, graphopt, lgl, circle, grid. [default: %default]"),
  make_option(c("-p", "--top_nodes_visualization_num"), type = "integer", default = 20, metavar = "INT",
              help = "Specifies the number of top nodes to visualize in the network graph."),
  # Omics names
  make_option(c("--omics1_name"), type = "character", default = "Proteomics", metavar = "STR",
              help = "Name for first omics dataset [default: %default]"),
  make_option(c("--omics2_name"), type = "character", default = "Phosphoproteomics", metavar = "STR",
              help = "Name for second omics dataset [default: %default]"),
  # color parameters
  make_option(c("--node_up"), type = "character", default = "#B83A2D", metavar = "character",
              help = "Color for upregulated nodes [default: %default]"),
  make_option(c("--node_down"), type = "character", default = "#657f68", metavar = "character",
              help = "Color for downregulated nodes [default: %default]"),
  make_option(c("--node_nonsig"), type = "character", default = "#E3C79F", metavar = "character",
              help = "Color for non-significant nodes [default: %default]"),
  make_option(c("--node_notdet"), type = "character", default = "grey90", metavar = "character",
              help = "Color for not detected nodes [default: %default]"),
  make_option(c("--node_disease_border_color"), type = "character", default = "#535c54", metavar = "character",
              help = "Border color for disease-related nodes [default: %default]"),
  make_option(c("--node_notdisease_border_color"), type = "character", default = "#acabab", metavar = "character",
              help = "Border color for non-disease nodes [default: %default]"),
  make_option(c("--function_enrich_color_gradient_low"), type = "character", default = "#175663", metavar = "character",
              help = "Low end color for functional enrichment gradient [default: %default]"),
  make_option(c("--function_enrich_color_gradient_high"), type = "character", default = "#90362d", metavar = "character",
              help = "high end color for functional enrichment gradient [default: %default]"),
  make_option(c("--edge_type_colors"), type = "character",
              default = "association:#d7d6d6;physical association:#838181;binding:#838181;direct interaction:#d7d6d6;activation:#d62c0c;catalysis:#7f00ff;proximity:#FF9301;reaction:#114335;phosphorylation:#2FBE95;dephosphorylation:#6B8E23;ptmod:#8C97D6;inhibition:#0cb6d6;expression:#FCF402;regulation:#e4dadd;colocalization:#4c95cd;covalent binding:#716F74;ubiquitination:#FF4500;multiRel:#9a6728",
              metavar = "character",
              help = "Edge type colors as semicolon-separated key:value pairs [default: %default]"),
  #general parameter
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE, metavar = "FLAG",
              help = "Print detailed output messages [Default: %default]")

)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                           description = "Visualize a protein-phosphoprotein network using igraph package.")
opt <- parse_args(opt_parser)



# #################temporate test########################
# # 模拟命令行参数
# opt <- list()
# opt$nodes <- "E:/pro_phos_test2/results/0.Data_Preprocessing/../3.Functional_Analysis/3.2KGbased_Functional_Network/Node_info.tsv"
# opt$edges <- "E:/pro_phos_test2/results/0.Data_Preprocessing/../3.Functional_Analysis/3.2KGbased_Functional_Network/Edge_info.tsv"
# opt$module_info <- "E:/pro_phos_test2/results/0.Data_Preprocessing/../3.Functional_Analysis/3.2KGbased_Functional_Network/DEP_community_detection/infomap_community_detection_modules.tsv"
# opt$module_id <- 11
# opt$outdir <- "E:/pro_phos_test2/results/0.Data_Preprocessing/../3.Functional_Analysis/3.2KGbased_Functional_Network/DEP_community_detection"
# opt$max_phosphoSite_displayed <- 5
# opt$network_mode <- "ALL"
# opt$SEEDNODEID <- "None"
# opt$node_filtering <- "betweenness"
# opt$network_layout <- "kk"
# opt$top_nodes_visualization_num <- 20
# opt$node_up <- "#ff7f00"
# opt$node_down <- "#012da7"
# opt$node_nonsig <- "#E3C79F"
# opt$node_notdet <- "grey90"
# opt$node_disease_border_color <- "#B83A2D"
# opt$node_notdisease_border_color <- "#bec389"
# opt$function_enrich_color_gradient_low <- "#012da7"
# opt$function_enrich_color_gradient_high <- "#ff7f00"
# opt$edge_type_colors <- "association:#d7d6d6;physical association:#838181;binding:#838181;direct interaction:#d7d6d6;activation:#d62c0c;catalysis:#7f00ff;proximity:#FF9301;reaction:#114335;phosphorylation:#2FBE95;dephosphorylation:#6B8E23;ptmod:#8C97D6;inhibition:#0cb6d6;expression:#FCF402;regulation:#e4dadd;colocalization:#4c95cd;covalent binding:#716F74;ubiquitination:#FF4500;multiRel:#9a6728"
# opt$function_enrich_rscript_path <- "E:/pro_phos_V2/scripts/3.2functional_enrichment_function.R"
# 
# #################temporate test########################


#parameter validation
parameter_validation <- function(opt){
  #path cleaning
  clean_path <- function(path) {
    # Trim leading and trailing whitespace and quotes
    path <- gsub("^['\" ]+|['\" ]+$","",path)
    # Replace backslashes (\) with forward slashes (/)
    path <- gsub("\\\\","/",path)
    return(path)
  }
  check_input_file <- function(file_path,file_name){
    if(!file.exists(file_path)){
      stop("❌ Error: The file '", file_name, "' does not exist at the specified path: ", file_path)
    }
    if(file.size(file_path) == 0){
      stop("❌ Error: The file '", file_name, "' is empty at the specified path: ", file_path)
    }
  }
  check_output_dir <- function(result_dir){
    result_dir <- paste0(result_dir,"/",tools::file_path_sans_ext(basename(opt$module_info)),'_',opt$module_id,'_visualization')
    if (!dir.exists(result_dir)) {
      # Recursively create directory and all parent directories
      dir.create(result_dir, recursive = TRUE)
      print(paste("✅ Output directory created:", result_dir))
    } else {
      print(paste("✅ Output directory already exists:", result_dir))
    }

    return(result_dir)
  }
  validate_option_choises <- function(value, allowed_values, option_name) {
    if (!value %in% allowed_values) {
      stop("❌ Error: Invalid", option_name, "value '", value, "'. Please use one of the following options: ",
           paste(allowed_values, collapse = ", "))
    }
  }
  validate_numeric_range <- function(value, min_val, max_val, option_name) {
    if (is.null(value)) {
      return()  # Skip validation if the value is NULL
    }
    if (value < min_val || value > max_val) {
      stop("❌ Error: ", option_name, " must be within the range [", min_val, ", ", max_val, "].")
    }
  }
  
  # 解析 edge_type_colors 函数
  parse_edge_colors <- function(color_str) {
    color_str <- gsub('"','',color_str)
    pairs <- strsplit(color_str, ";")[[1]]
    color_list <- list()
    for(pair in pairs) {
      key_val <- strsplit(pair, ":")[[1]]
      if(length(key_val) == 2) {
        color_list[[key_val[1]]] <- key_val[2]
      }
    }
    return(color_list)
  }
  
  
  # clean all pathway parameter
  opt$nodes <- clean_path(opt$nodes)
  opt$edges <- clean_path(opt$edges)
  opt$module_info <- clean_path(opt$module_info)
  opt$function_enrich_rscript_path <- clean_path(opt$function_enrich_rscript_path)
  opt$outdir <- clean_path(opt$outdir)

  #check input files
  check_input_file(opt$nodes,"nodes")
  check_input_file(opt$edges,"edges")
  check_input_file(opt$module_info,"module info")
  opt$outdir  <- check_output_dir(opt$outdir)


  #choices validation

  allowed_network_mode <- c("ALL",    # 3 catalogs
                            "KS","LR","TF",
                            "DEP", "DEPUP", "DEPDOWN", "SEEDNODE")
  allowed_node_filtering <- c("betweenness", "degree",  "closeness", "harmonic","eigenvector", "pagerank", "alpha", "hub",  "authority")
  allowed_network_layout_type <- c("fr","kk","dh","stress","tree","gem","graphopt", "lgl","circle","grid")
  allowed_enrich_fromType <- c('UNIPROT', 'SYMBOL')
  

  validate_option_choises(opt$network_mode,allowed_network_mode,"network mode parameter")
  validate_option_choises(opt$node_filtering,allowed_node_filtering,"node filtering parameter")
  validate_option_choises(opt$network_layout,allowed_network_layout_type,"network layout parameter")
  validate_option_choises(opt$enrich_fromType,allowed_enrich_fromType,"allowed enrich_fromType parameter")
  validate_numeric_range(opt$max_phosphoSite_displayed, 1,5, "max_phosphoSite_displayed")
  validate_numeric_range(opt$top_nodes_visualization_num, 3,100, "top_nodes_visualization_num")

  opt$SEEDNODEID <- clean_path(opt$SEEDNODEID)
  
  opt$edge_type_colors <- parse_edge_colors(opt$edge_type_colors)

  if (opt$verbose) { print(opt)}

  return(opt)
}

# 设置警告级别为忽略
options(warn = -1)



suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(graphlayouts)))
suppressWarnings(suppressPackageStartupMessages(library(ggforce)))
suppressWarnings(suppressPackageStartupMessages(library(scatterpie)))
suppressWarnings(suppressPackageStartupMessages(library(ggsci)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(tidygraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ggnewscale)))  # 用于处理多个颜色尺度
suppressWarnings(suppressPackageStartupMessages(library(svglite)))
# source("E:/pro_phos_V2/scripts/3.2functional_enrichment_function.R")
source(opt$function_enrich_rscript_path)


###parameter test
# opt <-NULL
# opt$nodes <- "E:/pro_phosphpro/results/Network_analysis/Node_info.tsv"
# opt$edges <- "E:/pro_phosphpro/results/Network_analysis/Edge_info.tsv"
# opt$module_info <- "E:/pro_phosphpro/results/Network_analysis/Community_detection/infomap_community_detection_modules.tsv"
# opt$module_id <- 1
# opt$max_phosphoSite_displayed <- 5
# opt$network_mode <- "ALL"
# opt$node_filtering <- "betweenness"
# opt$network_layout <- "kk"
# opt$outdir <- "E:/pro_phosphpro/results/Network_analysis/Community_detection/infomap_community_detection_modules_1_visualization"
top_nodes_visualization_num <- opt$top_nodes_visualization_num # Number of top nodes to visualize
# opt$network_mode <- 'SEEDNODE'
# opt$SEEDNODEID <- 'CFTR'


# subnetwork construction
subnetwork_construction <- function(opt,top_nodes_visualization_num){
  setwd(opt$outdir)
  ########################0. data loading#################################
  f_Node_path <- opt$nodes
  f_Edge_path <- opt$edges
  f_Module_info_path <- opt$module_info
  f_Module_id <- opt$module_id
  network_mode <- opt$network_mode
  SEEDNODEID <- opt$SEEDNODEID

  #
  node_filtering <- opt$node_filtering

  # Load nodes and edges
  f_Node_rawdata <- read.csv(f_Node_path,sep = '\t', stringsAsFactors = FALSE)
  f_Edge_rawdata <- read.csv(f_Edge_path,sep = '\t', stringsAsFactors = FALSE)
  f_Module_info <- read.csv(f_Module_info_path,sep = '\t',stringsAsFactors = FALSE)
  # print(f_Node_rawdata[1:6,1:6])
  # print(f_Edge_rawdata[1:6,])
  # print(colnames(f_Module_info))
  module_nodesID_list <- strsplit(f_Module_info[f_Module_info$ModuleID == f_Module_id,'ModuleItems'],';')[[1]]
  f_Node <- f_Node_rawdata[f_Node_rawdata$NodeName %in% module_nodesID_list,]
  # print(length(unique(f_Node$NodeName)))
  if(network_mode =='SEEDNODE'){
    if(length(SEEDNODEID) == 0){
      return("❌Error: SEEDNODEID value is empty!")
      # stop("❌Error: SEEDNODEID value is empty!")
    }else if(!any(SEEDNODEID %in% module_nodesID_list)){
      return(paste0("❌ Error: SEEDNODEID ",SEEDNODEID ," not in this module."))
      # stop(paste0("❌ Error: SEEDNODEID ",SEEDNODEID ," not in this module."))
    }
  }
  ########################1. nodes filtering#########################################
  # coping with the max_phosphoSite_num_displayed and DEP relative parameter
  max_phosphoSite_num_displayed <- opt$max_phosphoSite_displayed

  if( network_mode == 'DEP'){
    filter_node_data <-  f_Node[(f_Node$Phosphopro_class == 'UP' | f_Node$Phosphopro_class == 'DOWN') &(f_Node$Pro_class == 'UP'|f_Node$Pro_class == 'DOWN'),]%>%
      # Group processing
      group_by(NodeName) %>%
      # Calculate absolute values and sort
      arrange(
        desc(abs(Phosphopro_FC)),             # Primary priority: FC by absolute value (descending)
        .by_group = TRUE) %>%
      # Add row numbers for filtering
      mutate(phosphoSite_filter_index = row_number()) %>%
      # Filter top N (max_phosphoSite_num_displayed) rows per group
      filter(phosphoSite_filter_index <= max_phosphoSite_num_displayed) %>%
      # Ungroup the data
      ungroup()
    # print(filter_node_data[filter_node_data$Pro_class == 'UP',])
  }else if(network_mode == 'DEPUP'){
    filter_node_data <-  f_Node[f_Node$Phosphopro_class == 'UP' & f_Node$Pro_class == 'UP',]%>%
      # Group processing
      group_by(NodeName) %>%
      # Calculate absolute values and sort
      arrange(
        desc(abs(Phosphopro_FC)),             # Primary priority: FC by absolute value (descending)
        .by_group = TRUE) %>%
      # Add row numbers for filtering
      mutate(phosphoSite_filter_index = row_number()) %>%
      # Filter top N (max_phosphoSite_num_displayed) rows per group
      filter(phosphoSite_filter_index <= max_phosphoSite_num_displayed) %>%
      # Ungroup the data
      ungroup()
    # print(filter_node_data[filter_node_data$Pro_class == 'UP',])
  }else if(network_mode == 'DEPDOWN'){
    filter_node_data <-  f_Node[f_Node$Phosphopro_class == 'DOWN' & f_Node$Pro_class == 'DOWN',]%>%
      # Group processing
      group_by(NodeName) %>%
      # Calculate absolute values and sort
      arrange(
        desc(abs(Phosphopro_FC)),             # Primary priority: FC by absolute value (descending)
        .by_group = TRUE) %>%
      # Add row numbers for filtering
      mutate(phosphoSite_filter_index = row_number()) %>%
      # Filter top N (max_phosphoSite_num_displayed) rows per group
      filter(phosphoSite_filter_index <= max_phosphoSite_num_displayed) %>%
      # Ungroup the data
      ungroup()
  }else{
    filter_node_data <- f_Node %>%
      # Group processing
      group_by(NodeName) %>%
      # Calculate absolute values and sort
      arrange(
        desc(abs(Phosphopro_FC)),             # Primary priority: FC by absolute value (descending)
        .by_group = TRUE) %>%
      # Add row numbers for filtering
      mutate(phosphoSite_filter_index = row_number()) %>%
      # Filter top N (max_phosphoSite_num_displayed) rows per group
      filter(phosphoSite_filter_index <= max_phosphoSite_num_displayed) %>%
      # Ungroup the data
      ungroup()

  }

  if(all(is.na(filter_node_data$NodeName))){
    return(" ❌ Error: This module does not contain any nodes with the specified filtering criteria.")
    # stop(" ❌ Error: This module does not contain any nodes with the specified filtering criteria.")
  }
  ########################2. use filtered nodes to construct a graph#########################
  nodes <- unique(filter_node_data$NodeName)

  # print('???????????????????')
  # print(length(nodes))
  edges <- f_Edge_rawdata %>%
    filter(Source %in% nodes,Target %in% nodes)
  # dealing with edges

  clean_interaction_type <- function(types) {
    result <- types
    # Helper function: Safe replacement (skip if the result is empty after substitution)
    safe_replace <- function(pattern, replacement, x) {
      new_x <- gsub(pattern, replacement, x)
      # 只有当新值不为空时才替换
      ifelse(new_x == "", x, new_x)
    }
    result <- safe_replace("(^|;)association(;|$)", "\\1\\2", result)
    result <- safe_replace("/association", "", result)

    # Condition 1: Handling 'physical association'
    result <- ifelse(
      grepl("physical association", result) & grepl(";", result),
      safe_replace(";physical association|physical association;", "", result),
      result
    )

    # Condition 2: Handling 'binding'
    result <- ifelse(
      grepl("binding", result) & grepl(";", result),
      safe_replace(";binding|binding;", "", result),
      result
    )

    # Condition 3: Handling 'direct interaction'
    result <- ifelse(
      grepl("direct interaction", result) & grepl(";", result),
      safe_replace(";direct interaction|direct interaction;", "", result),
      result
    )

    # Condition 4: Handling 'expression' coexisting with 'activation' or 'inhibition'
    has_expression <- grepl("expression", result)
    has_activation <- grepl("activation", result)
    has_inhibition <- grepl("inhibition", result)

    result <- ifelse(
      has_expression & (has_activation | has_inhibition),
      safe_replace(";expression|expression;", "", result),
      result
    )


    # Condition 1: Handling 'expression' in the presence of 'activation' or 'inhibition'
    has_reaction <- grepl("(^reaction;)|(;reaction;)|(;reaction$)", result)
    has_catalysis <- grepl("catalysis", result)

    result <- ifelse(
      has_reaction & has_catalysis,
      safe_replace(";reaction|reaction;", "", result),
      result
    )

    # 清理多余的分号
    result <- safe_replace("^;|;$", "", result)  # Remove leading and trailing semicolons
    result <- safe_replace(";;+", ";", result)   # Replace multiple consecutive semicolons with a single semicolon

    # Ensure that an empty string is not returned - if empty, restore the original value
    result <- ifelse(result == "", types, result)

    return(result)

  }

  edges$Interaction_Type_cleaned <- clean_interaction_type(edges$Interaction_Type)
  edges$Interaction_Type_simplified <- ifelse(
    grepl(';',edges$Interaction_Type_cleaned),
    "multiRel",
    edges$Interaction_Type_cleaned
  )
  # Preliminary network construction
  g <- graph_from_data_frame(
    d = edges,
    directed = TRUE,
    vertices = nodes
  )
  # Attach edge attributes
  E(g)$Interaction_Type <- edges$Interaction_Type
  E(g)$Interaction_Type_cleaned <- edges$Interaction_Type_cleaned
  E(g)$Interaction_Type_simplified <- edges$Interaction_Type_simplified

  ########################3. filtering edges from the graph#########################
  # network_mode == 'ALL'    # No action taken

  if(network_mode == 'KS'){
    edges_phospho <- E(g)[grepl("phospho", E(g)$Interaction_Type)]
    g <- subgraph_from_edges(
      graph = g,
      eids = edges_phospho,  # Preserved edges
      delete.vertices = TRUE  # Automatically delete disconnected nodes
    )
  }else if(network_mode == 'LR'){
    edges_LR <- E(g)[grepl("binding", E(g)$Interaction_Type)]
    E(g)$Interaction_Type_simplified <- ifelse(E(g)$Interaction_Type_simplified == 'binding', 'binding','multiRel')
    g <- subgraph_from_edges(
      graph = g,
      eids = edges_LR,
      delete.vertices = TRUE
    )

  }else if(network_mode == 'TF'){
    edges_TF <- E(g)[
      grepl(
        pattern = "activation|regulation|repression|inhibition",
        x = E(g)$Interaction_Type,
        ignore.case = TRUE  # 忽略大小写（可选）
      )
    ]
    g <- subgraph_from_edges(
      graph = g,
      eids = edges_TF,
      delete.vertices = TRUE
    )
  }else if(network_mode == 'SEEDNODE'){
    focus_node_edges <- E(g)[.inc(V(g)[name == SEEDNODEID])]
    connected_vertices <- ends(g,focus_node_edges,names =TRUE)
    connected_vertices <- unique(as.vector((connected_vertices)))
    if(length(connected_vertices) >=10){
      final_node <- connected_vertices
    }else{
      # the second layer
      final_node <- NULL
      for (i in connected_vertices) {
        node_edges <- E(g)[.inc(V(g)[name == i ])]
        vertices <- ends(g,node_edges,names = TRUE)
      }
      final_node <- union(final_node,vertices)
    }

    g <- induced_subgraph(
      graph = g,
      vids = final_node
    )
  }

  if(length(V(g)) == 1){
    return(" ❌ Error: This module only contains one node with the specified filtering criteria.")
    # stop(" ❌ Error: This module only contains one node with the specified filtering criteria.")
  }


  ########################4. filtering nodes to visualise through centrality#########################
  if(node_filtering == 'betweenness'){
    # Calculate betweenness centrality
    values <- betweenness(g, directed = FALSE, normalized = FALSE)
  }else if(node_filtering == 'degree'){
    # Calculate degree centrality
    values <- degree(g, mode = "all", normalized = FALSE)
  }else if(node_filtering == 'closeness'){
    # Calculate closeness centrality
    values <- closeness(g, mode = "all", normalized = FALSE)
  }else if(node_filtering == 'harmonic'){
    # Calculate harmonic centrality
    values <- harmonic_centrality(g, mode = "all", normalized = FALSE)
  }else if(node_filtering == 'eigenvector'){
    # Calculate eigenvector centrality
    values <- evcent(g, directed = FALSE)$vector
  }else if(node_filtering == 'pagerank'){
    # Calculate PageRank centrality
    values <- page_rank(g, directed = FALSE)$vector
  }else if(node_filtering == 'alpha'){
    # Calculate alpha centrality
    values <- alpha_centrality(g, alpha = 0.85)
  }else if(node_filtering == 'hub'){
    # Calculate global clustering coefficient
    values <- hub_score(g)$vector
  }else if(node_filtering == 'authority'){
    # Calculate local clustering coefficient
    values <- authority_score(g)$vector
  }
  top_nodes_index <- order(values, decreasing = TRUE)[1:min(top_nodes_visualization_num,vcount(g))]  # Select top 100 nodes by degree centrality
  top_nodes_ids = V(g)$name[top_nodes_index]
  # print(top_nodes_ids)
  # print(vcount(g))
  # # find the top nodes' neighbors
  # all_neighbors <- unique(unlist(lapply(top_nodes_index, function(node) {
  #   neighbors(g, node, mode = "all")  #  Get all neighbors (in-degree and out-degree)
  # })))
  # neighbor_names <- V(g)$name[all_neighbors]
  # print(neighbor_names)
  # subgraph_nodes <- unique(c(top_nodes_ids, neighbor_names))
  # print(subgraph_nodes)
  # Create a subgraph with the top nodes and their neighbors

  ########################error
  if(all(is.na(top_nodes_ids))){
    return(" ❌ Error: This module does not contain any nodes with the specified filtering criteria.")
    # stop(" ❌ Error: This module does not contain any nodes with the specified filtering criteria.")
  }


  topnodes_subgraph <- induced_subgraph(
    g,
    vids = top_nodes_ids
  )
  # 删除孤立节点（度为0的节点）
  if(!(network_mode %in% c("DEP","DEPUP","DEPDOWN"))){
    topnodes_subgraph <- delete_vertices(
    topnodes_subgraph,
    v = which(degree(topnodes_subgraph, mode = "all") == 0)
  )
  }
  # print(V(topnodes_subgraph)$name)

  ########################5. adding nodes attributes#########################
  # Initialize vertex attributes
  vertex_count <- vcount(topnodes_subgraph)
  # Initialize disease_related attributes
  V(topnodes_subgraph)$DiseaseRelated <- rep("no", vertex_count)
  # Define all attribute prefixes
  prefixes <- c("", paste0("PhosphoPro_site", LETTERS[1:5],'_'))
  suffixes <- c("UP", "DOWN", "Non_significant", "Not_Detected")
  # Initialize each attribute combination
  for (prefix in prefixes) {
    for (suffix in suffixes) {
      attr_name <- paste0(prefix, suffix)
      # 正确使用 $ 操作符设置属性
      V(topnodes_subgraph)$attr_name <- rep(0, vertex_count)
      # 修复属性名引用问题
      names(vertex_attr(topnodes_subgraph))[length(vertex_attr(topnodes_subgraph))] <- attr_name
    }
  }
  # print(vertex_attr(topnodes_subgraph))  # View node attributes and their values

  # Assign values to these attributes after initialization
  # print(V(topnodes_subgraph)$name)
  for(nodeName_index in 1:length(V(topnodes_subgraph)$name)){
    # print(nodeName_index)
    nodeName = V(topnodes_subgraph)$name[nodeName_index]
    # print(nodeName)
    node_phosphoprot_info <- filter_node_data[filter_node_data$NodeName  == nodeName,]
    # print(names(node_phosphoprot_info))
    # print(node_phosphoprot_info)
    node_protein_class <- filter_node_data[filter_node_data$NodeName  == nodeName,]$Pro_class[1]
    Node_DiseaseRelated <- filter_node_data[filter_node_data$NodeName  == nodeName,]$DiseaseRelated[1]
    V(topnodes_subgraph)$DiseaseRelated[nodeName_index] = Node_DiseaseRelated
    # print(node_protein_class)    #DOWN
    # print(node_phosphoprot_info$phosphoSite_filter_index)   #int

    for(phosphoSite_filter_index in node_phosphoprot_info$phosphoSite_filter_index){
      # print(node_phosphoprot_info[node_phosphoprot_info$phosphoSite_filter_index == phosphoSite_filter_index,])
      phosphoprotein_ID = node_phosphoprot_info[node_phosphoprot_info$phosphoSite_filter_index == phosphoSite_filter_index,]$phosphoprotein_ID
      Phosphopro_type = node_phosphoprot_info[node_phosphoprot_info$phosphoSite_filter_index == phosphoSite_filter_index,]$Phosphopro_type
      Phosphopro_class = node_phosphoprot_info[node_phosphoprot_info$phosphoSite_filter_index == phosphoSite_filter_index,]$Phosphopro_class
      # Phosphopro_prodetectd_in_module = node_phosphoprot_info[node_phosphoprot_info$phosphoSite_filter_index == phosphoSite_filter_index,]$Phosphopro_prodetectd_in_module
      # print(Phosphopro_class)
      if (is.na(Phosphopro_class) || Phosphopro_class == ""){
        Phosphopro_class = "None"
      }
      # print(phosphoSite_filter_index)
      # print(print(phosphoSite_filter_index) == 1)
      # print( Phosphopro_type )

      if(phosphoSite_filter_index == 1){
        # print('here!!!111111')
        if( Phosphopro_class == 'UP'){
          V(topnodes_subgraph)$PhosphoPro_siteA_UP[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'DOWN'){
          V(topnodes_subgraph)$PhosphoPro_siteA_DOWN[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'Non-significant'){
          # print('here!!!22222')
          V(topnodes_subgraph)$PhosphoPro_siteA_Non_significant[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'None'){
          V(topnodes_subgraph)$PhosphoPro_siteA_Not_Detected[nodeName_index] <- 1
        }
      }
      if(phosphoSite_filter_index == 2){
        if( Phosphopro_class == 'UP'){
          V(topnodes_subgraph)$PhosphoPro_siteB_UP[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'DOWN'){
          V(topnodes_subgraph)$PhosphoPro_siteB_DOWN[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'Non-significant'){
          V(topnodes_subgraph)$PhosphoPro_siteB_Non_significant[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'None'){
          V(topnodes_subgraph)$PhosphoPro_siteB_Not_Detected[nodeName_index] <- 1
        }
      }
      if(phosphoSite_filter_index == 3){
        if( Phosphopro_class == 'UP'){
          V(topnodes_subgraph)$PhosphoPro_siteC_UP[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'DOWN'){
          V(topnodes_subgraph)$PhosphoPro_siteC_DOWN[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'Non-significant'){
          V(topnodes_subgraph)$PhosphoPro_siteC_Non_significant[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'None'){
          V(topnodes_subgraph)$PhosphoPro_siteC_Not_Detected[nodeName_index] <- 1
        }
      }
      if(phosphoSite_filter_index == 4){
        if( Phosphopro_class == 'UP'){
          V(topnodes_subgraph)$PhosphoPro_siteD_UP[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'DOWN'){
          V(topnodes_subgraph)$PhosphoPro_siteD_DOWN[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'Non-significant'){
          V(topnodes_subgraph)$PhosphoPro_siteD_Non_significant[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'None'){
          V(topnodes_subgraph)$PhosphoPro_siteD_Not_Detected[nodeName_index] <- 1
        }
      }
      if(phosphoSite_filter_index == 5){
        if( Phosphopro_class == 'UP'){
          V(topnodes_subgraph)$PhosphoPro_siteE_UP[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'DOWN'){
          V(topnodes_subgraph)$PhosphoPro_siteE_DOWN[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'Non-significant'){
          V(topnodes_subgraph)$PhosphoPro_siteE_Non_significant[nodeName_index] <- 1
        }else if ( Phosphopro_class == 'None'){
          V(topnodes_subgraph)$PhosphoPro_siteE_Not_Detected[nodeName_index] <- 1
        }
      }

    }


    if( node_protein_class == 'UP'){
      V(topnodes_subgraph)$UP[nodeName_index] <-
        V(topnodes_subgraph)$PhosphoPro_siteA_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_Not_Detected[nodeName_index]
    }else if ( node_protein_class == 'DOWN'){
      V(topnodes_subgraph)$DOWN[nodeName_index] <-
        V(topnodes_subgraph)$PhosphoPro_siteA_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_Not_Detected[nodeName_index]
    }else if ( node_protein_class == 'Non-significant'){

      V(topnodes_subgraph)$Non_significant[nodeName_index] <-
        V(topnodes_subgraph)$PhosphoPro_siteA_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_Not_Detected[nodeName_index]
    }else if ( node_protein_class == 'None'){
      V(topnodes_subgraph)$Not_Detected[nodeName_index] <-
        V(topnodes_subgraph)$PhosphoPro_siteA_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_UP[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_DOWN[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_Non_significant[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteA_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteB_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteC_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteD_Not_Detected[nodeName_index]+
        V(topnodes_subgraph)$PhosphoPro_siteE_Not_Detected[nodeName_index]
    }
  }
  # print(vertex_attr(topnodes_subgraph)) #查看节点属性和值
  return(topnodes_subgraph)
}

#network visualization
network_visualiztion <-function(g,opt){
  layout_type <- opt$network_layout
  network_mode <- opt$network_mode
  SEEDNODEID <- opt$SEEDNODEID
  if(layout_type == 'fr'){
    xy <- layout_with_fr(topnodes_subgraph)
  }else if(layout_type == 'kk'){
    xy <- layout_with_kk(topnodes_subgraph)
  }else if(layout_type == 'dh'){
    xy <- layout_with_dh(topnodes_subgraph)
  }else if(layout_type == 'stress'){
    xy <- layout_with_stress(topnodes_subgraph)
  }else if(layout_type == 'gem'){
    xy <- layout_with_gem(topnodes_subgraph) ###【但是我调好了】 没有箭头
  } else if(layout_type == 'graphopt'){
    xy <- layout_with_graphopt(topnodes_subgraph) ### 【但是我调好了】没有箭头
  }else if(layout_type == 'tree'){
    xy <- layout_as_tree(topnodes_subgraph)   ##很压缩，如果可以修改的话？？？【将个烂就】
  }else if(layout_type == 'lgl'){
    xy <- layout_with_lgl(topnodes_subgraph)  ##效果很好
  }else if(layout_type == 'circle'){
    xy <- layout_in_circle(topnodes_subgraph) ##还行吧
  }else if(layout_type == 'grid'){
    xy <- layout_on_grid(topnodes_subgraph)
  }
  # Assign layout coordinates to node attributes x and y
  V(topnodes_subgraph)$x <- xy[, 1]
  V(topnodes_subgraph)$y <- xy[, 2]
  #  Extract node and edge data
  nodes <- as_data_frame(topnodes_subgraph, what = "vertices")
  edges <- as_data_frame(topnodes_subgraph, what = "edges")
  # Add a column for border style to the node data about disease-related info
  nodes <- nodes %>%
    # mutate(node_border_color = ifelse(DiseaseRelated == "yes", "#535c54", "#acabab"),
    mutate(node_border_color = ifelse(DiseaseRelated == "yes", opt$node_disease_border_color, opt$node_notdisease_border_color),
           node_border_linetype = ifelse(DiseaseRelated == "yes", "longdash","solid"),
    )
  # Add coordinate information to the edge data
  edges <- edges %>%
    dplyr::left_join(nodes %>% dplyr::select(name, x, y), by = c("from" = "name")) %>%
    dplyr::rename(x_start = x, y_start = y) %>%
    dplyr::left_join(nodes %>% dplyr::select(name, x, y), by = c("to" = "name")) %>%
    dplyr::rename(x_end = x, y_end = y)
  #  Define column names for all pie chart segments
  pie_cols <- c(
    "PhosphoPro_siteA_UP", "PhosphoPro_siteA_DOWN", "PhosphoPro_siteA_Non_significant", "PhosphoPro_siteA_Not_Detected",
    "PhosphoPro_siteB_UP", "PhosphoPro_siteB_DOWN", "PhosphoPro_siteB_Non_significant", "PhosphoPro_siteB_Not_Detected",
    "PhosphoPro_siteC_UP", "PhosphoPro_siteC_DOWN", "PhosphoPro_siteC_Non_significant", "PhosphoPro_siteC_Not_Detected",
    "PhosphoPro_siteD_UP", "PhosphoPro_siteD_DOWN", "PhosphoPro_siteD_Non_significant", "PhosphoPro_siteD_Not_Detected",
    "PhosphoPro_siteE_UP", "PhosphoPro_siteE_DOWN", "PhosphoPro_siteE_Non_significant", "PhosphoPro_siteE_Not_Detected",
    "UP", "DOWN", "Non_significant", "Not_Detected"
  )
  # Create a color mapping vector
  # fill_values <- c(
  #   "UP" = "#B83A2D",
  #   "DOWN" = "#657f68",
  #   "Non_significant" = "#E3C79F",
  #   "Not_Detected" = "grey90"
  # )
  fill_values <- c(
    "UP" = opt$node_up,    #"#B83A2D",
    "DOWN" = opt$node_down,    #"#657f68",
    "Non_significant" = opt$node_nonsig,    #"#E3C79F",
    "Not_Detected" = opt$node_notdet   #"grey90"
  )
  #  Add color mappings for all sites
  # for(site in c("E", "D", "C", "B", "A")) {
  #   fill_values[paste0("PhosphoPro_site", site, "_UP")] = "#B83A2D"
  #   fill_values[paste0("PhosphoPro_site", site, "_DOWN")] = "#657f68"
  #   fill_values[paste0("PhosphoPro_site", site, "_Non_significant")] = "#E3C79F"
  #   fill_values[paste0("PhosphoPro_site", site, "_Not_Detected")] = "grey90"
  # }
  for(site in c("E", "D", "C", "B", "A")) {
    fill_values[paste0("PhosphoPro_site", site, "_UP")] =  opt$node_up  #"#B83A2D"
    fill_values[paste0("PhosphoPro_site", site, "_DOWN")] =  opt$node_down    #"#657f68"
    fill_values[paste0("PhosphoPro_site", site, "_Non_significant")] = opt$node_nonsig    #"#E3C79F"
    fill_values[paste0("PhosphoPro_site", site, "_Not_Detected")] = opt$node_notdet    #"grey90"
  }
  # Extend the color mapping to support more types
  # edge_type_colors <- c(
  #   "association" = "#d7d6d6",          # 浅灰
  #   "physical association" = "#838181", # 浅灰（与 binding 相同）
  #   "binding" = "#838181",              # 浅灰
  #   "direct interaction" = "#d7d6d6",   # 浅灰
  # 
  #   "activation" = "#d62c0c",           # 橙红色
  #   "catalysis" = "#7f00ff",            # 淡紫色
  #   "proximity" = "#FF9301",            # 橙色
  #   "reaction" = "#114335",             # 深绿色
  #   "phosphorylation" = "#2FBE95",      # 蓝绿色
  #   "dephosphorylation" = "#6B8E23",    # 橄榄绿
  # 
  #   "ptmod" = "#8C97D6",                # 淡紫色
  #   "inhibition" = "#0cb6d6",           # 深蓝色
  #   "expression" = "#FCF402",           # 黄色
  #   "regulation" = "#e4dadd",           # 紫灰色（新增）
  # 
  #   "colocalization" = "#4c95cd",       # 天蓝色
  #   "covalent binding" = "#716F74",      # 紫罗兰色
  #   "ubiquitination" = "#FF4500",       # 橙红色
  #   "multiRel" = '#9a6728'             # 棕色
  # )
  edge_type_colors <- opt$edge_type_colors

  #  Define a line type mapping
  edge_type_linetypes <- c(
    "association" = "dotted",
    "physical association" = "dotted",
    "binding" = "dotted",
    "direct interaction" = "dotted",
    "activation" = "solid",
    "catalysis" = "solid",
    "proximity" = "solid",
    "reaction" = "solid",
    "phosphorylation" = "solid",
    "ptmod" = "solid",
    "inhibition" = "solid",
    "expression" = "solid",
    "colocalization" = "solid",
    "covalent binding" = "solid",
    "dephosphorylation" = "solid",
    "ubiquitination" = "solid",
    "multiRel" = "solid",
    "regulation" = "solid"
  )
  
  # Calculate shortening ratio
  shorten_ratio <- 0.1*length(rownames(nodes)) /10
  loop_offset = 0.1*length(rownames(nodes))/12  # 自环偏移量（根据节点大小调整）
  if(layout_type == "gem"){
    shorten_ratio <- 0.8*length(rownames(nodes)) * 2
    loop_offset <- 0.1*length(rownames(nodes))/(top_nodes_visualization_num+5)*300
  }else if(layout_type == "graphopt"){
    shorten_ratio <- 0.2*length(rownames(nodes))/(top_nodes_visualization_num)*40
    loop_offset <- 0.1*length(rownames(nodes))/(top_nodes_visualization_num+5)*120
  }else if(layout_type == "tree"){
    shorten_ratio <- 0.1*length(rownames(nodes))/5
  }else if(layout_type == "circle"){
    shorten_ratio <- 0.1*length(rownames(nodes))/(top_nodes_visualization_num+1)
    loop_offset <- 0.1*length(rownames(nodes))/(top_nodes_visualization_num+5)
  }else if(layout_type == "dh"){
    shorten_ratio <- 0.1*length(rownames(nodes))/(top_nodes_visualization_num)*10
    # loop_offset <- 0.1*length(rownames(nodes))/(top_nodes_visualization_num+5)
  }else if(layout_type == "lgl"){
    shorten_ratio <- 0.1*length(rownames(nodes))/(top_nodes_visualization_num)*8
    loop_offset <- 0.1*length(rownames(nodes))/(top_nodes_visualization_num)*6
  }

  # Add self-loop edge identifier
  edges <- edges %>%
    mutate(is_self_loop = (from == to))
  # print(head(edges))
  # print(length(rownames(nodes)))
  # Calculate adjusted edge coordinates (to fix self-loop edge issues)
  edges <- edges %>%
    mutate(
      dx = x_end -x_start,
      dy = y_end -y_start,
      length =sqrt(dx^2 + dy^2),  # Calculate the length of the edge
      # Adjusted x and y coordinates for non self-loops
      x_start_adj = ifelse(!is_self_loop, x_start + dx * shorten_ratio / length, x_start),
      y_start_adj = ifelse(!is_self_loop, y_start + dy * shorten_ratio / length, y_start),
      x_end_adj = ifelse(!is_self_loop, x_end - dx * shorten_ratio / length, x_end),
      y_end_adj = ifelse(!is_self_loop, y_end - dy * shorten_ratio / length, y_end)
    ) %>%
    # Adjusted x and y coordinates for self-loops
    mutate(
      # loop_offset = 0.1*length(rownames(nodes))/12,  # 自环偏移量（根据节点大小调整）
      x_start_adj = ifelse(is_self_loop, x_start - loop_offset, x_start_adj),
      y_start_adj = ifelse(is_self_loop, y_start + loop_offset*0.5, y_start_adj),
      x_end_adj = ifelse(is_self_loop, x_start + loop_offset, x_end_adj),
      y_end_adj = ifelse(is_self_loop, y_start + loop_offset*0.5, y_end_adj)
    )

  # Separate self-loop edges from regular edges
  self_loops_all <- edges %>% filter(is_self_loop)
  non_self_edges_all <- edges %>% filter(!is_self_loop)

  # Add an epsilon value (a very small non-zero value) to ensure all states are visible
  epsilon <- 0.00001 # Small enough to avoid significantly affecting the scale
  nodes$UP <- nodes$UP + epsilon
  nodes$DOWN <- nodes$DOWN + epsilon
  nodes$Non_significant  <- nodes$Non_significant  +epsilon
  nodes$Not_Detected  <- nodes$Not_Detected  +epsilon
  # print(table(non_self_edges$Interaction_Type_simplified))
  
  if(network_mode == 'ALL'){
  for(network_mode in c('ALL','LR','TF','KS')){
  # network_mode = 'LR'

  nettype <- strsplit(opt$outdir, "/")[[1]]
  nettype_name <- strsplit(nettype[length(nettype) - 1],'_')[[1]][1]
  if(nettype_name == 'ALL'){
    nettype_name <- ''
  }
  if(network_mode == 'KS'){
    title_str <- paste0("Kinase-Substrate Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
  }else if(network_mode == 'LR'){
    title_str <- paste0("Ligand-Receptor Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
  }else if(network_mode == 'TF'){
    title_str <- paste0("Transcriptional Regulation Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
  }else if(network_mode == 'DEP'){
    title_str <- paste0("Differential Protein and Phosphosite Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
  }else if(network_mode == 'DEPUP'){
    title_str <- paste0("Upregulated Protein and Phosphosite Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
  }else if(network_mode == 'DEPDOWN'){
  title_str <- paste0("Downregulated Protein and Phosphosite Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
  }else{
    title_str <- paste0("Comprehensive Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
  }
  
  
  non_self_edges <- non_self_edges_all
  self_loops <- self_loops_all
  
  if(network_mode == 'KS'){
    
    mode_edges <- E(g)[grepl("phospho", E(g)$Interaction_Type)]
    sub_g <- subgraph_from_edges(
      graph = g,
      eids = mode_edges,  # Preserved edges
      delete.vertices = TRUE  # Automatically delete disconnected nodes
    )
    # non_self_edges <- non_self_edges_all[(grepl("phospho", non_self_edges_all$Interaction_Type)),]
    # self_loops <- self_loops_all[(grepl("phospho", self_loops_all$Interaction_Type)),]
    non_self_edges$alpha_mask <- ifelse(grepl(pattern = "phospho",x = non_self_edges_all$Interaction_Type,ignore.case = TRUE), 0.5,0)
    self_loops$alpha_mask <-  ifelse(grepl(pattern = "phospho",x = self_loops_all$Interaction_Type,ignore.case = TRUE), 0.5,0)
    
    drop_col <- c("alpha_mask","color_mask")
    filter_edges <- rbind(non_self_edges[non_self_edges$alpha_mask != 0,],self_loops[self_loops$alpha_mask != 0,])
    # write_tsv(nodes[nodes$name %in% V(sub_g)$name,][,!(colnames(nodes) %in% drop_col)], file = file.path(opt$outdir,paste0(SEEDNODEID,"_nodes_",network_mode,".txt")))
    # write_tsv(filter_edges[,!(colnames(filter_edges) %in% drop_col)], file = file.path(opt$outdir,paste0(SEEDNODEID,"_edges_",network_mode,".txt")))
    write_tsv(nodes[nodes$name %in% V(sub_g)$name,][,!(colnames(nodes) %in% drop_col)], file = file.path(opt$outdir,paste0("nodes_",network_mode,".tsv")))
    write_tsv(filter_edges[,!(colnames(filter_edges) %in% drop_col)], file = file.path(opt$outdir,paste0("edges_",network_mode,".tsv")))
    
    
  }else if(network_mode == 'LR'){
    mode_edges <- E(g)[grepl("binding", E(g)$Interaction_Type)]
    sub_g <- subgraph_from_edges(
      graph = g,
      eids = mode_edges,
      delete.vertices = TRUE
    )
    # non_self_edges <- non_self_edges_all[grepl(pattern = "binding",x = non_self_edges_all$Interaction_Type,ignore.case = TRUE),]
    # self_loops <- self_loops_all[grepl(pattern = "binding",x = self_loops_all$Interaction_Type,ignore.case = TRUE),]
    non_self_edges$alpha_mask <- ifelse(grepl(pattern = "binding",x = non_self_edges_all$Interaction_Type,ignore.case = TRUE), 0.5,0)
    self_loops$alpha_mask <-  ifelse(grepl(pattern = "binding",x = self_loops_all$Interaction_Type,ignore.case = TRUE), 0.5,0)
    
    non_self_edges$Interaction_Type_simplified <- ifelse(grepl(pattern = "binding",x = non_self_edges_all$Interaction_Type),
                                                         ifelse(non_self_edges$Interaction_Type_simplified == 'binding', 'binding','multiRel'),
                                                         non_self_edges$Interaction_Type_simplified)
    self_loops$Interaction_Type_simplified <- ifelse(grepl(pattern = "binding",x = self_loops$Interaction_Type),
                                                     ifelse(self_loops$Interaction_Type_simplified == 'binding', 'binding','multiRel'),
                                                     self_loops$Interaction_Type_simplified)
    drop_col <- c("alpha_mask","color_mask")
    filter_edges <- rbind(non_self_edges[non_self_edges$alpha_mask != 0,],self_loops[self_loops$alpha_mask != 0,])
    # write_tsv(nodes[nodes$name %in% V(sub_g)$name,][,!(colnames(nodes) %in% drop_col)], file = file.path(opt$outdir,paste0(SEEDNODEID,"_nodes_",network_mode,".txt")))
    # write_tsv(filter_edges[,!(colnames(filter_edges) %in% drop_col)], file = file.path(opt$outdir,paste0(SEEDNODEID,"_edges_",network_mode,".txt")))
    write_tsv(nodes[nodes$name %in% V(sub_g)$name,][,!(colnames(nodes) %in% drop_col)], file = file.path(opt$outdir,paste0("nodes_",network_mode,".tsv")))
    write_tsv(filter_edges[,!(colnames(filter_edges) %in% drop_col)], file = file.path(opt$outdir,paste0("edges_",network_mode,".tsv")))
    
    
  }else if(network_mode == 'TF'){
    mode_edges <- E(g)[
      grepl(
        pattern = "activation|regulation|repression|inhibition",
        x = E(g)$Interaction_Type,
        ignore.case = TRUE  # 忽略大小写（可选）
      )
    ]
    sub_g <- subgraph_from_edges(
      graph = g,
      eids = mode_edges,
      delete.vertices = TRUE
    )
    non_self_edges$alpha_mask <- ifelse(grepl(pattern = "activation|regulation|repression|inhibition",x = non_self_edges_all$Interaction_Type,ignore.case = TRUE), 0.5,0)
    self_loops$alpha_mask <-  ifelse(grepl(pattern = "activation|regulation|repression|inhibition",x = self_loops_all$Interaction_Type,ignore.case = TRUE), 0.5,0)
    
    drop_col <- c("alpha_mask","color_mask")
    filter_edges <- rbind(non_self_edges[non_self_edges$alpha_mask != 0,],self_loops[self_loops$alpha_mask != 0,])
    # write_tsv(nodes[nodes$name %in% V(sub_g)$name,][,!(colnames(nodes) %in% drop_col)], file = file.path(opt$outdir,paste0(SEEDNODEID,"_nodes_",network_mode,".txt")))
    # write_tsv(filter_edges[,!(colnames(filter_edges) %in% drop_col)], file = file.path(opt$outdir,paste0(SEEDNODEID,"_edges_",network_mode,".txt")))
    write_tsv(nodes[nodes$name %in% V(sub_g)$name,][,!(colnames(nodes) %in% drop_col)], file = file.path(opt$outdir,paste0("nodes_",network_mode,".tsv")))
    write_tsv(filter_edges[,!(colnames(filter_edges) %in% drop_col)], file = file.path(opt$outdir,paste0("edges_",network_mode,".tsv")))
    
  }else{
    non_self_edges <- non_self_edges_all
    self_loops <-  self_loops_all
    sub_g <- g
    # write_tsv(nodes, file = file.path(opt$outdir,paste0(SEEDNODEID,"_nodes_",network_mode,".txt")))
    # write_tsv(edges, file = file.path(opt$outdir,paste0(SEEDNODEID,"_edges_",network_mode,".txt")))
    write_tsv(nodes, file = file.path(opt$outdir,paste0("nodes_",network_mode,".tsv")))
    write_tsv(edges, file = file.path(opt$outdir,paste0("edges_",network_mode,".tsv")))
  }
  
  mode_nodes <- V(sub_g)$name
  
  if(length(mode_nodes) > 0){
    # nodes <- nodes[nodes$name %in% mode_nodes,]
    # non_self_edges <- non_self_edges[(non_self_edges$from %in% nodes$name & non_self_edges$to %in% nodes$name), ]
    # self_loops <- self_loops[self_loops$from %in% nodes$name, ]
    
    nodes$alpha_mask <- rep(0,nrow(nodes))
    nodes$alpha_mask <- ifelse(nodes$name %in% mode_nodes,1,0)
    nodes$color_mask <- rep(NA,nrow(nodes))
    nodes$color_mask <- ifelse(nodes$name %in% mode_nodes,nodes$node_border_color,NA)
    # print(head(nodes))
    

  if(network_mode == 'ALL' ){
    sub_p <- ggplot() +
      geom_segment(
        data = non_self_edges,
        aes(x = x_start_adj,
            y = y_start_adj,
            xend = x_end_adj,
            yend = y_end_adj,
            color = Interaction_Type_simplified,  # 关键修改：映射边类型到颜色
            linetype = Interaction_Type_simplified
        ),
        alpha = 0.5,
        linewidth = 0.9,
        # color = "#2FBE95",
        arrow = arrow(  # 添加箭头
          type = "closed",  # 封闭箭头
          length = unit(0.2, "cm"),  # 箭头长度
          angle = 30,  # 箭头角度
        ),
        # show.legend = FALSE
      ) +
      # 2. 添加自环边（使用贝塞尔曲线）
      geom_curve(
        data = self_loops,
        aes(
          x = x_start_adj,  # 起点左偏移
          y = y_start_adj ,  # 起点上偏移
          xend = x_end_adj ,  # 终点右偏移
          yend = y_end_adj,  # 终点上偏移
          color = Interaction_Type_simplified,
          linetype = Interaction_Type_simplified
        ),
        curvature = -1,  # 负值表示逆时针曲线
        angle = 90,        # 控制曲线平滑度
        alpha = 0.5,
        linewidth = 0.9,
        arrow = arrow(
          type = "closed",
          length = unit(0.2, "cm"),
          angle = 30
        ),
        # show.legend = FALSE
      ) +
      scale_color_manual(
        name = "Edge Type",  # 图例名称
        values = edge_type_colors,
        na.value = "#d7d6d6"  # 为缺失值指定颜色
      ) +
      scale_linetype_manual(
        name = "Edge Type",  # 图例名称（与颜色相同，这样会合并为一个图例）
        values = edge_type_linetypes
      ) +
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale("linetype") +
      # 添加饼图节点 - 使用aes映射边框属性
      scatterpie::geom_scatterpie(
        data = nodes,
        aes(x = x, y = y, group = name,
            colour  = node_border_color,
            linetype = node_border_linetype),
        cols = pie_cols,
        pie_scale = 1.1,
        lwd =0.85,
      ) +
      ggrepel::geom_text_repel(
        data = nodes,
        aes(x = x,
            y = y,
            label = name,
        ),
        size = 4,  # 标签字体大小
        color = "grey30",  # 标签颜色
        family = "sans",  # 字体类型
        fontface = "bold",
        box.padding = 1,  # 标签与点的最小间距
        max.overlaps = Inf,  # 允许无限重叠（根据实际情况调整）
        min.segment.length = 0,  # 总是显示连接线
        segment.color = "grey50",  # 连接线颜色
        segment.size = 0.5,  # 连接线粗细
        force = 2  # 排斥力强度
      ) +
      scale_size_identity() +
      scale_fill_manual(
        name = "Protein/PhosphoProtein Expression",
        values = fill_values,
        breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
        drop = FALSE  # 关键：强制显示所有breaks项，即使数据中不存在
      ) +
      scale_linetype_identity(
        name = "Node Border",  # 边框图例名称
        guide = "legend",
        breaks = c("longdash", "solid"),
        labels = c("Disease-related", "Disease-Notrelated")
      ) +
      scale_colour_identity(
        name = "Node Border",  # 边框图例名称
        guide = "legend",
        # breaks = c("#535c54", "#acabab"),
        breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
        labels = c("Disease-related", "Disease-Notrelated")
      ) +
      # 等比例坐标系
      coord_fixed() +
      # 清理主题
      theme_void() +
      labs(title = title_str) +
      # ggtitle(title_str)+
      theme(
      plot.title = element_text(
      size = 20,                  # 字体大小
      face = "bold",              # 粗体
      # color = "grey30",          # 字体颜色
      family = "sans",            # 字体类型
      hjust = 0,                # 水平居中 (0=左, 0.5=中, 1=右)
      vjust = 1,                  # 垂直位置调整
      margin = margin(b = 15, t = 10),  # 上下边距
      lineheight = 0.8  # 行高
    ),
    plot.title.position = "plot",  # 标题与绘图区域对齐
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
        legend.position = "right",
        text = element_text(family = "sans"),
        panel.background = element_blank(),  # 移除面板背景
        plot.background = element_rect(fill = "white", color = NA),  # 设置绘图区域背景为白色
        legend.background = element_rect(fill = "white", color = NA),  # 设置图例背景为白色
        legend.key = element_blank(),  # 移除图例键背景
      )
  }else if( network_mode == 'SEEDNODE'){
    sub_p <- ggplot() +
      geom_segment(
        data = non_self_edges,
        aes(x = x_start_adj,
            y = y_start_adj,
            xend = x_end_adj,
            yend = y_end_adj,
            color = Interaction_Type_simplified,  # 关键修改：映射边类型到颜色
            linetype = Interaction_Type_simplified
        ),
        alpha = 0.5,
        linewidth = 0.9,
        # color = "#2FBE95",
        arrow = arrow(  # 添加箭头
          type = "closed",  # 封闭箭头
          length = unit(0.2, "cm"),  # 箭头长度
          angle = 30,  # 箭头角度
        ),
        # show.legend = FALSE
      ) +
      # 2. 添加自环边（使用贝塞尔曲线）
      geom_curve(
        data = self_loops,
        aes(
          x = x_start_adj,  # 起点左偏移
          y = y_start_adj ,  # 起点上偏移
          xend = x_end_adj ,  # 终点右偏移
          yend = y_end_adj,  # 终点上偏移
          color = Interaction_Type_simplified,
          linetype = Interaction_Type_simplified
        ),
        curvature = -1,  # 负值表示逆时针曲线
        angle = 90,        # 控制曲线平滑度
        alpha = 0.5,
        linewidth = 0.9,
        arrow = arrow(
          type = "closed",
          length = unit(0.2, "cm"),
          angle = 30
        ),
        # show.legend = FALSE
      ) +
      scale_color_manual(
        name = "Edge Type",  # 图例名称
        values = edge_type_colors,
        na.value = "grey80"  # 为缺失值指定颜色
      ) +
      scale_linetype_manual(
        name = "Edge Type",  # 图例名称（与颜色相同，这样会合并为一个图例）
        values = edge_type_linetypes
      ) +
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale("linetype") +
      # 添加饼图节点 - 使用aes映射边框属性
      scatterpie::geom_scatterpie(
        data = nodes,
        aes(x = x, y = y, group = name,
            colour  = node_border_color,
            linetype = node_border_linetype),
        cols = pie_cols,
        pie_scale = 1.1,
        lwd =0.85,
      ) +
      ggrepel::geom_text_repel(
        data = nodes,
        aes(x = x,
            y = y,
            label = name,
            size = ifelse(name == SEEDNODEID, 6, 4)  # 种子节点字号更大
        ),
        # size = 4,  # 标签字体大小
        color = "grey30",  # 标签颜色
        family = "sans",  # 字体类型
        fontface = "bold",
        box.padding = 1,  # 标签与点的最小间距
        max.overlaps = Inf,  # 允许无限重叠（根据实际情况调整）
        min.segment.length = 0,  # 总是显示连接线
        segment.color = "grey50",  # 连接线颜色
        segment.size = 0.5,  # 连接线粗细
        force = 2  # 排斥力强度
      ) +
      scale_size_identity() +
      scale_fill_manual(
        name = "Protein/PhosphoProtein Expression",
        values = fill_values,
        breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
        drop = FALSE  # 关键：强制显示所有breaks项，即使数据中不存在
      ) +
      scale_linetype_identity(
        name = "Node Border",  # 边框图例名称
        guide = "legend",
        breaks = c("longdash", "solid"),
        labels = c("Disease-related", "Disease-Notrelated")
      ) +
      scale_colour_identity(
        name = "Node Border",  # 边框图例名称
        guide = "legend",
        # breaks = c("#535c54", "#808080"),
        breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
        labels = c("Disease-related", "Disease-Notrelated")
      ) +
      # 等比例坐标系
      coord_fixed() +
      # 清理主题
      theme_void() +
      labs(title = title_str) +
      # ggtitle(title_str)+
      theme(
      plot.title = element_text(
      size = 20,                  # 字体大小
      face = "bold",              # 粗体
      # color = "grey30",          # 字体颜色
      family = "sans",            # 字体类型
      hjust = 0,                # 水平居中 (0=左, 0.5=中, 1=右)
      vjust = 1,                  # 垂直位置调整
      margin = margin(b = 15, t = 10),  # 上下边距
      lineheight = 0.8  # 行高
    ),
    plot.title.position = "plot",  # 标题与绘图区域对齐
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
        legend.position = "right",
        text = element_text(family = "sans"),
        panel.background = element_blank(),  # 移除面板背景
        plot.background = element_rect(fill = "white", color = NA),  # 设置绘图区域背景为白色
        legend.background = element_rect(fill = "white", color = NA),  # 设置图例背景为白色
        legend.key = element_blank(),  # 移除图例键背景
      )
  }else{
    sub_p <- ggplot() +
      geom_segment(
        data = non_self_edges,
        aes(x = x_start_adj,
            y = y_start_adj,
            xend = x_end_adj,
            yend = y_end_adj,
            color = Interaction_Type_simplified,  # 关键修改：映射边类型到颜色
            alpha = alpha_mask
            ),
        # alpha = 0.5,
        linewidth = 0.9,
        # color = "#2FBE95",
        arrow = arrow(  # 添加箭头
          type = "closed",  # 封闭箭头
          length = unit(0.2, "cm"),  # 箭头长度
          angle = 30,  # 箭头角度

        ),
        # show.legend = FALSE
      ) +
      # 2. 添加自环边（使用贝塞尔曲线）
      geom_curve(
        data = self_loops,
        aes(
          x = x_start_adj,  # 起点左偏移
          y = y_start_adj ,  # 起点上偏移
          xend = x_end_adj ,  # 终点右偏移
          yend = y_end_adj,  # 终点上偏移
          color = Interaction_Type_simplified,
          alpha = alpha_mask
        ),
        curvature = -1,  # 负值表示逆时针曲线
        angle = 90,        # 控制曲线平滑度
        # alpha = 0.5,
        linewidth = 0.9,
        arrow = arrow(
          type = "closed",
          length = unit(0.2, "cm"),
          angle = 30
        ),
        # show.legend = FALSE
      ) +
      scale_color_manual(
        name = "Edge Type",  # 图例名称
        values = edge_type_colors,
        na.value = "grey80"  # 为缺失值指定颜色
      ) +
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale("linetype") +
      # 添加饼图节点 - 使用aes映射边框属性
      scatterpie::geom_scatterpie(
        data = nodes,
        aes(x = x, y = y, group = name,
            # colour  = node_border_color,
            colour  = color_mask,
            linetype = node_border_linetype,
            alpha = alpha_mask),
        cols = pie_cols,
        pie_scale = 1.1,
        lwd =0.85,
      ) +
      ggrepel::geom_text_repel(
        data = nodes,
        aes(x = x, y = y, label = name,
            alpha = alpha_mask),
        size = 4,  # 标签字体大小
        color = "grey30",  # 标签颜色
        family = "sans",  # 字体类型
        fontface = "bold",
        box.padding = 1,  # 标签与点的最小间距
        max.overlaps = Inf,  # 允许无限重叠（根据实际情况调整）
        min.segment.length = 0,  # 总是显示连接线
        segment.color = "grey50",  # 连接线颜色
        segment.size = 0.5,  # 连接线粗细
        force = 2  # 排斥力强度
      ) +
      scale_fill_manual(
        name = "Protein/PhosphoProtein Expression",
        values = fill_values,
        breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
        drop = FALSE  # 关键：强制显示所有breaks项，即使数据中不存在
      ) +
      scale_linetype_identity(
        name = "Node Border",  # 边框图例名称
        guide = "legend",
        breaks = c("longdash", "solid"),
        labels = c("Disease-related", "Disease-Notrelated")
      ) +
      scale_colour_identity(
        name = "Node Border",  # 边框图例名称
        guide = "legend",
        # breaks = c("#535c54", "#acabab"),
        breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
        labels = c("Disease-related", "Disease-Notrelated")
      ) +
      # 等比例坐标系
      coord_fixed() +
      # 清理主题
      theme_void() +
      labs(title = title_str) +
      # ggtitle(title_str)+
      theme(
      plot.title = element_text(
      size = 20,                  # 字体大小
      face = "bold",              # 粗体
      # color = "grey30",          # 字体颜色
      family = "sans",            # 字体类型
      hjust = 0,                # 水平居中 (0=左, 0.5=中, 1=右)
      vjust = 1,                  # 垂直位置调整
      margin = margin(b = 15, t = 10),  # 上下边距
      lineheight = 0.8  # 行高
    ),
    plot.title.position = "plot",  # 标题与绘图区域对齐
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
        legend.position = "right",
        text = element_text(family = "sans"),
        panel.background = element_blank(),  # 移除面板背景
        plot.background = element_rect(fill = "white", color = NA),  # 设置绘图区域背景为白色
        legend.background = element_rect(fill = "white", color = NA),  # 设置图例背景为白色
        legend.key = element_blank(),  # 移除图例键背景
      )+
      guides(alpha = "none")  # 隐藏 alpha 图例
  }
  
  
  if(layout_type == 'tree')  {
    ggsave(filename = paste0("pro_phosphopro_network_",network_mode,".png"), plot = sub_p, width = 16, height = 12, units = "in", dpi = 300)
  }else{
    ggsave(filename = paste0("pro_phosphopro_network_",network_mode,".png"), plot = sub_p, width = 15, height = 12, units = "in", dpi = 300)
  }
  }
  
}
  }else{
    non_self_edges <- non_self_edges_all
    self_loops <-  self_loops_all
    # write_tsv(nodes, file = file.path(opt$outdir,paste0(SEEDNODEID,"_nodes_",network_mode,".txt")))
    # write_tsv(edges, file = file.path(opt$outdir,paste0(SEEDNODEID,"_edges_",network_mode,".txt")))
    write_tsv(nodes, file = file.path(opt$outdir,paste0("nodes_",network_mode,".tsv")))
    write_tsv(edges, file = file.path(opt$outdir,paste0("edges_",network_mode,".tsv")))
    
    nettype <- strsplit(opt$outdir, "/")[[1]]
    nettype_name <- strsplit(nettype[length(nettype) - 1],'_')[[1]][1]
    if(nettype_name == 'ALL'){
      nettype_name <- ''
    }
    if(network_mode == 'KS'){
      title_str <- paste0("Kinase-Substrate Visualization Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
    }else if(network_mode == 'LR'){
      title_str <- paste0("Ligand-Receptor Visualization Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
    }else if(network_mode == 'TF'){
      title_str <- paste0("Transcriptional Regulation Visualization Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
    }else if(network_mode == 'DEP'){
      title_str <- paste0("Differential Protein and Phosphosite Visualization Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
    }else if(network_mode == 'DEPUP'){
      title_str <- paste0("Upregulated Protein and Phosphosite Visualization Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
    }else if(network_mode == 'DEPDOWN'){
      title_str <- paste0("Downregulated Protein and Phosphosite Visualization Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
    }else{
      title_str <- paste0("Comprehensive Visualization Subnetwork of Proteomics and Phosphoproteomics ",nettype_name," Networks")
    }
    
    
    if( network_mode == 'SEEDNODE'){
      sub_p <- ggplot() +
        geom_segment(
          data = non_self_edges,
          aes(x = x_start_adj,
              y = y_start_adj,
              xend = x_end_adj,
              yend = y_end_adj,
              color = Interaction_Type_simplified,  # 关键修改：映射边类型到颜色
              linetype = Interaction_Type_simplified
          ),
          alpha = 0.5,
          linewidth = 0.9,
          # color = "#2FBE95",
          arrow = arrow(  # 添加箭头
            type = "closed",  # 封闭箭头
            length = unit(0.2, "cm"),  # 箭头长度
            angle = 30,  # 箭头角度
          ),
          # show.legend = FALSE
        ) +
        # 2. 添加自环边（使用贝塞尔曲线）
        geom_curve(
          data = self_loops,
          aes(
            x = x_start_adj,  # 起点左偏移
            y = y_start_adj ,  # 起点上偏移
            xend = x_end_adj ,  # 终点右偏移
            yend = y_end_adj,  # 终点上偏移
            color = Interaction_Type_simplified,
            linetype = Interaction_Type_simplified
          ),
          curvature = -1,  # 负值表示逆时针曲线
          angle = 90,        # 控制曲线平滑度
          alpha = 0.5,
          linewidth = 0.9,
          arrow = arrow(
            type = "closed",
            length = unit(0.2, "cm"),
            angle = 30
          ),
          # show.legend = FALSE
        ) +
        scale_color_manual(
          name = "Edge Type",  # 图例名称
          values = edge_type_colors,
          na.value = "grey80"  # 为缺失值指定颜色
        ) +
        scale_linetype_manual(
          name = "Edge Type",  # 图例名称（与颜色相同，这样会合并为一个图例）
          values = edge_type_linetypes
        ) +
        ggnewscale::new_scale_color() +
        ggnewscale::new_scale("linetype") +
        # 添加饼图节点 - 使用aes映射边框属性
        scatterpie::geom_scatterpie(
          data = nodes,
          aes(x = x, y = y, group = name,
              colour  = node_border_color,
              linetype = node_border_linetype),
          cols = pie_cols,
          pie_scale = 1.1,
          lwd =0.85,
        ) +
        ggrepel::geom_text_repel(
          data = nodes,
          aes(x = x,
              y = y,
              label = name,
              size = ifelse(name == SEEDNODEID, 6, 4)  # 种子节点字号更大
          ),
          # size = 4,  # 标签字体大小
          color = "grey30",  # 标签颜色
          family = "sans",  # 字体类型
          fontface = "bold",
          box.padding = 1,  # 标签与点的最小间距
          max.overlaps = Inf,  # 允许无限重叠（根据实际情况调整）
          min.segment.length = 0,  # 总是显示连接线
          segment.color = "grey50",  # 连接线颜色
          segment.size = 0.5,  # 连接线粗细
          force = 2  # 排斥力强度
        ) +
        scale_size_identity() +
        scale_fill_manual(
          name = "Protein/PhosphoProtein Expression",
          values = fill_values,
          breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
          drop = FALSE  # 关键：强制显示所有breaks项，即使数据中不存在
        ) +
        scale_linetype_identity(
          name = "Node Border",  # 边框图例名称
          guide = "legend",
          breaks = c("longdash", "solid"),
          labels = c("Disease-related", "Disease-Notrelated")
        ) +
        scale_colour_identity(
          name = "Node Border",  # 边框图例名称
          guide = "legend",
          # breaks = c("#535c54", "#808080"),
          breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
          labels = c("Disease-related", "Disease-Notrelated")
        ) +
        # 等比例坐标系
        coord_fixed() +
        # 清理主题
        theme_void() +
        labs(title = title_str) +
        # ggtitle(title_str)+
        theme(
          plot.title = element_text(
            size = 20,                  # 字体大小
            face = "bold",              # 粗体
            # color = "grey30",          # 字体颜色
            family = "sans",            # 字体类型
            hjust = 0,                # 水平居中 (0=左, 0.5=中, 1=右)
            vjust = 1,                  # 垂直位置调整
            margin = margin(b = 15, t = 10),  # 上下边距
            lineheight = 0.8  # 行高
          ),
          plot.title.position = "plot",  # 标题与绘图区域对齐
          plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
          legend.position = "right",
          text = element_text(family = "sans"),
          panel.background = element_blank(),  # 移除面板背景
          plot.background = element_rect(fill = "white", color = NA),  # 设置绘图区域背景为白色
          legend.background = element_rect(fill = "white", color = NA),  # 设置图例背景为白色
          legend.key = element_blank(),  # 移除图例键背景
        )
    }else{
      sub_p <- ggplot() +
        geom_segment(
          data = non_self_edges,
          aes(x = x_start_adj,
              y = y_start_adj,
              xend = x_end_adj,
              yend = y_end_adj,
              color = Interaction_Type_simplified,  # 关键修改：映射边类型到颜色
          ),
          alpha = 0.5,
          linewidth = 0.9,
          # color = "#2FBE95",
          arrow = arrow(  # 添加箭头
            type = "closed",  # 封闭箭头
            length = unit(0.2, "cm"),  # 箭头长度
            angle = 30,  # 箭头角度
            
          ),
          # show.legend = FALSE
        ) +
        # 2. 添加自环边（使用贝塞尔曲线）
        geom_curve(
          data = self_loops,
          aes(
            x = x_start_adj,  # 起点左偏移
            y = y_start_adj ,  # 起点上偏移
            xend = x_end_adj ,  # 终点右偏移
            yend = y_end_adj,  # 终点上偏移
            color = Interaction_Type_simplified,
          ),
          curvature = -1,  # 负值表示逆时针曲线
          angle = 90,        # 控制曲线平滑度
          alpha = 0.5,
          linewidth = 0.9,
          arrow = arrow(
            type = "closed",
            length = unit(0.2, "cm"),
            angle = 30
          ),
          # show.legend = FALSE
        ) +
        scale_color_manual(
          name = "Edge Type",  # 图例名称
          values = edge_type_colors,
          na.value = "grey80"  # 为缺失值指定颜色
        ) +
        ggnewscale::new_scale_color() +
        ggnewscale::new_scale("linetype") +
        # 添加饼图节点 - 使用aes映射边框属性
        scatterpie::geom_scatterpie(
          data = nodes,
          aes(x = x, y = y, group = name,
              colour  = node_border_color,
              linetype = node_border_linetype),
          cols = pie_cols,
          pie_scale = 1.1,
          lwd =0.85,
        ) +
        ggrepel::geom_text_repel(
          data = nodes,
          aes(x = x, y = y, label = name),
          size = 4,  # 标签字体大小
          color = "grey30",  # 标签颜色
          family = "sans",  # 字体类型
          fontface = "bold",
          box.padding = 1,  # 标签与点的最小间距
          max.overlaps = Inf,  # 允许无限重叠（根据实际情况调整）
          min.segment.length = 0,  # 总是显示连接线
          segment.color = "grey50",  # 连接线颜色
          segment.size = 0.5,  # 连接线粗细
          force = 2  # 排斥力强度
        ) +
        scale_fill_manual(
          name = "Protein/PhosphoProtein Expression",
          values = fill_values,
          breaks = c("UP", "DOWN", "Non_significant", "Not_Detected"),
          drop = FALSE  # 关键：强制显示所有breaks项，即使数据中不存在
        ) +
        scale_linetype_identity(
          name = "Node Border",  # 边框图例名称
          guide = "legend",
          breaks = c("longdash", "solid"),
          labels = c("Disease-related", "Disease-Notrelated")
        ) +
        scale_colour_identity(
          name = "Node Border",  # 边框图例名称
          guide = "legend",
          # breaks = c("#535c54", "#acabab"),
          breaks = c(opt$node_disease_border_color, opt$node_notdisease_border_color),
          labels = c("Disease-related", "Disease-Notrelated")
        ) +
        # 等比例坐标系
        coord_fixed() +
        # 清理主题
        theme_void() +
        labs(title = title_str) +
        # ggtitle(title_str)+
        theme(
          plot.title = element_text(
            size = 20,                  # 字体大小
            face = "bold",              # 粗体
            # color = "grey30",          # 字体颜色
            family = "sans",            # 字体类型
            hjust = 0,                # 水平居中 (0=左, 0.5=中, 1=右)
            vjust = 1,                  # 垂直位置调整
            margin = margin(b = 15, t = 10),  # 上下边距
            lineheight = 0.8  # 行高
          ),
          plot.title.position = "plot",  # 标题与绘图区域对齐
          plot.margin = margin(t = 30, r = 10, b = 10, l = 10, unit = "pt"),
          legend.position = "right",
          text = element_text(family = "sans"),
          panel.background = element_blank(),  # 移除面板背景
          plot.background = element_rect(fill = "white", color = NA),  # 设置绘图区域背景为白色
          legend.background = element_rect(fill = "white", color = NA),  # 设置图例背景为白色
          legend.key = element_blank(),  # 移除图例键背景
        )
    }
    if(layout_type == 'tree')  {
      ggsave(filename = paste0("pro_phosphopro_network_",network_mode,".png"), plot = sub_p, width = 16, height = 12, units = "in", dpi = 300)
    }else{
      ggsave(filename = paste0("pro_phosphopro_network_",network_mode,".png"), plot = sub_p, width = 15, height = 12, units = "in", dpi = 300)
    }
    
  }
  return("✅ Success!")
}










opt <- parameter_validation(opt)
######## topnodes_subgraph <- subnetwork_construction(opt,top_nodes_visualization_num)
######## network_visualiztion(topnodes_subgraph,opt)

# node_up = "#B83A2D"
# node_down = "#657f68"
# node_nonsig = "#E3C79F"
# node_notdet = "grey90"
# node_disease_border_color = "#535c54"
# node_notdisease_border_color = "#acabab"
# # 边类型颜色映射
# edge_type_colors <- c(
#   "association" = "#d7d6d6",          # 浅灰
#   "physical association" = "#838181", # 浅灰（与 binding 相同）
#   "binding" = "#838181",              # 浅灰
#   "direct interaction" = "#d7d6d6",   # 浅灰
#   
#   "activation" = "#d62c0c",           # 橙红色
#   "catalysis" = "#7f00ff",            # 淡紫色
#   "proximity" = "#FF9301",            # 橙色
#   "reaction" = "#114335",             # 深绿色
#   "phosphorylation" = "#2FBE95",      # 蓝绿色
#   "dephosphorylation" = "#6B8E23",    # 橄榄绿
#   
#   "ptmod" = "#8C97D6",                # 淡紫色
#   "inhibition" = "#0cb6d6",           # 深蓝色
#   "expression" = "#FCF402",           # 黄色
#   "regulation" = "#e4dadd",           # 紫灰色（新增）
#   
#   "colocalization" = "#4c95cd",       # 天蓝色
#   "covalent binding" = "#716F74",      # 紫罗兰色
#   "ubiquitination" = "#FF4500",       # 橙红色
#   "multiRel" = '#9a6728'             # 棕色
# )



topnodes_subgraph <- subnetwork_construction(opt,top_nodes_visualization_num)
if(is.character(topnodes_subgraph)){
  print(paste0("❌ No subnetworks can be ifilterred under ",network_mode_new," mode."))
}else if(length(V(topnodes_subgraph)$name)>2 ){
  network_visualiztion(topnodes_subgraph,opt)
}else{
  print(paste0("❌ No subnetworks can be ifilterred under ",network_mode_new," mode."))
}
# print("==========功能富集===================")
f_Node_rawdata <- read.csv(opt$nodes,sep = '\t', stringsAsFactors = FALSE)
pro_data_raw <- f_Node_rawdata[f_Node_rawdata$NodeName %in% V(topnodes_subgraph)$name,]
pro_list <- unique(pro_data_raw[pro_data_raw$Pro_class != "None",]$NodeName)
phos_list <- unique(pro_data_raw[pro_data_raw$Phosphopro_class != "None",]$NodeName)
print(file.path(opt$outdir,"Functional_Enrichment"))
omics_enrichment_list(pro_list,phos_list,file.path(opt$outdir,"Functional_Enrichment"),
                      omics1_name=opt$omics1_name,omics2_name=opt$omics2_name,
                      enrich_fromType = opt$enrich_fromType,
                      pvalueCutoff = 0.05,GO_showCategory=6,KEGG_showCategory=15,
                      color_gradient_low = opt$function_enrich_color_gradient_low,
                      color_gradient_high = opt$function_enrich_color_gradient_high)
# print("==========功能富集===================")





