
#' Enrichment Analysis Class and Functions
#'
#' This module provides comprehensive enrichment analysis capabilities for
#' multi-omics data, supporting KEGG pathway enrichment across different
#' omics types (e.g., proteomics, phosphoproteomics). It includes:
#'
#'   - S4 class for storing enrichment results
#'   - Parallel processing support for multi-subnet analysis
#'   - KEGG database integration with local file support
#'   - Multi-omics comparison and result merging
#'   - P-value calculation using hypergeometric distribution
#'
#' The enrichment pipeline processes annotation tables containing KEGG IDs,
#' performs pathway enrichment analysis for specified omics types, and
#' returns structured results with pathway annotations, statistics, and
#' filtered results based on significance thresholds.
#'
#' @section Main Components:
#'   - Enrichment class: Stores group name and network enrichment results
#'   - run_enrichment(): Main function orchestrating the enrichment analysis
#'   - process_subnet(): Handles enrichment for individual subnetworks
#'   - process_annotation_table(): Annotates KEGG IDs with omics types
#'   - merge_enrich_tables(): Combines enrichment results from multiple omics
#'
#' @section Data Requirements:
#'   Input annotation tables must contain:
#'   - KEGG.ID column with pathway identifiers
#'   - node or feature_ID column for feature identification
#'   - Optional omics_name column for multi-omics separation
#'
#' @section KEGG Database Integration:
#'   The module supports local KEGG database files for offline analysis:
#'   - genes_ko.list: Gene to KO mapping
#'   - ko_map.tab: KO to pathway mapping
#'   - map_title.tab: Pathway descriptions
#'   - species-specific list files for target organisms
#'
#' @name EnrichmentAnalysis
#' @rdname EnrichmentAnalysis


setClass("Enrichment", slots = c(
  group_name="character",
  network = "ANY"
))

run_enrichment<-function(annotation_table_select=NULL,species=NULL,nCores=NULL,omics_name=NULL,group_name="data",enrichment_p_threshold=0.05,
database_path=NULL,glist_path=NULL){

  if(length(group_name)>1){
    group_name<-paste(group_name, collapse = "-vs-")
  }else{
    group_name<-gsub(":","-vs-", group_name)
  }

  e_omics_name<-NULL
  omics<-NULL
  if("KEGG.ID" %in% names(annotation_table_select)){
    node_table<-list(network=annotation_table_select)
  }else{
    node_table<-annotation_table_select 
  }

  node_table <- process_annotation_table(node_table) 


  all_have_omics_name <- all(sapply(node_table, function(df) "omics_name" %in% colnames(df)))
  if(all_have_omics_name){

    real_omics_name<-node_table |>
      purrr::map(~ .x$omics_name) |>
      unlist() |>
      as.character()
    real_omics_name<-unique(na.omit(real_omics_name))
  }else{
    real_omics_name<-NULL
  }
  if(!is.null(omics_name)){
    if (!all_have_omics_name) {
      stop("The annotation table does not contain an omics_name column, enrichment analysis terminated")
    }

    intersect_omics_name<-intersect(real_omics_name,omics_name)
    message(paste0("Detect the value of omics_name column from annotation_table_select: ", paste0(intersect_omics_name, collapse = ", ")))
    if(length(intersect_omics_name)==0){
      stop("Mismatch between the omics_name column of the annotation table and the omics_name parameter!")
    }
  }else{
    intersect_omics_name<-real_omics_name
  }


  omicslist<-node_table |>
    purrr::map(~ .x$omics_type) |>
    unlist() |>
    as.character()
  omicslist<-unique(na.omit(omicslist))
  path<-database_path#file.path(root_path,"Database") 
  if(all(is.na(omicslist))){
    message(paste0("Enrichment analysis of ",group_name," terminated due to missing KEGG IDs.") )
    #return(NULL)
    return(new("Enrichment",
               group_name=group_name,
               network =NULL
    ))
  }
closeAllConnections()
  database<-lapply(omicslist,function(omics){
   if(grepl("^Omics1",omics)){
     CPD_KO<-"cpd"
     ko_col<-2
     gff_data_filename<-"kegg_106.0_chart.txt"
     gff_data_grep<-"map"
   }else if(grepl("^Omics2",omics)){
     CPD_KO<-"ko"
     ko_col<-3
     gff_data_filename<-"ko_kegg_106.0_chart.txt"#第一列手动改成了KEGG_ID
     gff_data_grep<-"ko"
   }

   # df <- read_delim(file.path(path,"species.txt"),delim="\t") |> as.data.frame()
   sp<-species
   # sp<-"hsa"

   lines_file <-paste0(glist_path,"/Data_Analysis/Pro_network/species/", sp,"/",sp, ".list")

   # lines_file <- system.file("data", "komap", sp, paste0(sp, ".list"), package = "OmicsNetwork")
   lines <- readLines(lines_file)

   filtered_lines <- lines[grep(CPD_KO, lines, ignore.case = TRUE)]

   data <- strsplit(filtered_lines, "\t")

   result <- data.frame(CPD = character(), KO = character(), stringsAsFactors = FALSE)

   ko_dict <- list()

   for (line in data) {
     cpd <- line[1]
     ko <- line[ko_col]
     if (ko %in% names(ko_dict)) {
       ko_dict[[ko]] <- c(ko_dict[[ko]], cpd)
     } else {
       ko_dict[[ko]] <- cpd
     }
   }

   for (ko in names(ko_dict)) {
     result <- rbind(result, data.frame(CPD = ko, KO = paste(ko_dict[[ko]], collapse = "\t"), stringsAsFactors = FALSE))
   }

   result <- result[order(result$CPD), ]
   result <- unique(result)
   result$CPD <- gsub("cpd:", "",result$CPD)
  
   result$CPD <- gsub("ko:", "",result$CPD)
   result$CPD <- sub(" .*", "", result$CPD)
   result$KO <- gsub("\t", " ",result$KO)
   tab_data <- result

   #gff_data_file <- system.file("data", gff_data_filename, package = "OmicsNetwork")
   gff_data_file<-file.path(path,paste0("KEGG/",gff_data_filename))#"KEGG_106.txt"
   gff_data <- suppressMessages(readr::read_delim(gff_data_file, delim = "\t")) |> 
     as.data.frame()
 #  str(gff_data)
   if("ENTRY" %in% colnames(gff_data)){
    gff_data<- gff_data |>
         dplyr::rename("KEGG_ID"="ENTRY") |>
         dplyr::filter(!is.na(KEGG_ID)) |>
          dplyr::filter(KEGG_ID != "NULL")
   }

   CtoMap0 <- gff_data |> dplyr::mutate(Map = ifelse(grepl(gff_data_grep, PATHWAY), PATHWAY, NA))#"map"
   CtoMap <- setNames(CtoMap0$Map, CtoMap0$KEGG_ID)
   list(tab_data=tab_data, gff_data=gff_data, CtoMap=CtoMap)
 })
 names(database)<-omicslist
 #print(names(database[[1]]$CtoMap))
  ###################################
  #Define the function to process each subnet
  if(is.null(nCores)){
    nCores <- parallel::detectCores() - 1
  }else{
    nCores <- nCores#parallel::detectCores() - 1
  }
  if(length(intersect_omics_name)>0){ 
    e_omics_n=intersect_omics_name
  }else{
    e_omics_n=omicslist 
  }
 networks_Enrichment_list<-list()
 for (e_omics_name in e_omics_n) {
#lapply(e_omics_n, function(e_omics_name){
   cl <- parallel::makeCluster(nCores)
   parallel::clusterExport(cl, c("process_subnet","node_table", "database","e_omics_name","enrichment_p_threshold"), envir=environment())#指定不是全局数据而是当前函数里的环境数据
   networks_Enrichment_list[[e_omics_name]] <- parallel::parLapply(cl, names(node_table), function(subnet_num) {
     process_subnet(subnet_num=subnet_num, node_table=node_table,database=database,
                    e_omics_name=e_omics_name,enrichment_p_threshold=enrichment_p_threshold)
   })
   parallel::stopCluster(cl)
# })
 }
  if(length(e_omics_n)==1){
   # str(networks_Enrichment_list)
    networks_Enrichment<-networks_Enrichment_list[[1]]
    names(networks_Enrichment)<- names(node_table)
    networks_enrichment_result <- new("Enrichment",
                                      group_name=group_name,
                                      network =networks_Enrichment
    )
  }else{
    merge_subsets <- function(subsets) {
      max_subnets <- max(sapply(subsets, length))
      merged_nodes <- lapply(seq_len(max_subnets), function(subnet_index) {
        merge_enrich_tables(subsets,subnet_index, table_name="nodes")
      })
      merged_annotations <- lapply(seq_len(max_subnets), function(subnet_index) {
        merge_enrich_tables(subsets,subnet_index, table_name="annotations")
      })
      merged_annotations_filter <- lapply(seq_len(max_subnets), function(subnet_index) {
        annotations = merged_annotations[[subnet_index]]
        if (any(!is.na(annotations))) {
          plotdata <- annotations |> dplyr::mutate(RichFactor = Count / CountAll) |>
            dplyr::rename(Count = Count) |>
            dplyr::arrange(Pvalue) 
          plotdata$Count <- as.numeric(plotdata$Count)
          if(length(unique(plotdata$Types))==1){
            plotdata |>
            dplyr::filter(Pvalue < enrichment_p_threshold) 
          }else{
            plotdata |>
              dplyr::group_by(Pathway) |>
              dplyr::filter(dplyr::n() > 1) |>
              dplyr::ungroup()
          }
        } else {
          annotations
        }

      })
      
      network <- lapply(seq_along(merged_nodes), function(i) {
        list(
          nodes = merged_nodes[[i]],
          annotations = merged_annotations[[i]],
          annotations_filter = merged_annotations_filter[[i]]
        )
      })
      network
    }

    # 合并所有子集
    networks_Enrichment <- merge_subsets(subsets=networks_Enrichment_list)
    names(networks_Enrichment)<- names(node_table)
    networks_enrichment_result <- new("Enrichment",
                                      group_name=group_name,
                                      network =networks_Enrichment
    )
  }

  return(networks_enrichment_result)
}

####
process_subnet <- function(subnet_num=NULL, node_table=NULL, database=NULL,e_omics_name=NULL,enrichment_p_threshold=NULL) {
  subnet <- node_table[[subnet_num]] |>#data.fram(node,KEGG.ID)
  as.data.frame() |> 
  dplyr::filter(!is.na(KEGG.ID)) |>
  dplyr::filter(KEGG.ID!="NA")
  if("omics_name" %in% colnames(subnet)){
    subnet <-subnet |>
      dplyr::filter(omics_name==e_omics_name)
  }else{
    subnet <-subnet |>
      dplyr::filter(omics_type==e_omics_name)
  }
  if(nrow(subnet)>0){
 
  ###########################################################
  if(any(grepl("^Omics1",subnet$omics_type))){
    CK<-"^C"
    map_ko<-"map"
 #   list(tab_data=tab_data, gff_data=gff_data, CtoMap=CtoMap)
  }else if(any(grepl("^Omics2",subnet$omics_type))){
    CK<-"^K"
    map_ko<-"ko"
  }else{
list(nodes = NULL, annotations = NULL, annotations_filter = NULL)
}
 subnet<-subnet |>
  dplyr::filter(grepl(CK, KEGG.ID))
  sub_database_type<-unique(subnet$omics_type)
  sub_tab_data<-database[[sub_database_type]]
  tab_data=sub_tab_data$tab_data
  gff_data=sub_tab_data$gff_data
  CtoMap=sub_tab_data$CtoMap

 # Read selectCid file
  if(c("node") %in% colnames(subnet)){
    selectCid_data <- subnet |>
      dplyr::select(KEGG.ID, node)  
  }else if(c("feature_ID") %in% colnames(subnet)){
    selectCid_data <- subnet |>
      dplyr::select(KEGG.ID, feature_ID) |>
      dplyr::rename(node=feature_ID) 
  }else{
    stop("The column name of annotation_table_select must contain ‘node’ or ‘feature_ID’.")
  }


  # Initialize data structures
  c_ko <- list()
  path <- list()
  sum1 <- nrow(selectCid_data)
  sum2 <- nrow(gff_data)
  
  # Process tab file
  for (i in 1:nrow(tab_data)) {
    row <- tab_data[i, ]
    c_ko[[as.character(row$CPD)]] <- row$KO
  }
  
  # Process selectCid file
  for (i in 1:nrow(selectCid_data)) {
    row <- selectCid_data[i, ]
    if (grepl(CK, row$KEGG.ID)) {#"^C"
      if (length(grep(row$KEGG.ID, gff_data$KEGG_ID)) > 0) {
        if (grepl(map_ko, CtoMap[[row$KEGG.ID]])) {#"map"
          maps <- strsplit(CtoMap[[row$KEGG.ID]], ";")[[1]]
          for (map in maps) {
            map_match <- regmatches(map, regexec(paste0(map_ko,"(\\d+)\\s+(.*)"), map))[[1]]#"map(\\d+)\\s+(.*)"
            if (length(map_match) > 0) {
              mapId <- map_match[2]
              mapName <- map_match[3]
              if (exists(row$KEGG.ID, where = c_ko) && grepl(mapId, c_ko[[row$KEGG.ID]])) {
                path[[mapId]]$mapName <- mapName
                path[[mapId]]$cId <- unique(c(path[[mapId]]$cId, row[1]))
              }
            }
          }
        }
      }
    }
  }
  
  # Process all GFF data
  for (i in 1:nrow(gff_data)) {
    row <- gff_data[i, ]
    if (grepl(CK, row$KEGG_ID)) {#"^C"
      if (grepl(map_ko, row$PATHWAY)) {#"map"
        maps <- strsplit(row$PATHWAY, ";")[[1]]
        for (map in maps) {
          map_match <- regmatches(map, regexec(paste0(map_ko,"(\\d+)\\s+(.*)"), map))[[1]]#"map(\\d+)\\s+(.*)"
          if (length(map_match) > 0) {
            mapId <- map_match[2]
            if (exists(mapId, where = path)) {
              path[[mapId]]$cIdAll <- unique(c(path[[mapId]]$cIdAll, row$KEGG_ID))
            }
          }
        }
      }
    }
  }
  
  # Calculate P-values
  p_values <- numeric()
  for (pwID in sort(names(path))) {
    num1 <- length(path[[pwID]]$cId)
    num2 <- length(path[[pwID]]$cIdAll)
    p_value <- phyper(num1 - 1, num2, sum2 - num2, sum1, lower.tail = FALSE)
    p_values <- c(p_values, p_value)
  }
  
  # Initialize content variable
  content <- data.frame()
  
  # Add P-values to path list
  i <- 1
  for (pwID in sort(names(path))) {
    path[[pwID]]$pvalue <- p_values[i]
    i <- i + 1
  }
  
  # Sort by P-value
  sorted_path <- path[order(sapply(path, function(x) x$pvalue))]
  
  # Generate output content
  for (pwID in names(sorted_path)) {
    num1 <- length(path[[pwID]]$cId)
    num2 <- length(path[[pwID]]$cIdAll)
    mapName <- path[[pwID]]$mapName
    pvalue <- path[[pwID]]$pvalue
    cIds <- paste(path[[pwID]]$cId, collapse = "+")
    temp <- data.frame(Pathway = mapName, Count = num1, CountAll = num2, Pvalue = pvalue, PathwayID = paste0("map", pwID), `KEGG IDs` = cIds)
    content <- dplyr::bind_rows(content, temp)
  }
#  print(content)
  if (any(!is.na(content))) {
    content$Types<-e_omics_name
    plotdata <- content |> dplyr::mutate(RichFactor = Count / CountAll) |>
      dplyr::rename(Count = Count) |>
      dplyr::filter(Pvalue < enrichment_p_threshold) |>
      dplyr::arrange(Pvalue) 
    plotdata$Count <- as.numeric(plotdata$Count)
  } else {
    plotdata <- content
  }

  list(nodes = subnet, annotations = content, annotations_filter = plotdata)
  }else{

   list(nodes = NULL, annotations = NULL, annotations_filter = NULL)
  }
}
############
process_annotation_table <- function(annotation_table_select) {
  if (is.list(annotation_table_select)) {

    processed_omics_names <- annotation_table_select |>
      purrr::map(~ {
          .x |>
            dplyr::mutate(
              omics_type = dplyr::case_when(
                grepl("^C", KEGG.ID) ~ "Omics1",
                grepl("^K", KEGG.ID) ~ "Omics2",
                TRUE ~ NA_character_
              )) 
            

      }) 
  } else {
    stop("annotation_table_select must be a list")
  }
  
  return(processed_omics_names)
}

merge_enrich_tables <- function(subsets,subnet_index,table_name) {
  # Initialize an empty list to store non-NULL node_tables
  tables <- list()
  
  # Loop through each omics type in subsets
  for (omics_type in names(subsets)) {
    # Check if the subnet_index exists in the current omics type
    if (subnet_index <= length(subsets[[omics_type]])) {
      # Get the node_table for the current subnet_index
      node_table <- subsets[[omics_type]][[subnet_index]][[table_name]]
      
      # Check if the node_table is not NULL
      if (!is.null(node_table)) {
        tables[[omics_type]] <- node_table
      }
    }
  }
  
  # Merge all non-NULL tables
  if (length(tables) > 0) {
    # Combine all tables into a single data frame
    merged_table <- do.call(rbind, tables)
  } else {
    merged_table <- data.frame()  # Return an empty data frame if no non-NULL tables
  }
  if(table_name=="nodes" && "node" %in% colnames(merged_table)){
    merged_table<-merged_table |>
      dplyr::distinct(node, .keep_all = TRUE)
  }else if(table_name=="annotations_filter" && "Pathway" %in% colnames(merged_table)){
    
    merged_table<-merged_table |>
      dplyr::group_by(Pathway) |>
      dplyr::filter(dplyr::n() > 1) |>
      dplyr::ungroup()
  }
  
  return(merged_table)
}

