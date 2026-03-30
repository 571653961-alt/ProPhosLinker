#' @title Differential Network Class
#'
#' @description An S4 class for storing differential network analysis results.
#' Contains comparison group names, differential nodes, and differential edge information.
#'
#' @slot group_name Character vector, comparison group name in the format "experimental_group-control_group"
#' @slot diff_nodes Data frame containing all differential node information
#' @slot diff_edges Data frame containing all differential edge information and their statistical metrics
#'
#' @details
#' Differential network analysis results include:
#' \itemize{
#' \item Comparison group name: Identifies the experimental conditions being compared
#' \item Differential nodes: Nodes that show differences between the two network groups
#' \item Differential edges: Edges that show differences between the two network groups, including:
#' \itemize{
#' \item Edges present in only one group
#' \item Edges present in both groups but with significant strength differences
#' \item Statistical metrics of edges (fold change, p-value, etc.)
#' }
#' }
#'
#' @seealso \code{\link{run_diff_network}} The function that generates this class
#' @export
setClass("Differential_network", slots = c(
  group_name="character",
  # Conditional_network_layout="ANY",
  diff_nodes = "data.frame",
  diff_edges = "data.frame"
))
#' @title Differential Network Analysis
#'
#' @description This function performs differential network analysis to compare network differences between two conditions, including changes in nodes and edges.
#'
#' @param Conditional_network Conditional network object containing case and control group networks
#' @param Conditional_multiplexnetwork Multiplex network object, can be used as an alternative input to Conditional_network
#' @param edge_FC_threshold Edge fold change threshold, default: 1.2
#' @param edge_p_threshold Edge significance p-value threshold, default: 0.05
#' @param compare_group Comparison groups in the format "experimental_group:control_group"
#'
#' @return Returns a Differential_network object containing:
#' \itemize{
#' \item group_name: Comparison group name
#' \item diff_nodes: Differential nodes data frame
#' \item diff_edges: Differential edges data frame
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#' \item Extract node and edge information for case and control groups from input objects
#' \item Parse comparison group names
#' \item Calculate edge differences using vectorized operations for high performance, identifying edges present in only one group and edges common to both groups
#' \item For edges common to both groups, use t-test to calculate significance of correlation differences
#' \item Filter significantly different edges based on thresholds
#' \item Merge node information
#' \item Return differential network object
#' }
#'
#' @examples
#' \dontrun{
#' # Perform differential analysis using conditional network object
#' # Create example network data
#' data("differential_network_result")
#' Conditional_network=differential_network_result@Conditional_network
#'
#' #Run differential network analysis
#' Differential_network<-run_diff_network(Conditional_network=Conditional_network,
#'                                       edge_FC_threshold=1.2,edge_p_threshold=0.05,compare_group="T:N")
#' # Visualization
#' network_show(Network = Differential_network,plot_type = 'diff_network',show_node_legend = TRUE,show_edge_legend = TRUE,node_colortype = 'Log2FC')
#'
#' # Perform differential analysis using multiplex network object
#' #' # Create example network data
#' data("multiplex_network_result")
#' Conditional_network=multiplex_network_result@Conditional_network
#' Conditional_multiplexnetwork=multiplex_network_result@Conditional_multiplexnetwork
#' Differential_multiplexnetwork<-run_diff_network(Conditional_network=Conditional_network,
#'                                                 Conditional_multiplexnetwork=Conditional_multiplexnetwork,
#'                                                 edge_FC_threshold=1.2,edge_p_threshold=0.05,compare_group="T:N")
#' #' # Visualization
#' network_show(Network = Differential_multiplexnetwork,plot_type = 'diff_network',show_node_legend = TRUE,show_edge_legend = TRUE,node_colortype = 'Log2FC')
#' }
#'
#' @export
run_diff_network<-function(Conditional_network=NULL,Conditional_multiplexnetwork=NULL,
                           edge_FC_threshold=1.2,edge_p_threshold=0.05,compare_group=NULL){
  if(any(!is.null(Conditional_multiplexnetwork))){
    case_node<-Conditional_multiplexnetwork@network_case@nodes
    case_edge<-Conditional_multiplexnetwork@network_case@edges
    control_node<-Conditional_multiplexnetwork@network_control@nodes
    control_edge<-Conditional_multiplexnetwork@network_control@edges
    # Conditional_network_layout<-NULL
  }else{
    case_node<-Conditional_network@network_case@bootnet_result_filter@bootnet_node
    case_edge<-Conditional_network@network_case@bootnet_result_filter@bootnet_edge
    control_node<-Conditional_network@network_control@bootnet_result_filter@bootnet_node
    control_edge<-Conditional_network@network_control@bootnet_result_filter@bootnet_edge
    # Conditional_network_layout<-Conditional_network@Conditional_network_layout
  }

  comparison<-gsub(":","-vs-", compare_group)
  split_result <- strsplit(comparison, "-vs-")[[1]]
  if(length(split_result)>1){
    casename <- split_result[1]
    controlname <- split_result[2]
  }else{
    casename <- "case"
    controlname <- "control"
  }
  Only_in_control<-paste0("Only in ",controlname)
  Only_in_case<-paste0("Only in ",casename)
  Enhanced_in_case<-paste0("Enhanced in ",casename)
  Enhanced_in_control<-paste0("Enhanced in ",controlname)

  # Process Case: Ensure Wide Format for fast calculation
  case_bootnet_stable_raw <- Conditional_network@network_case@bootnet_result_filter@bootnet_stable
  if (!any(grepl("^boot ", colnames(case_bootnet_stable_raw)))) {
      # Convert legacy Long format to Wide
      case_bootnet_stable <- case_bootnet_stable_raw |>
          tidyr::pivot_wider(id_cols = c(node1, node2), names_from = name, values_from = value)
  } else {
      case_bootnet_stable <- case_bootnet_stable_raw
  }
  # Rename for run_compare_edge
  colnames(case_bootnet_stable)[colnames(case_bootnet_stable) == "node1"] <- "from"
  colnames(case_bootnet_stable)[colnames(case_bootnet_stable) == "node2"] <- "to"
  # Add from_to ID
  case_bootnet_stable <- run_compare_edge(case_bootnet_stable)

  # Process Control: Ensure Wide Format
  control_bootnet_stable_raw <- Conditional_network@network_control@bootnet_result_filter@bootnet_stable
  if (!any(grepl("^boot ", colnames(control_bootnet_stable_raw)))) {
      control_bootnet_stable <- control_bootnet_stable_raw |>
          tidyr::pivot_wider(id_cols = c(node1, node2), names_from = name, values_from = value)
  } else {
      control_bootnet_stable <- control_bootnet_stable_raw
  }
  colnames(control_bootnet_stable)[colnames(control_bootnet_stable) == "node1"] <- "from"
  colnames(control_bootnet_stable)[colnames(control_bootnet_stable) == "node2"] <- "to"
  control_bootnet_stable <- run_compare_edge(control_bootnet_stable)

  ###edges diff
  cor_data_all1<-run_compare_edge(case_edge) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::rename(cor_case=cor) |>
    dplyr::rename(CIrange_case=CIrange) |>
    dplyr::rename(p_adjust_case=p_adjust)
  if("multiplex_status" %in% colnames(cor_data_all1)){
    cor_data_all1 <- cor_data_all1 |>
      dplyr::mutate(
        relationship = ifelse(
          grepl("_cor$", multiplex_status),
          sub("_cor$", "", multiplex_status),
          multiplex_status
        )
      ) |>
      dplyr::select(-multiplex_status)
  }

  cor_data_all2 <-run_compare_edge(control_edge) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::rename(cor_control=cor) |>
    dplyr::rename(CIrange_control=CIrange) |>
    dplyr::rename(p_adjust_control=p_adjust)
  if("multiplex_status" %in% colnames(cor_data_all2)){
    cor_data_all2 <- cor_data_all2 |>
      dplyr::mutate(
        relationship = ifelse(
          grepl("_cor$", multiplex_status),
          sub("_cor$", "", multiplex_status),
          multiplex_status
        )
      ) |>
      dplyr::select(-multiplex_status)
  }

  if("score" %in% colnames(cor_data_all1) && "score" %in% colnames(cor_data_all2)){
    cor_data_all <- dplyr::full_join(cor_data_all1, cor_data_all2, by = c("from", "to","score", "from_to","relationship")) |>
      dplyr::mutate(multiplex_status=dplyr::case_when(
        (!is.na(cor_case) | !is.na(cor_control)) & !is.na(score) ~ paste(relationship, "cor", sep = "_"),
        (!is.na(cor_case) | !is.na(cor_control)) & is.na(score) ~ "cor",
        (is.na(cor_case) & is.na(cor_control)) & !is.na(score) ~ relationship
      )) |>
      dplyr::select(-relationship)

  }else{
    cor_data_all <- dplyr::full_join(cor_data_all1, cor_data_all2, by = c("from", "to", "from_to"))
  }

  cor_data_all<-cor_data_all |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    as.data.frame() |>
    dplyr::mutate(
      case_control=cor_case*cor_control,
      cor_status = dplyr::case_when(
        (!is.na(case_control)) & (case_control<0) ~ "Conflict relation",
        (!is.na(case_control)) & (case_control>0) ~ "Both group",
        (is.na(cor_case)) & (!is.na(cor_control)) ~ Only_in_control,#"Only in control",
        (!is.na(cor_case)) & (is.na(cor_control)) ~ Only_in_case,#"Only in case",
        (is.na(cor_case)) & (is.na(cor_control)) ~ "Neither group"
      )) |>
    dplyr::select(-case_control)


  cor_data_both<-cor_data_all |>
    dplyr::filter(cor_status=="Both group")
  
  if(nrow(cor_data_both)>0){
    # Target edges
    target_edges <- cor_data_both$from_to
    
    # Filter and align data rows
    case_mat_sub <- case_bootnet_stable[case_bootnet_stable$from_to %in% target_edges, ]
    control_mat_sub <- control_bootnet_stable[control_bootnet_stable$from_to %in% target_edges, ]
    
    case_mat_sub <- case_mat_sub[match(target_edges, case_mat_sub$from_to), ]
    control_mat_sub <- control_mat_sub[match(target_edges, control_mat_sub$from_to), ]
    
    # Extract matrices and take absolute values (as per original logic)
    mat_case <- abs(as.matrix(case_mat_sub[, grepl("^boot ", colnames(case_mat_sub))]))
    mat_control <- abs(as.matrix(control_mat_sub[, grepl("^boot ", colnames(control_mat_sub))]))
    
    # 1. Calculate FC (Mean / Mean)
    mean_case <- rowMeans(mat_case, na.rm=TRUE)
    mean_control <- rowMeans(mat_control, na.rm=TRUE)
    fc_vec <- ifelse(mean_control == 0, Inf, mean_case / mean_control)
    
    # 2. Calculate P-value (Row-wise T-test on Log2)
    log_case <- log2(mat_case)
    log_control <- log2(mat_control)
    # Handle -Inf
    log_case[!is.finite(log_case)] <- NA
    log_control[!is.finite(log_control)] <- NA
    
    p_vec <- sapply(seq_len(nrow(log_case)), function(i) {
      tryCatch({
        t.test(log_case[i, ], log_control[i, ])$p.value
      }, error = function(e) NA)
    })
    
    results_df <- data.frame(cor_FC = fc_vec, cor_p_value = p_vec)
    
    results_df <- results_df |>
      dplyr::mutate(cor_status= dplyr::case_when(
        cor_p_value < edge_p_threshold & cor_FC > edge_FC_threshold ~ Enhanced_in_case,
        cor_p_value < edge_p_threshold & cor_FC < 1 / edge_FC_threshold ~ Enhanced_in_control,
        TRUE ~ "Non-significant"
      ))
    results_df$from_to <- target_edges

    cor_data_diff <- cor_data_all |>
      dplyr::left_join(results_df, by = "from_to")

    # Update cor_status in cor_data_diff where there is a match
    cor_data_diff$cor_status <- ifelse(!is.na(cor_data_diff$cor_status.y),
                                       cor_data_diff$cor_status.y,
                                       cor_data_diff$cor_status.x)

    # Remove the temporary columns used for joining
    cor_data_diff <- cor_data_diff |>
      dplyr::select(-cor_status.x, -cor_status.y)

  }else{
    cor_data_diff<- cor_data_all |>
      dplyr::mutate(cor_status=cor_status,
                    cor_FC = NA,
                    cor_p_value = NA)
  }
  cor_data_diff<- cor_data_diff |>
    dplyr::mutate(cor = cor_FC)

  ##nodes diff
  merge_node <- rbind(case_node,control_node) |>
    dplyr::distinct(node, .keep_all = TRUE) #|>

  Differential_network_result <- new("Differential_network",
                                     group_name=comparison,
                                     diff_nodes=merge_node,
                                     diff_edges=cor_data_diff
  )

  return(Differential_network_result)
}



#' @title Edge Comparison Preprocessing Function
#'
#' @description Preprocesses network edge data to generate standardized edge identifiers for comparison.
#' By sorting node pairs and creating unique edge identifiers, ensures that undirected edge comparisons
#' are not affected by node order.
#'
#' @param edges Data frame containing network edge information, must include "from" and "to" columns
#'
#' @return Returns processed edge data frame with an additional "from_to" column serving as unique edge identifier
#' (format: "larger_node_smaller_node"), with temporary sorting columns removed
#'
#' @details
#' This function performs the following operations:
#' \enumerate{
#' \item Sorts the two nodes of each edge to ensure uniform representation
#' \item Creates unique edge identifiers in the format "larger_node_smaller_node"
#' \item Removes temporary sorting columns to maintain data cleanliness
#' }
#'
#' @examples
#' \dontrun{
#' # Create sample edge data
#' edges <- data.frame(
#' from = c("A", "B", "C"),
#' to = c("B", "A", "D"),
#' cor = c(0.5, 0.8, 0.3)
#' )
#'
#' # Preprocess edge data
#' processed_edges <- run_compare_edge(edges)
#' }
#'
#' @keywords internal
run_compare_edge<-function(edges){
  edges <-transform(edges, sorted_from = pmin(from, to), sorted_to = pmax(from, to))
  edges$from_to<-paste0(edges$sorted_to,"_",edges$sorted_from)
  edges<-edges |>
    dplyr::select(-sorted_from, -sorted_to)
  return(edges)
}


