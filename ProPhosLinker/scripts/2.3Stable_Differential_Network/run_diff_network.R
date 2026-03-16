#' Perform differential network analysis between two conditions (e.g., case vs control).
#'
#' This function compares two conditional networks (either from `Conditional_network` or 
#' `Conditional_multiplexnetwork` objects) to identify edges that are:
#' - present only in one condition,
#' - shared but significantly different in strength (via fold-change and t-test on bootstrapped correlations),
#' - or conflicting in sign.
#' 
#' It also merges node information from both networks. The result is stored in a `Differential_network` S4 object.
#'
#' @param Conditional_network An object of class containing `network_case` and `network_control`,
#'        each with slots `bootnet_result_filter` holding `bootnet_node`, `bootnet_edge`, and `bootnet_stable`.
#'        Used when multiplex network is not provided.
#' @param Conditional_multiplexnetwork An alternative input object containing precomputed `network_case` and
#'        `network_control` with direct `nodes` and `edges` slots. If provided, it takes precedence.
#' @param edge_FC_threshold Fold-change threshold for absolute correlation difference between groups (default: 1.2).
#' @param edge_p_threshold P-value threshold for significance of correlation difference (default: 0.05).
#' @param compare_group Character string specifying the comparison label (e.g., "Case:Control").
#'        Colons (`:`) are automatically replaced with "-vs-".
#'
#' @return An object of S4 class `Differential_network` with slots:
#'   \describe{
#'     \item{group_name}{Formatted comparison name (e.g., "Case-vs-Control").}
#'     \item{diff_nodes}{Combined node data frame from both conditions.}
#'     \item{diff_edges}{Edge-level comparison results, including status ("Only in case", "Enhanced in control", etc.), 
#'                       fold-change (`cor_FC`), p-value, and original correlation values.}
#'   }
#'
#' @importFrom dplyr select rename distinct filter mutate case_when full_join left_join pull
#' @importFrom stats t.test pmin pmax
#' @importFrom base transform
#' @import dplyr
#'
#' @examples
#' # This function is typically used internally after constructing conditional networks.
#' # See package vignettes for full workflow.
#'
#' @export

setClass("Differential_network", slots = c(
  group_name="character",
  diff_nodes = "data.frame",
  diff_edges = "data.frame"
))

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
  split_result <- strsplit(comparison, "-vs-")
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
  
  case_bootnet_stable<-Conditional_network@network_case@bootnet_result_filter@bootnet_stable |>
    dplyr::select(node1,node2,value)
  colnames(case_bootnet_stable)<-c("from","to","cor")
  
  control_bootnet_stable<-Conditional_network@network_control@bootnet_result_filter@bootnet_stable |>
    dplyr::select(node1,node2,value)
  colnames(control_bootnet_stable)<-c("from","to","cor")

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
    cor_data_boot1 <-run_compare_edge(case_bootnet_stable) |>
      dplyr::filter(from_to  %in% cor_data_both$from_to) 
    cor_data_boot2 <-run_compare_edge(control_bootnet_stable) |>
      dplyr::filter(from_to  %in% cor_data_both$from_to)
    results<-lapply(cor_data_both$from_to,function(x){
      casedata<-cor_data_boot1 |>
        dplyr::filter(from_to==x) |>
        dplyr::pull(cor) |>
        abs()
      condata<-cor_data_boot2 |>
        dplyr::filter(from_to==x) |>
        dplyr::pull(cor) |>
        abs()
      ca_mean <- mean(casedata)
      co_mean <- mean(condata)
      fd <- ca_mean/co_mean 
      p <- t.test(log2(casedata),log2(condata))
      k <- c(fd,p$p.value)
      names(k) <- c("cor_FC","cor_p_value") 
      return(k)
    })
    results_df <- do.call(rbind, results)
    results_df<-as.data.frame(results_df) |>   
      dplyr::mutate(cor_status= dplyr::case_when(
        cor_p_value < edge_p_threshold & cor_FC > edge_FC_threshold ~ Enhanced_in_case,#"Enhanced in case",#"Up",
        cor_p_value < edge_p_threshold & cor_FC < 1 / edge_FC_threshold ~Enhanced_in_control,#"Enhanced in control",# "Down",
        TRUE ~ "Non-significant"
      ))
    results_df$from_to<-cor_data_both$from_to
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
  #  dplyr::filter(node %in% union(cor_data_diff$from,cor_data_diff$to))
  Differential_network_result <- new("Differential_network",
                                     group_name=comparison,
                                   #  Conditional_network_layout=Conditional_network_layout,
                                     diff_nodes=merge_node,
                                     diff_edges=cor_data_diff
  )
  
  return(Differential_network_result)
}


#' Standardize edge representation by sorting node pairs alphabetically to enable comparison.
#'
#' This helper function ensures that an edge from A to B is treated identically to B to A
#' by creating a canonical `from_to` identifier (e.g., "GeneB_GeneA" → always larger first).
#' It adds a `from_to` column and removes temporary sorting columns.
#'
#' @param edges A data frame with columns `from` and `to` representing undirected edges.
#' @return The same data frame with an additional `from_to` column in the format "larger_smaller".
#'
#' @importFrom dplyr select
#' @importFrom base transform
#'
#' @keywords internal
run_compare_edge<-function(edges){
  edges <-transform(edges, 
                    sorted_from = pmin(from, to), 
                    sorted_to = pmax(from, to)) 
  edges$from_to<-paste0(edges$sorted_to,"_",edges$sorted_from)
  edges<-edges |>
    dplyr::select(-sorted_from, -sorted_to)
  return(edges)
}