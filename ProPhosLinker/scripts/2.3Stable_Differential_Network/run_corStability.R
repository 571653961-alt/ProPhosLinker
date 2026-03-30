#' Stability Network Analysis and Related Classes
#'
#' This module provides classes and functions for performing stability analysis of correlation networks, including
#' deterministic bootstrap methods to evaluate the stability of network edges and filter out stable network structures.
#'
#' @name stable_network
NULL

#' Stability Network Class
#'
#' An S4 class for storing the results of network stability analysis.
#'
#' @slot group_name A string, the name of the analyzed group
#' @slot bootnet_result Complete bootstrap analysis results
#' @slot bootnet_result_filter Filtered stable network analysis results
#'
#' @export
setClass("StableNetwork",
         slots = c(
           group_name = "character",
           bootnet_result = "ANY",
           bootnet_result_filter = "ANY"
         )
)
#' Bootstrap Result Class
#'
#' An S4 class for storing complete bootstrap network analysis results.
#'
#' @slot bootnet_stable A data frame of bootstrap stability results
#' @slot bootnet_list A list containing results from all bootstrap iterations
#'
#' @export
setClass("BootnetResult",
         slots = c(
           bootnet_stable = "data.frame",
           bootnet_list = "list"
         )
)
#' Filtered Bootstrap Result Class
#'
#' An S4 class for storing filtered stable network analysis results.
#'
#' @slot bootnet_summary A filtered data frame summarizing the bootstrap results
#' @slot bootnet_stable A filtered data frame of bootstrap stability results
#' @slot bootnet_edge A data frame of filtered network edges
#' @slot bootnet_node A data frame of filtered network nodes
#'
#' @export
setClass("BootnetResultFilter",
         slots = c(
           bootnet_summary = "data.frame",
           bootnet_stable = "data.frame",
           bootnet_edge = "data.frame",
           bootnet_node = "data.frame"#,
           # bootnet_list = "list"
         )
)


#' Perform Correlation Network Stability Analysis
#'
#' This function uses bootstrap methods to evaluate the stability of edges in a correlation network by performing multiple resamplings to calculate
#' correlation coefficients and their confidence intervals, identifying stable network structures.
#'
#' @param count_table A data frame containing feature expression counts, must include a feature_ID column
#' @param group_name A string, the name of the analyzed group
#' @param annotation_table An optional node annotation table
#' @param cor_method The correlation method used ("pearson")
#' @param nBoots The number of bootstrap iterations, default is 50
#' @param nCores The number of CPU cores used for parallel processing
#' @param bootnet_R_threshold The threshold for correlation coefficients, used to filter weak edges
#' @param stability_threshold The confidence interval range threshold for stability evaluation
#' @param p_filter_table An optional precomputed p-value filtering table
#' @param bootnet_p_threshold The p-value threshold for statistical significance filtering
#' @param node_list An optional list of nodes to be analyzed
#' @param uniform_layout A logical value indicating whether to use a uniform node layout
#' @param keep_boot_num Numeric, number of bootstrap networks to keep for visualization (default 10). Set to reduce memory usage when nBoots is large.
#'
#' @return Returns a StableNetwork object containing complete bootstrap results and filtered stable networks
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example network data
#' data("stable_subnetwork_result")
#' data("metabolite_data")
#' filter_table <- stable_subnetwork_result@PreCor@filter_table
#' p_filter_table<-stable_subnetwork_result@PreCor@precor$edges |> dplyr::select(from,to,padj)
#'
#' # Perform stability analysis (In order to speed up the computation, nBoots is set to 5)
#' StableNetwork <- run_corStability(
#' count_table = filter_table, annotation_table = metabolite_data$annotation_table,
#' p_filter_table=p_filter_table,group_name = "T",
#' nBoots = 5,  stability_threshold = 0.4,bootnet_R_threshold=0.6)
#'
#' # Visualization
#' network_show(Network = StableNetwork,plot_type = 'overall_network',show_node_legend =TRUE)
#'
#' }
run_corStability<-function(count_table=NULL,group_name="data",annotation_table=NULL,
                           cor_method="spearman", nBoots=50, nCores=NULL,bootnet_R_threshold=0,
                           stability_threshold=0.3,p_filter_table=NULL,bootnet_p_threshold=0.05,
                           node_list=NULL,uniform_layout=FALSE, keep_boot_num = 10){#'cor_method' should be one of “cor”, “cov”, “cor_auto”, “npn”, “spearman”
  set.seed(123) # Ensure reproducibility for bootstrap resampling
  if(nrow(count_table)<2){
    message("Insufficient data rows for stability analysis. Computation terminated.\n")
    return(NULL)
  }

  #######################group_name
  if(length(group_name)>1){
    precor_group_name<-paste(group_name, collapse = "-vs-")
  }else{
    precor_group_name<-group_name
  }
  ##########################
  if(all(is.null(p_filter_table))){
    PreCor <- run_prenetwork(count_table = count_table, group_name = group_name,
                             cor_method = cor_method, R_threshold = 0, p_threshold = bootnet_p_threshold)
    p_filter_table=PreCor@precor$edges |> dplyr::select(from,to,padj)
  }

  #######################BootnetResult
  # # Fix: Force sort by feature_ID to ensure deterministic input for bootnet
  # count_table <- count_table[order(count_table$feature_ID), ]

  rownames(count_table)<-count_table$feature_ID
  count_matrix<-count_table |>
    dplyr::select(-feature_ID)
  count_matrix<-as.matrix(t(count_matrix))
    if(is.null(nCores)){
      nCores <- parallel::detectCores() - 1
    }else{
      nCores <- nCores
    }

    # Use new deterministic bootnet function
    bootnetresult <- run_bootnet(
        data = count_matrix,
        nBoots = nBoots,
        nCores = nCores,
        corMethod = cor_method,
        verbose = FALSE
    )

    # Full bootstrap table for statistics (T-tests, etc.)
    bootnet_sta <- bootnetresult$bootTable |>
      dplyr::filter(type=="edge") |>
      dplyr::filter(is.finite(value)) |> # Fix: Remove Infinite/NA values
      dplyr::mutate(num = as.numeric(gsub("boot ", "", name))) |>
      dplyr::arrange(num) |>
      dplyr::select(-num)

    # We need to manually compute summary since we don't have bootnet::summary anymore
    # Group by node1, node2 and compute mean, CI
    bootnet_sum <- bootnet_sta |>
      dplyr::group_by(node1, node2) |>
      dplyr::summarise(
          mean = mean(value, na.rm=TRUE),
          sd = sd(value, na.rm=TRUE),
          CIlower = quantile(value, 0.025, na.rm=TRUE),
          CIupper = quantile(value, 0.975, na.rm=TRUE),
          .groups = 'drop'
      ) |>
      dplyr::filter(is.finite(mean)) |>
          dplyr::arrange(node1, node2) # Ensure deterministic order

        bootnet_sum<-as.data.frame(bootnet_sum)
        # Fix: Add id column for matching
        bootnet_sum$id <- paste0(bootnet_sum$node1, "--", bootnet_sum$node2)

        bootnet_sum$CIrange<-bootnet_sum$CIupper-bootnet_sum$CIlower
        bootnet_sum$type <- "edge" # Add type column to match previous structure

  #R filter
  bootnet_sta_R<- bootnet_sta
  # Slim down for LIST only: Filter edges for visualization list to save memory
  bootnet_list_e0<- bootnet_sta_R |>
    dplyr::mutate(num = as.numeric(gsub("boot ", "", name))) |>
    dplyr::filter(num <= keep_boot_num) |> # Only keep limited networks in the list object
    dplyr::select(node1,node2,value, name) # Keep name for splitting

  colnames(bootnet_list_e0)<-c("from","to","cor", "name")
  bootnet_list_e0<-split(bootnet_list_e0 |> dplyr::select(-name),
                         factor(bootnet_list_e0$name, levels = unique(bootnet_list_e0$name)))

  # Optimization: Create the full nodes table ONCE, outside the loop
  # The union of all nodes in summary is sufficient to cover all edges
  all_nodes_vector <- union(bootnet_sum$node1, bootnet_sum$node2)
  base_nodes_df <- data.frame(node = all_nodes_vector) |>
      dplyr::left_join(annotation_table, by = c("node" = "feature_ID"))

  bootnet_list0 <- lapply(bootnet_list_e0, function(edges) {
    # Filter nodes for this specific bootstrap iteration if needed,
    # or just return the relevant subset.
    # Actually, the original code returned a node table specific to the edges of that iteration?
    # Original: nodes <-data.frame(node =union(bootnet_sum$node1, bootnet_sum$node2))
    # Wait, the original code used 'bootnet_sum' inside the loop, which is CONSTANT.
    # So it was creating the EXACT SAME node table in every iteration.
    # This optimization is definitely correct.

    list(nodes = base_nodes_df, edges = edges)
  })

  CI_cut<-stability_threshold

  bootnet_sum$CI_Status <- ifelse(bootnet_sum$CIrange < CI_cut, "stable", "Unstable")

  bootnet_sum_f<-bootnet_sum |>
    dplyr::filter(abs(mean)>bootnet_R_threshold) |>
    dplyr::filter(CI_Status =="stable")

  if(any(!is.null(node_list))){
    bootnet_sum_f<-bootnet_sum_f |>
      dplyr::filter(node1 %in% node_list) |>
      dplyr::filter(node2 %in% node_list)
  }
  if(nrow(bootnet_sum_f)==0){
    message("No stable edges, the operation is terminated.\n")
    return(NULL)
  }


  bootnet_sum_f<- bootnet_sum_f |>
    dplyr::mutate(from=node1,to=node2)

  bootnet_sum_f <-run_compare_edge(bootnet_sum_f) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::select(-from,-to)

  if(nrow(p_filter_table)==0){
    message("No stable edges, the operation is terminated.\n")
    return(NULL)
  }
  p_edges<-run_compare_edge(p_filter_table) |>
    dplyr::distinct(from_to, .keep_all = TRUE) |>
    dplyr::select(-from,-to)

  bootnet_sum_f <- dplyr::inner_join(bootnet_sum_f, p_edges, by = c("from_to"))
  if(nrow(bootnet_sum_f)==0){
    message("No significant edges, the operation is terminated.\n")
    return(NULL)
  }
  #######################

  bootnet_sum_e<- bootnet_sum_f |>#"CIrange"    "CI_Status"
    dplyr::select(node1,node2,mean,CIrange,padj)
  colnames(bootnet_sum_e)<-c("from","to","cor","CIrange","p_adjust")

  # Fix: Force sort edges to ensure deterministic graph construction downstream
  bootnet_sum_e <- bootnet_sum_e |>
    dplyr::arrange(from, to)

  if(uniform_layout){

    bootnet_sum_n<-data.frame(node =union(bootnet_sum$node1, bootnet_sum$node2))  |>
      dplyr::left_join(annotation_table,by=c("node"="feature_ID"))
    if(any(!is.null(node_list))){
      bootnet_sum_n1<-data.frame(node=node_list) |>
        dplyr::left_join(bootnet_sum_n,by=c("node"="node"))
      bootnet_sum_n=bootnet_sum_n1
    }
  }else{
    bootnet_sum_n<-data.frame(node = union(bootnet_sum_e$from,bootnet_sum_e$to)) |>
      dplyr::left_join(annotation_table,by=c("node"="feature_ID"))
  }

  # Fix: Force sort nodes
  bootnet_sum_n <- bootnet_sum_n |>
    dplyr::arrange(node)



  bootnet_sta_f<-  bootnet_sta |>
    dplyr::filter(id %in% bootnet_sum_f$id) |>
    dplyr::filter(abs(value)>bootnet_R_threshold)

  # Optimization: Convert to Wide Format to save space (reduces rows by factor of nBoots)
  bootnet_sta_wide <- bootnet_sta |>
    tidyr::pivot_wider(id_cols = c(graph, type, node1, node2, id), names_from = name, values_from = value)
    
  bootnet_sta_f_wide <- bootnet_sta_f |>
    tidyr::pivot_wider(id_cols = c(graph, type, node1, node2, id), names_from = name, values_from = value)

  bootnet_result <- new("BootnetResult",
                        bootnet_stable = bootnet_sta_wide,
                        bootnet_list= bootnet_list0
  )
  bootnet_result_filter <- new("BootnetResultFilter",
                               bootnet_summary = bootnet_sum_f,
                               bootnet_stable = bootnet_sta_f_wide,
                               bootnet_edge = bootnet_sum_e,
                               bootnet_node = bootnet_sum_n
  )

  stable_network <- new("StableNetwork",
                        group_name = precor_group_name,
                        bootnet_result = bootnet_result,
                        bootnet_result_filter = bootnet_result_filter
  )
  invisible(stable_network)
}


#' Deterministic Network Estimation (Correlation)
#'
#' A simplified and deterministic wrapper for correlation network estimation.
#'
#' @param data Data matrix (Samples x Features)
#' @param corMethod Correlation method ("spearman", "pearson")
#' @param verbose Logical
#' @export
run_estimateNetwork <- function(data, corMethod = "spearman", verbose = FALSE) {
  
  # 1. Compute Correlation
  # Handle "cor_auto" if needed, otherwise standard cor
  if (corMethod == "cor_auto") {
      corMat <- qgraph::cor_auto(data, verbose = verbose)
  } else {
      corMat <- stats::cor(data, method = corMethod, use = "pairwise.complete.obs")
  }
  
  # 2. Safety Cleanup
  # Replace NA/Inf with 0 to prevent downstream crashes (e.g. eigen error)
  if (any(!is.finite(corMat))) {
      corMat[!is.finite(corMat)] <- 0
  }
  
  # Return the correlation matrix as the network weights
  return(corMat)
}

#' Deterministic Bootstrap for Networks
#'
#' A fully reproducible bootstrap function.
#'
#' @param data Data matrix (Samples x Features)
#' @param nBoots Number of bootstraps
#' @param nCores Number of cores
#' @param seed Base random seed
#' @param corMethod Correlation method
#' @export
run_bootnet <- function(data, nBoots = 100, nCores = 1, seed = 123, corMethod = "spearman", verbose = TRUE) {
  
  nSample <- nrow(data)
  nNode <- ncol(data)
  labels <- colnames(data)
  
  # 1. Estimate Original Network
  if(verbose) message("Estimating original network...")
  sample_graph <- run_estimateNetwork(data, corMethod = corMethod)
  
  # 2. Prepare Seeds
  # Crucial for reproducibility
  if (!exists(".Random.seed", envir = .GlobalEnv)) set.seed(NULL)
  old_seed <- .Random.seed
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  
  set.seed(seed)
  boot_seeds <- sample.int(1000000, nBoots)
  
  # 3. Bootstrap Loop
  if(verbose) message(paste0("Bootstrapping (", nBoots, " iterations)..."))
  
  # Define single run function
  run_one_boot <- function(i) {
      set.seed(boot_seeds[i]) # Local seed
      
      # Resample
      indices <- sample(1:nSample, nSample, replace = TRUE)
      boot_data <- data[indices, ]
      
      # Estimate
      tryCatch({
          net <- run_estimateNetwork(boot_data, corMethod = corMethod)
          # We only need the upper triangle to save space/time for symmetry
          # But keeping full matrix is easier for now
          return(list(success = TRUE, graph = net))
      }, error = function(e) {
          return(list(success = FALSE, error = e$message))
      })
  }
  
  # Execution
  if (nCores > 1) {
      cl <- parallel::makeCluster(nCores)
      # Export needed variables
      parallel::clusterExport(cl, c("run_estimateNetwork", "boot_seeds", "data", "nSample", "corMethod"), envir = environment())
      # If qgraph/stats needed
      res_list <- parallel::parLapply(cl, 1:nBoots, run_one_boot)
      parallel::stopCluster(cl)
  } else {
      res_list <- lapply(1:nBoots, run_one_boot)
  }
  
  # 4. Process Results into Bootnet-like format
  # We need to construct the 'bootTable' (long format)
  
  success_indices <- which(sapply(res_list, function(x) x$success))
  n_success <- length(success_indices)
  
  if (n_success < nBoots) {
      warning(paste("Only", n_success, "out of", nBoots, "bootstraps succeeded."))
  }
  
  # Optimization: Pre-calculate static vectors
  # Since all networks have the same nodes in the same order, the upper triangle indices
  # and the resulting node labels/IDs are identical for every bootstrap.
  
  # Use the sample graph to define the structure
  mat_template <- sample_graph
  ut <- upper.tri(mat_template)
  
  # Get indices
  row_inds <- row(mat_template)[ut]
  col_inds <- col(mat_template)[ut]
  
  # Pre-calculate node labels and IDs
  node1_vec <- labels[row_inds]
  node2_vec <- labels[col_inds]
  id_vec <- paste0(node1_vec, "--", node2_vec)
  
  # Pre-allocate list for speed
  boot_table_list <- vector("list", n_success)
  
  for (k in seq_along(success_indices)) {
      idx <- success_indices[k]
      mat <- res_list[[idx]]$graph
      
      # Extract values only
      values <- mat[ut]
      
      # Create dataframe efficiently
      # Re-use the pre-calculated vectors
      # Use rep.int for the name column which is much faster than paste0 on every row
      n_edges <- length(values)
      name_val <- paste0("boot ", idx)
      
      df <- data.frame(
          graph = "1",
          type = "edge",
          node1 = node1_vec,
          node2 = node2_vec,
          value = values,
          name = name_val, # R will automatically recycle this
          id = id_vec
      )
      boot_table_list[[k]] <- df
  }
  
  bootTable <- dplyr::bind_rows(boot_table_list)
  
  # Create sampleTable (Original network)
  # Re-use vectors again
  sampleTable <- data.frame(
      graph = "1",
      type = "edge",
      node1 = node1_vec,
      node2 = node2_vec,
      value = sample_graph[ut],
      name = "sample",
      id = id_vec
  )
  
  # Return a structured list that mocks a bootnet object
  return(list(
      bootTable = bootTable,
      sampleTable = sampleTable,
      sample = list(graph = sample_graph, nNode = nNode, nPerson = nSample, labels = labels),
      metric = "edge",
      nBoots = nBoots,
      nCores = nCores
  ))
}

