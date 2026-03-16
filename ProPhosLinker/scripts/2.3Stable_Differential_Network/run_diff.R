#' @title analysis of variance
#'
#' @param quantitative_table Quantitative table  feature_ID sample1 sample2...
#' @param samplelist The sample table contains two columns: sample and group
#' @param compare_group
#' @param p_threshold
#' @param FC_threshold
#' @param paired
#' @return Analysis of variance table 
#' @export
#'
#' @examples

run_diff<-function(quantitative_table=NULL,samplelist=NULL,compare_group=NULL,p_threshold=0.05,FC_threshold=1.2,p_value_type="q_value",paired = FALSE){
  # quantitative_table
  quantitative_table<-as.data.frame(quantitative_table)

  if("feature_ID" %in% colnames(quantitative_table)){
  core_table <- quantitative_table |>
    dplyr::distinct(!!dplyr::sym("feature_ID"), .keep_all = TRUE) |>
    tidyr::drop_na()  
  }else{
    stop("No feature_ID column in the quantitative_table, no differential analysis will be performed.")
  }
 
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
  
  # samplelist
  if(all(c("group","sample") %in% colnames(samplelist))){
  compare_table <- tibble::as_tibble_col(cgroup) |>  
    dplyr::group_by(.data$value) |>   
    dplyr::mutate(
      sample_list = list(
        samplelist[samplelist$group %in% unlist(stringr::str_split(.data$value, "\\+")), ]$sample
      )
    ) |>
    dplyr::mutate(sample_count = purrr::map_int(.data$sample_list, length))
  }else{
    stop("No group column or sample column in the samplelist, no differential analysis will be performed.")
  }
  diffdata<-statistic(compare_table,mean_table, core_table, paired = paired)
  # Volcano
  if (all(is.na(diffdata$p_value))) {
    diffdata <- diffdata |>
      dplyr::mutate(
        State = dplyr::case_when(
          FC > FC_threshold ~ "Up",
          FC < 1 / FC_threshold ~ "Down",
          TRUE ~ "Non-significant"
        )
      )
  } else {
    diffdata <- diffdata |>
      dplyr::mutate(
        State = dplyr::case_when(
          !!dplyr::sym(p_value_type) < p_threshold & FC > FC_threshold ~ "Up",
          !!dplyr::sym(p_value_type) < p_threshold & FC < 1 / FC_threshold ~ "Down",
          TRUE ~ "Non-significant"
        )
      )
  }
  diffdata<-diffdata |>
    dplyr::relocate(feature_ID, FC,Log2FC, p_value, q_value, State, .before = 1) |>
    dplyr::arrange(!!dplyr::sym(p_value_type)) 
  return(diffdata)
  
}


statistic <- function(compare_table,mean_table, core_table, paired = FALSE) {
  if (paired & compare_table$sample_count[1] == compare_table$sample_count[2]) {
    if (compare_table$sample_count[1] != 0 | compare_table$sample_count[2] != 0) {
      res <- core_table |>
        dplyr::rowwise() |> 
        dplyr::mutate(
          p_value = t.test(
            x = dplyr::c_across(all_of(compare_table$sample_list[[1]])),
            y = dplyr::c_across(all_of(compare_table$sample_list[[2]])),
            paired = TRUE  
          )$p.value,
          .keep = "none", 
          feature_ID = .data$feature_ID  
        ) |>
        dplyr::mutate(
          q_value = p.adjust(.data$p_value, method = "BH")
        )
    } else {
      stop("No sample in the group, please check the group information.")
    }
  } else {
  
    if (compare_table$sample_count[1] == 0 | compare_table$sample_count[2] == 0) {
      stop("No sample in the group, please check the group information.")
    } else if (
    
      compare_table$sample_count[1] == 1 && compare_table$sample_count[2] == 1
    ) {
      message("Only FC")   
      res <- core_table |> dplyr::mutate(
        p_value = NA, q_value = NA, .keep = "none",
        feature_ID = .data$feature_ID
      )
     
    } else if (
      compare_table$sample_count[1] == 1 && compare_table$sample_count[2] > 1
    ) {
      res <- core_table |>
        dplyr::rowwise() |>
        dplyr::mutate(
          p_value = t.test(
            x = dplyr::c_across(all_of(compare_table$sample_list[[2]])),
            mu = .data[[compare_table$sample_list[[1]]]],
          )$p.value,
          .keep = "none",
          feature_ID = .data$feature_ID
        ) |>
        dplyr::mutate(
          q_value = p.adjust(.data$p_value, method = "BH")
        )
    } else if (compare_table$sample_count[1] > 1 && compare_table$sample_count[2] == 1) {
      res <- core_table |>
        dplyr::rowwise() |>
        dplyr::mutate(
          p_value = t.test(
            x = dplyr::c_across(all_of(compare_table$sample_list[[1]])),
            mu = .data[[compare_table$sample_list[[2]]]],
          )$p.value,
          .keep = "none",
          feature_ID = .data$feature_ID
        ) |>
        dplyr::mutate(
          q_value = p.adjust(.data$p_value, method = "BH")
        )
    } else {
      res <- core_table |>
        dplyr::rowwise() |>
        dplyr::mutate(
          p_value = t.test(
            x = dplyr::c_across(all_of(compare_table$sample_list[[1]])),
            y = dplyr::c_across(all_of(compare_table$sample_list[[2]])),
          )$p.value,
          .keep = "none",
          feature_ID = .data$feature_ID
        ) |>
        dplyr::mutate(
          q_value = p.adjust(.data$p_value, method = "BH")
        )
    }
  }

  result <- core_table |>
    dplyr::mutate(
      FC = rowMeans(pick(compare_table$sample_list[[1]])) /
        rowMeans(pick(compare_table$sample_list[[2]])),
      .keep = "none"
    ) |>
    dplyr::bind_cols(res) |> 
    dplyr::mutate(Log2FC=as.numeric(log2(FC)))

 
  return(result)
}


run_mean<-function(count_table=NULL,samplelist=NULL,compare_group=NULL){
  count_table<-as.data.frame(count_table)
  if("feature_ID" %in% colnames(count_table)){
    core_table <- count_table |>
      dplyr::distinct(!!dplyr::sym("feature_ID"), .keep_all = TRUE) |>
      tidyr::drop_na()
  }else{
    stop("No feature_ID column in the count_table, no differential analysis will be performed.")
  }
  
  ###mean_table

  mean_dat<-as.data.frame(core_table)
  rownames(mean_dat)<-mean_dat$feature_ID  
  mean_dat<-mean_dat |>
    dplyr::select(-1) |>
    t() |>
    as.data.frame(rownames.force = TRUE,colnames.force = TRUE)
  mean_table<- apply(mean_dat,2,function(x){
    if(max(x)!=0){
      x/max(x)
    }else{
      x 
    }
    
    })
  mean_table<-mean_table |>
    t() |> 
    as.data.frame(rownames.force = TRUE,colnames.force = TRUE)
  mean_table$feature_ID<-rownames(mean_table)
  mean_table<-mean_table |>
    dplyr::select(feature_ID, everything())
  
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
  
  #samplelist
  if(all(c("group","sample") %in% colnames(samplelist))){
    compare_table <- tibble::as_tibble_col(cgroup) |>
      dplyr::group_by(.data$value) |>
      dplyr::mutate(
        sample_list = list(
          samplelist[samplelist$group %in% unlist(stringr::str_split(.data$value, "\\+")), ]$sample
        )
      ) |>
      dplyr::mutate(sample_count = purrr::map_int(.data$sample_list, length)) 
  }else{
    stop("No group column or sample column in the samplelist, no differential analysis will be performed.")
  }
  
  result_mean<-mean_table |>
    dplyr::mutate(
      case_norm_mean=rowMeans(pick(compare_table$sample_list[[1]])) ,
      control_norm_mean=rowMeans(pick(compare_table$sample_list[[2]]))
      )  |>
    dplyr::select(feature_ID,case_norm_mean,control_norm_mean)
  return(result_mean)
}

process_intensity <- function(intensity1_str, intensity2_str) {
  vec1 <- suppressWarnings(as.numeric(str_split(intensity1_str, ";", simplify = TRUE)))
  vec2 <- suppressWarnings(as.numeric(str_split(intensity2_str, ";", simplify = TRUE)))
  
  max_total <- max(c(vec1,vec2), na.rm = TRUE)
  
  if (is.finite(max_total) && max_total > 0) {
    norm1 <- vec1 / max_total
    norm2 <- vec2 / max_total
  } else {
    norm1 <- rep(NA, length(vec1))
    norm2 <- rep(NA, length(vec2))
  }
  
  mean1 <- mean(norm1, na.rm = TRUE)
  mean2 <- mean(norm2, na.rm = TRUE)
  
  if (is.nan(mean1)) mean1 <- NA
  if (is.nan(mean2)) mean2 <- NA
  
  return(list(mean_intensity1 = mean1, mean_intensity2 = mean2))
}


run_process_intensity<-function(diff_table){
  diff_table <- diff_table %>%
    rowwise() %>%  
    mutate(
      result = list(process_intensity(`Intensity(1)`, `Intensity(2)`)),
      case_norm_mean = result$mean_intensity2,
      control_norm_mean = result$mean_intensity1
    ) %>%
    select(-result)
}


