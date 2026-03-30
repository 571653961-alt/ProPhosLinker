#' @title analysis of variance

#' Perform differential analysis between experimental and control groups
#'
#' This function conducts differential expression/abundance analysis between
#' two groups (experimental vs control) using t-tests. It calculates fold change,
#' p-values, and adjusted q-values, then classifies features as up-regulated,
#' down-regulated, or non-significant based on user-defined thresholds.
#'
#' @param quantitative_table data.frame, containing a 'feature_ID' column and
#'                           abundance values for each sample as additional columns.
#' @param samplelist data.frame, sample metadata containing 'sample' and 'group'
#'                   columns for group assignment.
#' @param compare_group character, comparison specification in format
#'                      "Experimental_group:Control_group". Groups can be combined
#'                      using "+" (e.g., "Treatment1+Treatment2:Control").
#' @param p_threshold numeric, significance threshold for differential features.
#'                     Default 0.05.
#' @param FC_threshold numeric, fold change threshold. Features with |FC| > threshold
#'                      are considered differentially expressed. Default 1.2.
#' @param p_value_type character, which p-value to use for filtering: "p_value" or
#'                      "q_value". Default "q_value".
#' @param paired logical, whether to perform paired t-test. Default FALSE.
#'
#' @return A data.frame containing differential analysis results with columns:
#'   \item{feature_ID}{Feature identifiers}
#'   \item{FC}{Fold change (Experimental/Control)}
#'   \item{Log2FC}{Log2 transformed fold change}
#'   \item{p_value}{Raw p-value from t-test}
#'   \item{q_value}{BH-adjusted p-value}
#'   \item{State}{Classification: "Up", "Down", or "Non-significant"}
#'
#' @details
#' Analysis steps:
#' 1. Validate and deduplicate quantitative_table by feature_ID.
#' 2. Parse compare_group to extract experimental and control group definitions.
#' 3. Filter samplelist to include samples from specified groups.
#' 4. Perform appropriate t-test (paired, one-sample, or two-sample) based on
#'    sample sizes and paired parameter.
#' 5. Calculate fold change as ratio of group means.
#' 6. Classify features based on p/q-value and FC thresholds.
#' 7. Return sorted results by significance.
#'
#' @importFrom dplyr distinct mutate relocate arrange rowwise
#' @importFrom tidyr drop_na
#' @importFrom tibble as_tibble_col
#' @importFrom purrr map_int
#' @importFrom stringr str_split
#'
#' @examples
#' \dontrun{
#' # Differential analysis between Treatment and Control groups
#' diff_results <- run_diff(
#'   quantitative_table = abundance_data,
#'   samplelist = sample_metadata,
#'   compare_group = "Treatment:Control",
#'   p_threshold = 0.05,
#'   FC_threshold = 1.5,
#'   p_value_type = "q_value"
#' )
#' 
#' # View up-regulated features
#' subset(diff_results, State == "Up")
#' }
#'
#' @export

run_diff<-function(quantitative_table=NULL,samplelist=NULL,compare_group=NULL,p_threshold=0.05,FC_threshold=1.2,p_value_type="q_value",paired = FALSE){
  ###quantitative_table
  quantitative_table<-as.data.frame(quantitative_table)
  if("feature_ID" %in% colnames(quantitative_table)){
  core_table <- quantitative_table |>
    dplyr::distinct(!!dplyr::sym("feature_ID"), .keep_all = TRUE) |>
    tidyr::drop_na()
  }else{
    stop("No feature_ID column in the quantitative_table, no differential analysis will be performed.")
  }
 
  #compare_group:Experimental group:control group
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
  
  ##samplelist
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

#' Perform statistical testing for differential analysis
#'
#' Internal function that performs t-tests between groups based on sample sizes
#' and experimental design. Handles various scenarios including paired tests,
#' one-sample tests, and standard two-sample t-tests.
#'
#' @param compare_table data.frame, containing group information with columns:
#'                      'value' (group name), 'sample_list' (list of sample names),
#'                      and 'sample_count' (number of samples per group).
#' @param mean_table data.frame, not used in current implementation (reserved).
#' @param core_table data.frame, feature abundance table with 'feature_ID' column.
#' @param paired logical, whether to perform paired t-test. Default FALSE.
#'
#' @return A data.frame with statistical results containing:
#'   \item{feature_ID}{Feature identifiers}
#'   \item{p_value}{Raw p-values from t-tests}
#'   \item{q_value}{BH-adjusted p-values}
#'   \item{FC}{Fold change (Experimental/Control)}
#'   \item{Log2FC}{Log2 transformed fold change}
#'
#' @details
#' Testing scenarios:
#' \itemize{
#'   \item Paired t-test: When paired = TRUE and both groups have equal sample sizes
#'   \item One-sample t-test: When one group has single sample (testing against that value)
#'   \item Two-sample t-test: Standard Welch's t-test for independent samples
#'   \item NA values: When either group has no samples, or both groups have single sample
#' }
#'
#' @keywords internal
#'
#' @importFrom dplyr rowwise mutate bind_cols c_across
#' @importFrom rlang .data
#'
#' @noRd
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
    } else if (
      compare_table$sample_count[1] > 1 && compare_table$sample_count[2] == 1
    ) {
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
    # dplyr::bind_cols(
    #   core_table |>
    #     dplyr::select(
    #       compare_table$sample_list[[1]],
    #       compare_table$sample_list[[2]]
    #     )
    # )

 
  return(result)
}

#' Calculate normalized mean abundances for comparison groups
#'
#' This function normalizes abundance values by scaling each feature to a maximum
#' of 1 (feature-wise normalization) and computes mean normalized abundances for
#' experimental and control groups.
#'
#' @param count_table data.frame, containing a 'feature_ID' column and abundance
#'                    values for each sample as additional columns.
#' @param samplelist data.frame, sample metadata containing 'sample' and 'group'
#'                   columns for group assignment.
#' @param compare_group character, comparison specification in format
#'                      "Experimental_group:Control_group". Groups can be combined
#'                      using "+" (e.g., "Treatment1+Treatment2:Control").
#'
#' @return A data.frame with normalized mean abundances:
#'   \item{feature_ID}{Feature identifiers}
#'   \item{case_norm_mean}{Mean normalized abundance in experimental group}
#'   \item{control_norm_mean}{Mean normalized abundance in control group}
#'
#' @details
#' Normalization steps:
#' 1. Validate and deduplicate count_table by feature_ID.
#' 2. Transpose data to sample × feature matrix.
#' 3. Apply feature-wise normalization: scale each feature to [0,1] range by
#'    dividing by its maximum value.
#' 4. Transpose back to feature × sample format.
#' 5. Parse compare_group to identify experimental and control group samples.
#' 6. Calculate mean normalized abundance for each group.
#'
#' @importFrom dplyr distinct mutate select rowwise pick
#' @importFrom tidyr drop_na
#' @importFrom tibble as_tibble_col
#' @importFrom purrr map_int
#' @importFrom stringr str_split
#'
#' @examples
#' \dontrun{
#' # Calculate normalized means for Treatment vs Control
#' mean_abundance <- run_mean(
#'   count_table = abundance_data,
#'   samplelist = sample_metadata,
#'   compare_group = "Treatment:Control"
#' )
#' 
#' # View results
#' head(mean_abundance)
#' }
#'
#' @export
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
    dplyr::select(-1) |>  # 去掉第一列（如化合物ID等）
    t() |>   # 转置
    as.data.frame(rownames.force = TRUE,colnames.force = TRUE)   # 强制保留行名为数据框的行名
  mean_table<- apply(mean_dat,2,function(x){
    #x[x==0.00000000000001]<-NA #如果count_table做了fix缺失值填充
    #scales::rescale(x, to = c(0, 1))
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
    dplyr::select(feature_ID, everything())     # 把 feature_ID 放在第一列
  
  #compare_group:Experimental group:control group
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
  
  ##samplelist
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


#' Process intensity values from semicolon-separated strings
#'
#' Internal function that converts semicolon-separated intensity strings to numeric
#' vectors, normalizes by the maximum value across both vectors, and computes mean
#' normalized intensities.
#'
#' @param intensity1_str character, semicolon-separated intensity values (e.g., "1.2;2.3;3.4")
#' @param intensity2_str character, semicolon-separated intensity values
#'
#' @return A list containing:
#'   \item{mean_intensity1}{Mean normalized intensity for first vector}
#'   \item{mean_intensity2}{Mean normalized intensity for second vector}
#'
#' @keywords internal
#'
#' @importFrom stringr str_split
#'
#' @noRd
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

#' Apply intensity processing to differential analysis results
#'
#' This function processes semicolon-separated intensity columns (Intensity(1) and
#' Intensity(2)) in a differential analysis table, converting them to normalized
#' mean intensities for case and control groups.
#'
#' @param diff_table data.frame, differential analysis results containing columns
#'                   'Intensity(1)' and 'Intensity(2)' with semicolon-separated values.
#'
#' @return A data.frame with original columns plus two additional columns:
#'   \item{case_norm_mean}{Mean normalized intensity for case/experimental group}
#'   \item{control_norm_mean}{Mean normalized intensity for control group}
#'   The original 'Intensity(1)' and 'Intensity(2)' columns are preserved.
#'
#' @details
#' Processing steps:
#' 1. For each row, split Intensity(1) and Intensity(2) strings into numeric vectors.
#' 2. Find maximum value across both vectors.
#' 3. Normalize both vectors by dividing by the maximum value.
#' 4. Calculate mean of normalized values for each vector.
#' 5. Add results as new columns to the original table.
#'
#' @importFrom dplyr rowwise mutate select
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Process intensity columns in differential results
#' diff_with_intensity <- run_process_intensity(diff_table)
#' 
#' # View normalized means
#' diff_with_intensity %>%
#'   select(feature_ID, case_norm_mean, control_norm_mean)
#' }
#'
#' @export
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


