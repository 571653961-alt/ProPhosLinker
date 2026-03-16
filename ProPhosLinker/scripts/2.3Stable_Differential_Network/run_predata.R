#' Prepare and validate input data for downstream analysis by aligning count table with sample metadata.
#'
#' This function integrates a feature count table (e.g., gene expression, microbial abundance)
#' with a sample annotation list (`samplelist`), performs basic validation, removes duplicates,
#' ensures column consistency, and returns a standardized S4 object of class `PreData`.
#' It also checks for minimal sample size and verifies required columns.
#'
#' @param count_table A data frame where rows represent features (must contain a column named `feature_ID`)
#'                    and columns represent samples. All sample columns should be numeric or coercible to numeric.
#' @param samplelist A data frame containing at least two columns: `sample` (sample identifiers) and `group`
#'                   (group labels). Duplicate samples are removed, keeping the first occurrence.
#' @param group_name Optional character vector specifying which group(s) to retain from `samplelist`.
#'                   If multiple names are provided, they are collapsed into a single label using "-vs-".
#'                   If `NULL`, all samples in `samplelist` are used and the group name defaults to `"data"`.
#'
#' @return An object of S4 class `PreData` with slots:
#'   \describe{
#'     \item{group_name}{A formatted character string representing the selected group(s).}
#'     \item{samplelist}{The filtered and deduplicated sample metadata data frame.}
#'     \item{count_table}{The cleaned count table containing only `feature_ID` and selected samples, with numeric values.}
#'   }
#'
#' @importFrom dplyr distinct filter select mutate across all_of sym
#' @import dplyr
#'
#' @examples
#' # Example count table
#' count_table <- data.frame(
#'   feature_ID = c("Gene1", "Gene2", "Gene3"),
#'   S1 = c(10, 20, 30),
#'   S2 = c(15, 25, 35),
#'   S3 = c(12, 22, 32)
#' )
#' 
#' # Example sample list
#' samplelist <- data.frame(
#'   sample = c("S1", "S2", "S3"),
#'   group = c("Control", "Control", "Treatment")
#' )
#' 
#' predata_obj <- run_predata(count_table, samplelist, group_name = "Control")
#' print(predata_obj@group_name)
#' head(predata_obj@count_table)
#'
#' @export
setClass("PreData", slots = c(
  group_name = "character",
  samplelist = "data.frame",
  count_table = "data.frame"
))

run_predata<-function(count_table=NULL,samplelist=NULL,group_name=NULL){
  if(is.null(group_name)){
    precor_group_name<-"data"
  }else{
    if(length(group_name)>1){
      precor_group_name<-paste(group_name, collapse = "-vs-")
    }else{
      precor_group_name<-group_name
    }
  }
  
  if(all(c("group","sample") %in% colnames(samplelist))){
    samplelist<-samplelist |>
      dplyr::distinct(!!dplyr::sym("sample"), .keep_all = TRUE) 
  }else{
    stop("No group column or sample column in the samplelist, please check it.")
  }
  if(!is.null(group_name)){
    samplelist<-samplelist |> #Remove duplicate samples
      dplyr::filter(group %in% group_name)
  }
  if(nrow(samplelist)==0){
    stop("No sample in the group, please check the group information.")
  }
  if(nrow(samplelist)<5){
    warning("The sample size is too small and may have an impact on the stability test.")
  }
  raw_table<-as.data.frame(count_table)
  if("feature_ID" %in% colnames(raw_table)){
    tryCatch({
      raw_table <- raw_table |>
        dplyr::distinct(!!dplyr::sym("feature_ID"), .keep_all = TRUE) |> #Remove duplicate ids
        dplyr::select("feature_ID",all_of(samplelist$sample))  |>
        dplyr::mutate(across(-!!dplyr::sym("feature_ID"), as.numeric))
    }, error = function(e) {
      message("Error in matching quantitative tables and sample tables: ", e$message)
    })
  }else{
    stop("No feature_ID column in the count_table, please check it.")
  }
  ################################
  Predata<-new("PreData", 
               group_name = precor_group_name,
               samplelist = samplelist,
               count_table = raw_table
               
  )
  return(Predata)#unique_classes,colors,tits_nodesizes,
}