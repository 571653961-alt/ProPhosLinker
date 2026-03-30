#' Prepare and validate data for pre-network analysis
#'
#' This function validates and preprocesses count tables and sample metadata
#' for downstream network analysis. It filters samples based on group information,
#' removes duplicates, ensures proper data types, and returns a structured PreData
#' S4 object.
#'
#' @param count_table data.frame, containing a 'feature_ID' column (feature identifiers)
#'                    and abundance values for each sample as additional columns.
#' @param samplelist data.frame, sample metadata containing at least two columns:
#'                    'sample' (sample identifiers matching count_table column names)
#'                    and 'group' (group assignments for each sample).
#' @param group_name character, optional. Group name(s) to filter samples.
#'                    If NULL, all samples are kept. If length > 1, group names
#'                    are joined with "-vs-" for labeling.
#'
#' @return A PreData S4 object with the following slots:
#'   \item{group_name}{character, processed group label for result identification}
#'   \item{samplelist}{data.frame, filtered and deduplicated sample metadata}
#'   \item{count_table}{data.frame, filtered abundance table with feature_ID column
#'                      and sample columns matching the filtered samplelist}
#'
#' @details 
#' Processing steps:
#' 1. If group_name is provided, filter samplelist to include only specified groups.
#'    Otherwise, set default group name to "data".
#' 2. Validate samplelist contains required 'sample' and 'group' columns.
#' 3. Remove duplicate samples from samplelist (keeping first occurrence).
#' 4. Check for empty samplelist after filtering; stop if no samples remain.
#' 5. Issue warning if sample size < 5 (may affect statistical stability).
#' 6. Validate count_table contains 'feature_ID' column.
#' 7. Remove duplicate feature IDs from count_table (keeping first occurrence).
#' 8. Subset count_table to only include samples present in filtered samplelist.
#' 9. Convert all abundance columns to numeric type.
#'
#' @importFrom dplyr distinct filter select mutate across sym all_of
#' @importFrom rlang !! sym
#'
#' @examples
#' \dontrun{
#' # Prepare data with group filtering
#' sample_meta <- data.frame(
#'   sample = c("S1", "S2", "S3", "S4"),
#'   group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' 
#' count_data <- data.frame(
#'   feature_ID = c("Gene1", "Gene2", "Gene3"),
#'   S1 = c(10, 20, 30),
#'   S2 = c(12, 22, 32),
#'   S3 = c(15, 25, 35),
#'   S4 = c(18, 28, 38)
#' )
#' 
#' # Keep only Control group
#' pre_data <- run_predata(
#'   count_table = count_data,
#'   samplelist = sample_meta,
#'   group_name = "Control"
#' )
#' 
#' # Access filtered data
#' head(pre_data@count_table)
#' pre_data@samplelist
#' }
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
