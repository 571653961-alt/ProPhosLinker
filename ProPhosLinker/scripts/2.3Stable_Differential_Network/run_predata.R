#dplyr
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
