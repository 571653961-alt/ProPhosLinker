
run_mediation<-function(count_table=NULL,samplelist=NULL,x=NULL,y=NULL,m=NULL,group_name="data",
                        mediation_p_threshold=0.05,outdir="./"){
  tryCatch({
    if(length(group_name)>1){
      group_name<-paste(group_name, collapse = "-vs-")
    }else{
      group_name<-gsub(":","-vs-", group_name)#diff
    }

    ##samplelist
    if(all(c("group","sample") %in% colnames(samplelist))){
      samplelist<-samplelist |>
        dplyr::filter(group %in% group_name)
    }else{
      stop("No group column or sample column in the samplelist, no mediation analysis will be performed.")
    }
    if(nrow(samplelist)==0){
      stop("Please check that the group column in the samplelist matches the group_name parameter.")
    }
    if("feature_ID" %in% colnames(count_table)){
      count_table <- as.data.frame(count_table) |>
        dplyr::distinct(!!dplyr::sym("feature_ID"), .keep_all = TRUE)  |>
        dplyr::select("feature_ID",all_of(samplelist$sample))
    }else{
      stop("No feature_ID column in the count_table, no mediation analysis will be performed.")
    }
    if(nrow(count_table)==0){
      stop("Please check that sample names in the count_table matches the sample column in the samplelist")
    }
    count_table<-count_table |>
      dplyr::filter(feature_ID %in% c(x,m,y))
    if(nrow(count_table)!=3){
      stop("Please check that feature_ID column in the count_table matches the x,m,y parameter.")
    }
    data<-  count_table |>
      dplyr::select(-feature_ID)
    # 只保留数值型数据（去除 feature_ID 列）
    data_matrix <- as.matrix(t(data))
    colnames(data_matrix)<-count_table$feature_ID

    # 拟合模型
    f <- paste(y, "~", x, "+ (",m,")")
    model_names <-paste("mediation",x,  m, y,group_name, sep = "_")
    png(file.path(outdir,paste0(model_names,".png")),width = 2000,height = 2000,res=300)
    model <-psych::mediate(as.formula(f), data = data_matrix)
    dev.off()
    # 设置筛选阈值
  #  R2_threshold<-mediation_R2_threshold
        df<- data.frame(
          Model =paste(x,m, y, sep = "_") ,
          a=round(as.numeric(model[["a"]][1]),2),
          b=round(as.numeric(model[["b"]][1]),2),
          ab=round(as.numeric(model[["ab"]][1]),2),
          c=round(as.numeric(model[["c"]][1]),2),
          c1= round(as.numeric(model[["indirect"]][1]),2),
          model_r2=round(as.numeric(model[["cprime.reg"]][["R2"]]),2),#模型解释力
          c_p   = round(as.numeric(model[["total.reg"]][["prob"]][2]),4),  # 总效应 p 值(X 对 Y 的总影响（不考虑中介变量）)
          c1_p  = round(as.numeric(model[["cprime.reg"]][["prob"]][2]),4),    # 直接效应 p 值(X 和 M 同时对 Y 的影响,2是X，3是M)
          a_p=  round(as.numeric(model[["a.reg"]][["prob"]][2]),4),
          b_p=  round(as.numeric(model[["b.reg"]][["prob"]][2]),4),
          indirect_lower   = round(as.numeric(model[["boot"]][["ci.ab"]][1]),2),        # 间接效应下限
          indirect_upper   =round(as.numeric(model[["boot"]][["ci.ab"]][2]),2)       # 间接效应上限
        )
    mediation_table<-add_mediation(df=df,p_threshold=mediation_p_threshold)
    readr::write_delim(mediation_table,file.path(outdir,paste0(model_names,"_data.txt")),delim = "\t")
    return(mediation_table)
  }, error = function(e) {
    message("Error in run_mediation: ", y, "\nMessage: ", conditionMessage(e))
    return(NULL)
  })
  
}
# mode="HMDB0062735_HMDB0240773_HMDB0000201"
# data<-differential_network1@Mediation_subnetwork@count_table
# metabolites<-data$feature_ID
# data_matrix <- as.matrix(t(data[,-1]))
# colnames(data_matrix)=metabolites
# if(nrow(resultdata)>0){
#   lapply(resultdata$Model, function(mode){
#     vector<-unlist(strsplit(mode, split = "_"))
#     x=vector[1]
#     m=vector[2]
#     y=vector[3]
#     dat <- tibble(
#       x = data_matrix[, x],
#       m = data_matrix[, m],
#       y = data_matrix[, y]
#     )
#     colnames(dat)=c(x,m,y)
#     png(paste0(mode,".png"),width = 2000,height = 2000,res=300)
#     p<-psych::mediate(as.formula(paste(y, "~", x, "+ (",m,")")), data = dat)
#     dev.off()
#   })
# }else{
#   message("No eligible mediation analysis models")
# }
# 
# readr::write_delim(resultdata,"Mediation_data.txt",delim = "\t")


add_mediation<-function(df,p_threshold = 0.05){
  # 添加中介状态列（英文）
  df <- df |>
    dplyr::mutate(
      mediation_status =  dplyr::case_when(
        # Condition 1: a and b are significant, c' not significant → Full Mediation
        (a_p <= p_threshold & b_p <= p_threshold) & c1_p > p_threshold ~ "Full Mediation",
        
        # Condition 2: a and b significant, c1 significant, ab and c1 same sign → Partial Mediation
        (a_p <= p_threshold & b_p <= p_threshold) &
          c1_p <= p_threshold &
          ((ab > 0 & c1 > 0) | (ab < 0 & c1 < 0)) ~ "Partial Mediation",
        
        # Condition 3: a and b significant, c1 significant, ab and c1 opposite signs → Suppression
        (a_p <= p_threshold & b_p <= p_threshold) &
          c1_p <= p_threshold &
          ((ab > 0 & c1 < 0) | (ab < 0 & c1 > 0)) ~ "Suppression",
        
        # Condition 4: a or b not significant, indirect CI includes 0 → No Significant Mediation
        (a_p > p_threshold | b_p > p_threshold) &
          indirect_lower <= 0 & indirect_upper >= 0 ~ "No Significant Mediation",
        
        # Condition 5: a or b not significant, indirect CI does NOT include 0, c1 not significant → Full Mediation
        (a_p > p_threshold | b_p > p_threshold) &
          ((indirect_lower > 0) | (indirect_upper < 0)) &
          c1_p > p_threshold ~ "Full Mediation",
        
        # Condition 6: a or b not significant, indirect CI does NOT include 0, c1 significant, same sign → Partial Mediation
        (a_p > p_threshold | b_p > p_threshold) &
          ((indirect_lower > 0) | (indirect_upper < 0)) &
          c1_p <= p_threshold &
          ((ab > 0 & c1 > 0) | (ab < 0 & c1 < 0)) ~ "Partial Mediation",
        
        # Condition 7: a or b not significant, indirect CI does NOT include 0, c1 significant, opposite sign → Suppression
        (a_p > p_threshold | b_p > p_threshold) &
          ((indirect_lower > 0) | (indirect_upper < 0)) &
          c1_p <= p_threshold &
          ((ab > 0 & c1 < 0) | (ab < 0 & c1 > 0)) ~ "Suppression",
        
        # Default
        TRUE ~ "No Significant Mediation"
      )
    )
  return(df)
}

