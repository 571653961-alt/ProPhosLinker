#' Perform a simple mediation analysis using three specified features (X → M → Y).
#'
#' This function conducts a mediation analysis based on a user-specified independent variable (`x`),
#' mediator (`m`), and dependent variable (`y`) from a count/abundance table, using the `psych::mediate`
#' function. It validates input data alignment with sample metadata, fits the mediation model,
#' generates a high-resolution path diagram (PNG), computes key statistics (a, b, c, c', indirect effect),
#' classifies the mediation type via bootstrapped confidence intervals and p-values,
#' and saves results to disk.
#'
#' @param count_table A data frame with rows as features (must include column `feature_ID`)
#'                    and columns as samples. Numeric values represent feature abundances.
#' @param samplelist A data frame with columns `sample` (sample IDs) and `group` (group labels).
#' @param x Character string: name of the independent variable (predictor) in `feature_ID`.
#' @param y Character string: name of the dependent variable (outcome) in `feature_ID`.
#' @param m Character string: name of the mediator variable in `feature_ID`.
#' @param group_name Character vector specifying which group(s) to subset from `samplelist`.
#'                   If length > 1, collapsed with "-vs-"; if contains ":", replaced with "-vs-".
#' @param mediation_p_threshold Numeric: significance threshold for p-values in mediation inference (default: 0.05).
#' @param outdir Output directory path to save the mediation plot and result table (default: `"./"`).
#'
#' @return A data frame with mediation statistics and classification (e.g., "Full Mediation", "Partial Mediation"),
#'         or `NULL` if an error occurs (with error message printed).
#'         The function also produces:
#'         \itemize{
#'           \item A PNG file: `{mediation_x_m_y_group}.png` (2000×2000, 300 DPI)
#'           \item A tab-delimited text file: `{mediation_x_m_y_group}_data.txt`
#'         }
#'
#' @importFrom dplyr filter select distinct mutate case_when sym all_of
#' @importFrom readr write_delim
#' @importFrom graphics png dev.off
#' @importFrom utils file.path
#' @import psych
#' @import dplyr
#'
#' @examples
#' # Example setup
#' count_table <- data.frame(
#'   feature_ID = c("X", "M", "Y", "Z"),
#'   S1 = c(1, 2, 3, 4),
#'   S2 = c(2, 3, 4, 5),
#'   S3 = c(1.5, 2.5, 3.5, 4.5)
#' )
#' samplelist <- data.frame(
#'   sample = c("S1", "S2", "S3"),
#'   group = "Control"
#' )
#' result <- run_mediation(count_table, samplelist, x = "X", m = "M", y = "Y", group_name = "Control")
#'
#' @export
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

    data_matrix <- as.matrix(t(data))
    colnames(data_matrix)<-count_table$feature_ID

    f <- paste(y, "~", x, "+ (",m,")")
    model_names <-paste("mediation",x,  m, y,group_name, sep = "_")
    png(file.path(outdir,paste0(model_names,".png")),width = 2000,height = 2000,res=300)
    model <-psych::mediate(as.formula(f), data = data_matrix)
    dev.off()

        df<- data.frame(
          Model =paste(x,m, y, sep = "_") ,
          a=round(as.numeric(model[["a"]][1]),2),
          b=round(as.numeric(model[["b"]][1]),2),
          ab=round(as.numeric(model[["ab"]][1]),2),
          c=round(as.numeric(model[["c"]][1]),2),
          c1= round(as.numeric(model[["indirect"]][1]),2),
          model_r2=round(as.numeric(model[["cprime.reg"]][["R2"]]),2),
          c_p   = round(as.numeric(model[["total.reg"]][["prob"]][2]),4),
          c1_p  = round(as.numeric(model[["cprime.reg"]][["prob"]][2]),4),    
          a_p=  round(as.numeric(model[["a.reg"]][["prob"]][2]),4),
          b_p=  round(as.numeric(model[["b.reg"]][["prob"]][2]),4),
          indirect_lower   = round(as.numeric(model[["boot"]][["ci.ab"]][1]),2),  
          indirect_upper   =round(as.numeric(model[["boot"]][["ci.ab"]][2]),2)
        )
    mediation_table<-add_mediation(df=df,p_threshold=mediation_p_threshold)
    readr::write_delim(mediation_table,file.path(outdir,paste0(model_names,"_data.txt")),delim = "\t")
    return(mediation_table)
  }, error = function(e) {
    message("Error in run_mediation: ", y, "\nMessage: ", conditionMessage(e))
    return(NULL)
  })
  
}


#' Classify mediation effects based on statistical criteria and bootstrapped confidence intervals.
#'
#' This helper function takes a data frame of mediation statistics (from `psych::mediate`)
#' and assigns a mediation type label using conventional rules based on:
#' (1) significance of paths a, b, and c';
#' (2) sign consistency between indirect (ab) and direct (c') effects;
#' (3) whether the bootstrapped CI for the indirect effect excludes zero.
#'
#' Classification includes: "Full Mediation", "Partial Mediation", "Suppression", or "No Significant Mediation".
#'
#' @param df A data frame containing columns:
#'   `a_p`, `b_p`, `c1_p`, `ab`, `c1`, `indirect_lower`, `indirect_upper`.
#' @param p_threshold Significance threshold for p-values (default: 0.05).
#'
#' @return The input `df` with an additional column `mediation_status`.
#'
#' @importFrom dplyr mutate case_when
#' @import dplyr
#'
#' @keywords internal
add_mediation<-function(df,p_threshold = 0.05){
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