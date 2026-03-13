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
#tidyr dplyr tibble  purrr stringr
run_diff<-function(quantitative_table=NULL,samplelist=NULL,compare_group=NULL,p_threshold=0.05,FC_threshold=1.2,p_value_type="q_value",paired = FALSE){
  ###quantitative_table
  quantitative_table<-as.data.frame(quantitative_table) # 将定量表转换为数据框格式，确保数据格式统一
  # 检查定量表是否包含必需的"feature_ID"列
  if("feature_ID" %in% colnames(quantitative_table)){
    # 数据清洗和处理：
  core_table <- quantitative_table |>
    dplyr::distinct(!!dplyr::sym("feature_ID"), .keep_all = TRUE) |> # 去除重复的特征ID，每个特征只保留第一条记录
    tidyr::drop_na()  # 删除包含NA值的行，确保数据完整性
  }else{
    # 如果没有feature_ID列，停止执行并报错
    stop("No feature_ID column in the quantitative_table, no differential analysis will be performed.")
  }
 
  #compare_group:Experimental group:control group
  # 解析比较组信息：将"Experimental:Control"格式的字符串拆分为两组
  # 例如：如果compare_group = "Treatment:Control"，则cgroup = c("Treatment", "Control")
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
  
  ##samplelist
  # 验证样本信息表并准备比较分析
  if(all(c("group","sample") %in% colnames(samplelist))){
  compare_table <- tibble::as_tibble_col(cgroup) |>  # 将比较组转换为tibble列
    dplyr::group_by(.data$value) |>                    # 按组进行分组（Treatment和Control分别处理）
    dplyr::mutate(
      # 为每个组创建样本列表：支持多组合并（如"GroupA+GroupB"格式）
      sample_list = list(
        samplelist[samplelist$group %in% unlist(stringr::str_split(.data$value, "\\+")), ]$sample
        # stringr::str_split: 按"+"分割组名（支持多组合并，如"GroupA+GroupB"）
        # 获取属于这些组的所有样本名称
      )
    ) |>
    dplyr::mutate(sample_count = purrr::map_int(.data$sample_list, length))   # 计算每个组的样本数量
  }else{
    # 如果样本信息表缺少必要列，停止执行并报错
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
  # 情况1：配对t检验（需要样本数量相等且非零）
  if (paired & compare_table$sample_count[1] == compare_table$sample_count[2]) {
    if (compare_table$sample_count[1] != 0 | compare_table$sample_count[2] != 0) {
      res <- core_table |>
        dplyr::rowwise() |> 
        dplyr::mutate(
          p_value = t.test(
            x = dplyr::c_across(all_of(compare_table$sample_list[[1]])),  # 实验组数据
            y = dplyr::c_across(all_of(compare_table$sample_list[[2]])),  # 对照组数据
            paired = TRUE   # 配对t检验
          )$p.value,
          .keep = "none",  # 不保留其他列
          feature_ID = .data$feature_ID  # 保留特征ID
        ) |>
        dplyr::mutate(
          q_value = p.adjust(.data$p_value, method = "BH") # BH方法校正p值（FDR）
        )
    } else {
      stop("No sample in the group, please check the group information.")
    }
  } else {
    # 情况2：组内无样本
    if (compare_table$sample_count[1] == 0 | compare_table$sample_count[2] == 0) {
      stop("No sample in the group, please check the group information.")
    } else if (
      # 情况3：两组都只有1个样本（无法进行t检验）
      compare_table$sample_count[1] == 1 && compare_table$sample_count[2] == 1
    ) {
      message("Only FC")   # 只能计算fold change
      res <- core_table |> dplyr::mutate(
        p_value = NA, q_value = NA, .keep = "none",   # p值和q值设为NA
        feature_ID = .data$feature_ID
      )
      # 情况4：实验组1个样本，对照组多个样本（单样本t检验）
    } else if (
      compare_table$sample_count[1] == 1 && compare_table$sample_count[2] > 1
    ) {
      res <- core_table |>
        dplyr::rowwise() |>
        dplyr::mutate(
          p_value = t.test(
            x = dplyr::c_across(all_of(compare_table$sample_list[[2]])), # 对照组数据
            mu = .data[[compare_table$sample_list[[1]]]],     # 与实验组单个值比较
          )$p.value,
          .keep = "none",
          feature_ID = .data$feature_ID
        ) |>
        dplyr::mutate(
          q_value = p.adjust(.data$p_value, method = "BH")
        )
    } else if (compare_table$sample_count[1] > 1 && compare_table$sample_count[2] == 1) {
      # 情况5：实验组多个样本，对照组1个样本（单样本t检验）
      res <- core_table |>
        dplyr::rowwise() |>
        dplyr::mutate(
          p_value = t.test(
            x = dplyr::c_across(all_of(compare_table$sample_list[[1]])),# 实验组数据
            mu = .data[[compare_table$sample_list[[2]]]],# 与对照组单个值比较
          )$p.value,
          .keep = "none",
          feature_ID = .data$feature_ID
        ) |>
        dplyr::mutate(
          q_value = p.adjust(.data$p_value, method = "BH")
        )
    } else {
      # 情况6：标准的两独立样本t检验（两组都有多个样本）
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
  # 计算Fold Change和Log2FC
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

#####################################################################计算每组每个代谢物的归一化均值
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
  # 2. 数据标准化处理 ----
  mean_dat<-as.data.frame(core_table)
  rownames(mean_dat)<-mean_dat$feature_ID  
  mean_dat<-mean_dat |>
    dplyr::select(-1) |>  # 去掉第一列（如化合物ID等）
    t() |>   # 转置
    as.data.frame(rownames.force = TRUE,colnames.force = TRUE)   # 强制保留行名为数据框的行名
  # 3. 最大值标准化：将每个特征的值除以其最大值，缩放到0-1范围
  mean_table<- apply(mean_dat,2,function(x){
    #x[x==0.00000000000001]<-NA #如果count_table做了fix缺失值填充
    #scales::rescale(x, to = c(0, 1))
    if(max(x)!=0){
      x/max(x)
    }else{
      x # 如果最大值为0，保持原值（通常都是0）
    }
    
    })
  # 4. 数据重构 ----
  mean_table<-mean_table |>
    t() |> # 转置回原始方向（特征为行，样本为列）
    as.data.frame(rownames.force = TRUE,colnames.force = TRUE)
  mean_table$feature_ID<-rownames(mean_table) # 添加feature_ID列
  mean_table<-mean_table |>
    dplyr::select(feature_ID, everything())     # 把 feature_ID 放在第一列
  
  #compare_group:Experimental group:control group
  # 5. 比较组解析 ----
  # 将"Experimental:Control"格式的字符串拆分为两组
  cgroup <-unlist(strsplit(x = compare_group, split = ":"))
  
  ##samplelist
  # 6. 样本信息验证和处理 ----
  if(all(c("group","sample") %in% colnames(samplelist))){
    # 创建比较分析表
    compare_table <- tibble::as_tibble_col(cgroup) |>
      dplyr::group_by(.data$value) |>  # 按组分组
      dplyr::mutate(
        # 为每个组创建样本列表（支持多组合并，如"GroupA+GroupB"）
        sample_list = list(
          samplelist[samplelist$group %in% unlist(stringr::str_split(.data$value, "\\+")), ]$sample
        )
      ) |>
      # 计算每组样本数量
      dplyr::mutate(sample_count = purrr::map_int(.data$sample_list, length)) 
  }else{
    stop("No group column or sample column in the samplelist, no differential analysis will be performed.")
  }
  
  # 7. 计算分组标准化平均值 ----
  result_mean<-mean_table |>
    dplyr::mutate(
      # 计算实验组的标准化平均值
      case_norm_mean=rowMeans(pick(compare_table$sample_list[[1]])) ,
      # 计算对照组的标准化平均值
      control_norm_mean=rowMeans(pick(compare_table$sample_list[[2]]))
      )  |>
    dplyr::select(feature_ID,case_norm_mean,control_norm_mean)  # 只保留需要的列
  return(result_mean)
}
#check
# mean_table=as.data.frame(mean_table)
# mean_table$sample=rownames(mean_dat)
# mean_table<-mean_table |>
#   dplyr::left_join(samplelist,by=c("sample"="sample"))
# ggplot(mean_table, aes(x = 1, y = HMDB0000138, color = group)) +
#   geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
#   geom_jitter(width = 0.1, size = 3) +
#   labs(x = "", y = "HMDB0000138 Value") +
#   theme_minimal() + 
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())

# 处理函数：将字符串转换为数值向量，计算归一化后的平均值
process_intensity <- function(intensity1_str, intensity2_str) {
  # 将分号分隔的字符串转换为数值向量
  # str_split: 将字符串按分号分割
  # simplify = TRUE: 返回矩阵而不是列表
  # as.numeric: 转换为数值，无法转换的会变成NA
  # suppressWarnings: 抑制转换过程中产生的警告（如无法转换的字符）
  vec1 <- suppressWarnings(as.numeric(str_split(intensity1_str, ";", simplify = TRUE)))
  vec2 <- suppressWarnings(as.numeric(str_split(intensity2_str, ";", simplify = TRUE)))
  
  # 找到总和的最大值（忽略NA）
  max_total <- max(c(vec1,vec2), na.rm = TRUE)
  # 找到两个向量中所有数值的最大值（忽略NA）
  # 用于后续的归一化处理
  
  # 归一化处理
  # 归一化处理：将每个值除以最大值，缩放到0-1范围
  if (is.finite(max_total) && max_total > 0) {
    norm1 <- vec1 / max_total
    norm2 <- vec2 / max_total
  } else {
    # 如果最大值无效（全NA或零），则返回NA
    norm1 <- rep(NA, length(vec1))
    norm2 <- rep(NA, length(vec2))
  }
  
  # 计算归一化后的平均值（忽略NA）
  mean1 <- mean(norm1, na.rm = TRUE)
  mean2 <- mean(norm2, na.rm = TRUE)
  
  # 处理全NA的情况
  # 处理全NA的情况：如果mean返回NaN（所有值都是NA），转换为NA
  if (is.nan(mean1)) mean1 <- NA
  if (is.nan(mean2)) mean2 <- NA
  
  return(list(mean_intensity1 = mean1, mean_intensity2 = mean2))
}

# 包装函数：对差异分析表中的每一行应用处理函数
run_process_intensity<-function(diff_table){
  # 应用处理函数到每一行
  diff_table <- diff_table %>%
    rowwise() %>%  # 按行操作
    mutate(
      # 对每行应用process_intensity函数
      result = list(process_intensity(`Intensity(1)`, `Intensity(2)`)),
      # 提取处理结果：通常Intensity(2)是实验组(case)
      case_norm_mean = result$mean_intensity2,
      # Intensity(1)是对照组(control)
      control_norm_mean = result$mean_intensity1
    ) %>%
    select(-result)  # 移除临时结果列 
}


