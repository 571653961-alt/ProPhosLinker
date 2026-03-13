run_samplelist <- function(input_str) {
  # 1. 分割成向量
  comparisons <- unlist(strsplit(input_str, ";"))
  
  # 2. 创建空列表存储结果
  samplelist <- list()
  
  # 3. 处理每个比较组
  for (comp in comparisons) {
    # 检查格式是否正确
    if (!grepl("@", comp) || !grepl("-VS-", comp)) {
      warning(paste("跳过格式不正确的比较组:", comp))
      next
    }
    
    # 分割组名和样本信息
    parts <- unlist(strsplit(comp, "@"))
    
    # 获取组名部分 (如 "PCOS-VS-Control")
    group_str <- parts[1]
    
    # 分割组名
    group_parts <- unlist(strsplit(group_str, "-VS-"))
    group1 <- group_parts[1]
    group2 <- group_parts[2]
    
    # 处理样本部分 (如 "(PCOS_1,...)/(Control_1,...)")
    sample_str <- parts[2]
    
    # 检查样本字符串格式
    if (!grepl("^\\(", sample_str) || !grepl("\\)$", sample_str)) {
      warning(paste("跳过格式不正确的样本字符串:", sample_str))
      next
    }
    
    # 移除开头和结尾的括号
    sample_str <- gsub("^\\(", "", sample_str)
    sample_str <- gsub("\\)$", "", sample_str)
    
    # 分割样本字符串
    samples <- unlist(strsplit(sample_str, "/"))
    
    # 检查是否正好有两组样本
    if (length(samples) != 2) {
      warning(paste("样本格式不正确，应有两组样本:", sample_str))
      next
    }
    
    # 清理样本名称 - 移除所有括号和空格
    clean_samples <- function(sample_list) {
      # 分割样本
      samples <- unlist(strsplit(sample_list, ","))
      # 移除括号和空格
      samples <- gsub("[\\(\\)\\s]", "", samples)
      # 移除空值
      samples <- samples[samples != ""]
      return(samples)
    }
    
    # 提取并清理组1的样本
    group1_samples <- clean_samples(samples[1])
    
    # 提取并清理组2的样本
    group2_samples <- clean_samples(samples[2])
    
    # 检查样本是否为空
    if (length(group1_samples) == 0 || length(group2_samples) == 0) {
      warning(paste("一组或多组样本为空:", sample_str))
      next
    }
    
    # 创建当前比较组的数据框
    df_group1 <- data.frame(
      group = group1,
      sample = group1_samples,
      comparison = group_str,
      stringsAsFactors = FALSE
    )
    
    df_group2 <- data.frame(
      group = group2,
      sample = group2_samples,
      comparison = group_str,
      stringsAsFactors = FALSE
    )
    
    # 合并两个组
    samplelist[[group_str]] <- rbind(df_group1, df_group2)
  }
  
  # 4. 检查是否有有效数据
  if (length(samplelist) == 0) {
    stop("没有有效的比较组数据")
  }
  
  # 5. 合并所有比较组的结果
  final_samplelist <- do.call(rbind, samplelist)
  
  # 6. 重置行名
  rownames(final_samplelist) <- NULL
  
  return(final_samplelist)
}
