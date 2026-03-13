# functional_enrichment_analysis

# suppressWarnings(library(optparse))
# option_list <- list(
#   # Input parameters
#   make_option(c("--pro_diff"), type = "character", metavar = "FILE",
#               help = "Differential protein expression file (CSV)"),
#   make_option(c("--phos_diff"), type = "character", metavar = "FILE",
#               help = "Differential phosphoprotein expression file (CSV)"),
#   make_option(c("--phos_pro"), type = "character", metavar = "FILE",
#               help = "Phosphosite to protein mapping file (CSV)"),
# 
#   # Output parameters
#   make_option(c("--outdir"), type = "character", default = getwd(), metavar = "DIR",
#               help = "Output directory [default: %default]"),
# 
#   # Analysis parameters
#   make_option(c("--log2FC"), type = "numeric", default = 1.2, metavar = "NUM",
#               help = "Log2 fold change cutoff for differential expression [default: %default]"),
#   make_option(c("--diff_p_adj"), type = "numeric", default = 0.05, metavar = "NUM",
#               help = "Adjusted p-value cutoff for significance [default: %default]"),
# 
# 
#   # Omics names
#   make_option(c("--omics_name1"), type = "character", default = "Proteomics", metavar = "STR",
#               help = "Name for first omics dataset [default: %default]"),
#   make_option(c("--omics_name2"), type = "character", default = "Phosphoproteomics", metavar = "STR",
#               help = "Name for second omics dataset [default: %default]"),
# 
# 
#   make_option(c("--GO_showCategory"), type = "integer", default = 6, metavar = "INT",
#               help = "Number of GO terms to display [default: %default]"),
#   make_option(c("--KEGG_showCategory"), type = "integer", default = 15, metavar = "INT",
#               help = "Number of KEGG pathways to display [default: %default]"),
#   make_option(c("--pvalueCutoff"), type = "numeric", default = 0.05, metavar = "NUM",
#               help = "P-value cutoff for enrichment analysis [default: %default]"),
# 
#   #general parameter
#   make_option(c("-v", "--verbose"), action = "store_true", default = FALSE, metavar = "FLAG",
#               help = "Print detailed output messages [Default: %default]")
# )
# 
# # Parse command line arguments
# opt_parser <- OptionParser(option_list = option_list,
#                            description = "Perform comparative enrichment analysis (GO and KEGG) between protein and phosphoprotein datasets.")
# opt <- parse_args(opt_parser)
# 
# #parameter validation
# parameter_validation <- function(opt){
#   #path cleaning
#   clean_path <- function(path) {
#     # Trim leading and trailing whitespace and quotes
#     path <- gsub("^['\" ]+|['\" ]+$","",path)
#     # Replace backslashes (\) with forward slashes (/)
#     path <- gsub("\\\\","/",path)
#     return(path)
#   }
#   check_input_file <- function(file_path,file_name){
#     if(!file.exists(file_path)){
#       stop("❌ Error: The file '", file_name, "' does not exist at the specified path: ", file_path)
#     }
#     if(file.size(file_path) == 0){
#       stop("❌ Error: The file '", file_name, "' is empty at the specified path: ", file_path)
#     }
#   }
#   check_output_dir <- function(result_dir){
#     if (!dir.exists(result_dir)) {
#       # Recursively create directory and all parent directories
#       dir.create(result_dir, recursive = TRUE)
#       print(paste("✅ Output directory created:", result_dir))
#     } else {
#       print(paste("✅ Output directory already exists:", result_dir))
#     }
# 
#     return(result_dir)
#   }
#   validate_option_choises <- function(value, allowed_values, option_name) {
#     if (!value %in% allowed_values) {
#       stop("❌ Error: Invalid", option_name, "value '", value, "'. Please use one of the following options: ",
#            paste(allowed_values, collapse = ", "))
#     }
#   }
#   validate_numeric_range <- function(value, min_val, max_val, option_name) {
#     if (is.null(value)) {
#       return()  # Skip validation if the value is NULL
#     }
#     if (value < min_val || value > max_val) {
#       stop("❌ Error: ", option_name, " must be within the range [", min_val, ", ", max_val, "].")
#     }
#   }
# 
#   # clean all pathway parameter
#   opt$pro_diff <- clean_path(opt$pro_diff)
#   opt$phos_diff <- clean_path(opt$phos_diff)
#   opt$phos_pro <- clean_path(opt$phos_pro)
#   # opt$module_info <- clean_path(opt$module_info)
#   opt$outdir <- clean_path(opt$outdir)
# 
#   #check input files
#   check_input_file(opt$pro_diff,"pro_diff")
#   check_input_file(opt$phos_diff,"phos_diff")
#   check_input_file(opt$phos_pro,"phos_pro")
#   opt$outdir  <- check_output_dir(opt$outdir)
# 
# 
#   #choices validation
#   validate_numeric_range(opt$log2FC, 0,10, "log2FC")
#   validate_numeric_range(opt$diff_p_adj, 0,1, "diff_p_adj")
#   validate_numeric_range(opt$GO_showCategory, 1,20, "GO_showCategory")
#   validate_numeric_range(opt$KEGG_showCategory, 1,20, "KEGG_showCategory")
#   validate_numeric_range(opt$pvalueCutoff, 0,1, "pvalueCutoff")
# 
#   if (opt$verbose) { print(opt)}
# 
#   return(opt)
# }



suppressWarnings(suppressPackageStartupMessages(library(clusterProfiler)))
suppressWarnings(suppressPackageStartupMessages(library(org.Hs.eg.db)))
suppressWarnings(suppressPackageStartupMessages(library(enrichplot)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
options(timeout = 1000)  # 默认是60秒


###======================================================parameters

# 0.
enrichment_predata <- function(pro_diff_path,phos_diff_path,phos_pro_path,outdir="./",log2FC = 1.2,p_adj=0.05){
  
  pro_diff <- read.csv(pro_diff_path,sep = '\t')
  phos_diff <- read.csv(phos_diff_path,sep = '\t')
  phos_pro <- read.csv(phos_pro_path,sep = '\t')
  
  pro_diff_filter_list <- pro_diff[((pro_diff$logFC >= log2FC | pro_diff$logFC <= - log2FC) & pro_diff$adj.P.Val < p_adj),]$Protein
  phos_diff_filter_list <- phos_diff[((phos_diff$logFC >= log2FC | phos_diff$logFC <= - log2FC) & phos_diff$adj.P.Val < p_adj),]$Protein
  
  phos_diff_filter_list <- phos_pro[phos_pro$Site  %in% phos_diff_filter_list,]$Protein
  
  return(list("protein_list"=pro_diff_filter_list,"phosphoprotein_list" =phos_diff_filter_list))
}

# 1.单组学——GO富集分析：
single_omics_GO_enrichment <- function(omics_list,omics_name,outdir ="./",pvalueCutoff = 0.05,showCategory = 6){
diff_pro_go_all<- enrichGO(gene          = omics_list,
                           OrgDb         = org.Hs.eg.db, # 指定物种数据库
                           keyType       = 'SYMBOL',     # 输入ID类型
                           ont           = "ALL",         # 本体： BP, MF, CC 或 "ALL"
                           pAdjustMethod = "BH",         # p值校正方法
                           pvalueCutoff  = pvalueCutoff,         # p值阈值
                           qvalueCutoff  = 0.2,          # q值阈值
                           readable      = FALSE)

if(!is.null(diff_pro_go_all)){
# 提取富集结果数据
go_data <- as.data.frame(diff_pro_go_all)
write.table(go_data, 
            file = paste0(outdir, "/", omics_name, "_GO_enrichment.tsv"), 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

# 数据预处理
go_data_clean <- go_data %>%
  # 计算GeneRatio的数值
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
    BgRatio_num = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
    FoldEnrichment = GeneRatio_num / BgRatio_num
  ) %>%
  # 按本体和p值排序
  group_by(ONTOLOGY) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  # 每个本体取前6个
  slice_head(n = showCategory) %>%
  ungroup() %>%
  # 创建排序因子，确保正确的显示顺序
  mutate(
    Description = factor(Description, levels = rev(unique(Description))),
    ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC"))
  )%>%
  arrange(desc(GeneRatio_num), .by_group = TRUE)%>%
  # 创建排序因子，确保正确的显示顺序
  mutate(
    Description = factor(Description, levels = rev(unique(Description))),
    ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC")),
    # 如果需要将GeneRatio_num转换为因子，在这里进行
    GeneRatio_factor = as.factor(GeneRatio_num)  # 可选：如果需要因子形式
  )

# 创建ggplot2图形
p_GO <- ggplot(go_data_clean, aes(x = GeneRatio_num, y = Description)) +
  # 添加点，大小表示Count，颜色表示p.adjust
  geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
  # 分面显示三个本体
  facet_grid(ONTOLOGY ~ ., 
             scales = "free_y", 
             space = "free_y"
             ) +
  
  # 颜色渐变设置
  scale_color_gradient(low = "#175663", high = "#90362d", 
                       name = "Adjusted P-value",
                       trans = "log10") +
  
  # 点的大小范围设置
  scale_size_continuous(range = c(2, 8), 
                        name = "Protein Count",
                        guide = guide_legend(
                          override.aes = list(color = "grey60")  # 图例圆点改为灰色
                        )) +
  
  # 坐标轴和标题
  labs(x = "Gene Ratio", 
       y = "",
       title = paste0("GO Enrichment Analysis of ",omics_name)) +
  
  # 主题设置
  theme_bw() +
  theme(
    # 标题设置
    plot.title = element_text(hjust = 0.5, face = "bold", size = 50,
                              margin = margin(b = 20)),
    
    # 分面标题设置
    strip.text = element_text(size = 35, face = "bold",
                              margin = margin(t = 5, r = 5, b = 5, l = 5)),
    strip.background = element_rect(fill = "lightgray"),
    
    # 坐标轴设置
    axis.text.x = element_text(size = 30, color = "grey20", face = "bold"),
    axis.text.y = element_text(size = 30, color = "grey20", face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold",
                                margin = margin(t = 15)),
    
    # 图例设置
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 25),
    legend.key.size = unit(0.8, "cm"),
    
    # 网格线和边距
    panel.grid.major = element_line(color = "grey90", linewidth  = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )+
  guides(
    fill = guide_colorbar(barheight = unit(4, "cm")),
    shape = guide_legend(
      override.aes = list(
        fill = "grey80",
        # color = "black",    # 关键：覆盖透明边框
        size = 6,
        stroke = 0.5
      ))
  )

ggsave(paste0(outdir,"/",omics_name,"_GO_enrichment.png"),p_GO,width = 12, height =8, dpi = 300)
}
}

# 1.单组学—— KEGG富集分析：
single_omics_KEGG_enrichment <- function(omics_entrez_list,omics_name,outdir ="./",pvalueCutoff = 0.05,showCategory = 15){
  # 现在用 Entrez ID 进行 KEGG 富集分析
  kk <- enrichKEGG(
    gene = omics_entrez_list,
    organism = "hsa",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = 0.2
  )
  # print(kk)
  # print(is.null(kk))
  if(!is.null(kk)){
  # 提取富集结果数据
  kk_data <- as.data.frame(kk)
  write.table(kk_data, 
              file = paste0(outdir, "/", omics_name, "_KEGG_enrichment.tsv"), 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE, 
              quote = FALSE)
  # print(kk_data)
  # print(nrow(kk_data))
  if(nrow(kk_data) > 0){
  
  # 数据预处理
  kk_data_clean <- kk_data %>%
    # 计算GeneRatio的数值
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
      BgRatio_num = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
      FoldEnrichment = GeneRatio_num / BgRatio_num
    ) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    # 每个本体取前n个
    slice_head(n = showCategory) %>%
    arrange(desc(Count), .by_group = TRUE)%>%
    # # 创建排序因子，确保正确的显示顺序
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      GeneRatio_factor = as.factor(GeneRatio_num)  # 可选：如果需要因子形式
    )
    
  p_KEGG <- ggplot(kk_data_clean, aes(x = GeneRatio_num, y = Description)) +
    # 添加点，大小表示Count，颜色表示p.adjust
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
    # 颜色渐变设置
    scale_color_gradient(low = "#175663", high = "#90362d", 
                         name = "Adjusted P-value",
                         trans = "log10") +
    
    # 点的大小范围设置
    scale_size_continuous(range = c(2, 8), 
                          name = "Protein Count",
                          guide = guide_legend(
                            override.aes = list(color = "grey60")  # 图例圆点改为灰色
                          )) +
    
    # 坐标轴和标题
    labs(x = "Gene Ratio", 
         y = "",
         title = paste0("KEGG Enrichment Analysis of ",omics_name)) +
    
    # 主题设置
    theme_bw() +
    theme(
      # 标题设置
      plot.title = element_text(hjust = 0.5, face = "bold", size = 50,
                                margin = margin(b = 20)),
      # 坐标轴设置
      axis.text.x = element_text(size = 30, color = "grey20", face = "bold"),
      axis.text.y = element_text(size = 30, color = "grey20", face = "bold"),
      axis.title.x = element_text(size = 30, face = "bold",
                                  margin = margin(t = 15)),
      
      # 图例设置
      legend.title = element_text(size = 30, face = "bold"),
      legend.text = element_text(size = 25),
      legend.key.size = unit(0.8, "cm"),
      
      # 网格线和边距
      panel.grid.major = element_line(color = "grey90", linewidth  = 0.5),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )+
    guides(
      fill = guide_colorbar(barheight = unit(4, "cm")),
      shape = guide_legend(
        override.aes = list(
          fill = "grey80",
          # color = "black",    # 关键：覆盖透明边框
          size = 6,
          stroke = 0.5
        ))
    )
  ggsave(paste0(outdir,"/",omics_name,"_KEGG_enrichment.png"),p_KEGG,width = 10, height =6, dpi = 300)
  }
  }
  
}

# 2.创建对比分析的数据框
two_omicses_GO_enrichment <- function(omics1_list,omics2_list,omics_name1,omics_name2,outdir ="./",pvalueCutoff = 0.05,showCategory = 15){
  compare_lists <- list(Protein = omics1_list,
                   Phosphoprotein = omics2_list)

# 对比GO富集分析
go_compare <- compareCluster(geneClusters = compare_lists,
                             fun = "enrichGO",
                             OrgDb = org.Hs.eg.db,
                             keyType = "SYMBOL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             pvalueCutoff = pvalueCutoff,
                             qvalueCutoff = 0.2)

if(!is.null(go_compare)){

# 提取富集结果数据
go_compare_data <- as.data.frame(go_compare)


write.table(go_compare_data, 
            file = paste0(outdir,"/",omics_name1,"_",omics_name2,"_GO_enrichment.tsv"), 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

# 数据预处理
go_compare_clean <- go_compare_data %>%
  # 计算GeneRatio的数值
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
    BgRatio_num = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
    FoldEnrichment = GeneRatio_num / BgRatio_num
  ) %>%
  # 按本体和p值排序
  group_by(ONTOLOGY) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  # 每个本体取前showCategory个
  slice_head(n = showCategory) %>%
  ungroup() %>%
  # 创建排序因子，确保正确的显示顺序
  mutate(
    Description = factor(Description, levels = rev(unique(Description))),
    ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC"))
  )%>%
  arrange(desc(GeneRatio_num), .by_group = TRUE)%>%
  # 创建排序因子，确保正确的显示顺序
  mutate(
    Description = factor(Description, levels = rev(unique(Description))),
    ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC")),
    # 如果需要将GeneRatio_num转换为因子，在这里进行
    GeneRatio_factor = as.factor(GeneRatio_num)  # 可选：如果需要因子形式
  )


# 创建ggplot2图形
go_compare_GO <- ggplot(go_compare_clean, aes(x = Cluster, y = Description)) +
  # 添加点，大小表示Count，颜色表示p.adjust
  geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
  # 分面显示三个本体
  facet_grid(ONTOLOGY ~ .,
             scales = "free_y",
             space = "free_y"
  ) +
  # 颜色渐变设置
  scale_color_gradient(low = "#175663", high = "#90362d",
                       name = "Adjusted P-value",
                       trans = "log10") +
  # 点的大小范围设置
  scale_size_continuous(range = c(2, 8),
                        name = "Protein Count",
                        guide = guide_legend(
                          override.aes = list(color = "grey60")  # 图例圆点改为灰色
                        )) +

  # 坐标轴和标题
  labs(x = "Gene Ratio",
       y = "",
       title = paste0("GO Enrichment Comparison: ",omics_name1,"  vs ", omics_name2)) +

  # 主题设置
  theme_bw() +
  theme(
    # 标题设置
    plot.title = element_text(hjust = 0.5, face = "bold", size = 50,
                              margin = margin(b = 20)),

    # 分面标题设置
    strip.text = element_text(size = 35, face = "bold",
                              margin = margin(t = 5, r = 5, b = 5, l = 5)),
    strip.background = element_rect(fill = "lightgray"),

    # 坐标轴设置
    axis.text.x = element_text(size = 30, color = "grey20", face = "bold"),
    axis.text.y = element_text(size = 30, color = "grey20", face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold",
                                margin = margin(t = 15)),

    # 图例设置
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 25),
    legend.key.size = unit(0.8, "cm"),

    # 网格线和边距
    panel.grid.major = element_line(color = "grey90", linewidth  = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )+
  guides(
    fill = guide_colorbar(barheight = unit(4, "cm")),
    shape = guide_legend(
      override.aes = list(
        fill = "grey80",
        # color = "black",    # 关键：覆盖透明边框
        size = 6,
        stroke = 0.5
      ))
  )

ggsave(paste0(outdir,"/",omics_name1,"_",omics_name2,"_GO_enrichment.png"),go_compare_GO,width = 12, height =8, dpi = 300)
}
}

two_omicses_KEGG_enrichment <- function(omics1_list,omics2_list,omics_name1,omics_name2,outdir ="./",pvalueCutoff = 0.05, showCategory = 15){
  # omics1_list_entrez <- bitr(omics1_list, 
  #                            fromType = "SYMBOL",      # 输入ID类型
  #                            toType = "ENTREZID",      # 输出ID类型
  #                            OrgDb = org.Hs.eg.db)$ENTREZID     # 物种数据库
  # 
  # omics2_list_entrez <- bitr(omics2_list, 
  #                            fromType = "SYMBOL",      # 输入ID类型
  #                            toType = "ENTREZID",      # 输出ID类型
  #                            OrgDb = org.Hs.eg.db)$ENTREZID     # 物种数据库
  compare_lists <- list(Protein = omics1_list,
                        Phosphoprotein = omics2_list)
  
  kk_compare <- compareCluster(geneClusters = compare_lists,
                               fun="enrichKEGG",
                               organism="hsa",
                               pvalueCutoff = pvalueCutoff)

  if(!is.null(kk_compare)){
  # 提取富集结果数据
  kk_compare_data <- as.data.frame(kk_compare)
  if(nrow(kk_compare_data) > 0){
  write.table(kk_compare_data, 
              file = paste0(outdir,"/",omics_name1,"_",omics_name2,"_KEGG_enrichment.tsv"), 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE, 
              quote = FALSE)
  
  # 数据预处理
  kk_compare_data_clean <- kk_compare_data %>%
    # 计算GeneRatio的数值
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
      BgRatio_num = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
      FoldEnrichment = GeneRatio_num / BgRatio_num
    ) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    # 每个本体取前n个
    slice_head(n = showCategory) %>%
    arrange(desc(Count), .by_group = TRUE)%>%
    # # 创建排序因子，确保正确的显示顺序
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      GeneRatio_factor = as.factor(GeneRatio_num)  # 可选：如果需要因子形式
    )
  
  p_kk_compare <- ggplot(kk_compare_data_clean, aes(x = Cluster, y = Description)) +
    # 添加点，大小表示Count，颜色表示p.adjust
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
    # 颜色渐变设置
    scale_color_gradient(low = "#175663", high = "#90362d", 
                         name = "Adjusted P-value",
                         trans = "log10") +
    
    # 点的大小范围设置
    scale_size_continuous(range = c(2, 8), 
                          name = "Protein Count",
                          guide = guide_legend(
                            override.aes = list(color = "grey60")  # 图例圆点改为灰色
                          )) +
    
    # 坐标轴和标题
    labs(x = "Gene Ratio", 
         y = "",
         title = paste0("KEGG Enrichment Comparison: ",omics_name1,"  vs ", omics_name2)) +
    
    # 主题设置
    theme_bw() +
    theme(
      # 标题设置
      plot.title = element_text(hjust = 0.5, face = "bold", size = 50,
                                margin = margin(b = 20)),
      # 坐标轴设置
      axis.text.x = element_text(size = 30, color = "grey20", face = "bold"),
      axis.text.y = element_text(size = 30, color = "grey20", face = "bold"),
      axis.title.x = element_text(size = 30, face = "bold",
                                  margin = margin(t = 15)),
      
      # 图例设置
      legend.title = element_text(size = 30, face = "bold"),
      legend.text = element_text(size = 25),
      legend.key.size = unit(0.8, "cm"),
      
      # 网格线和边距
      panel.grid.major = element_line(color = "grey90", linewidth  = 0.5),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )+
    guides(
      fill = guide_colorbar(barheight = unit(4, "cm")),
      shape = guide_legend(
        override.aes = list(
          fill = "grey80",
          # color = "black",    # 关键：覆盖透明边框
          size = 6,
          stroke = 0.5
        ))
    )
  
  ggsave(paste0(outdir,"/",omics_name1,"_",omics_name2,"_KEGG_enrichment.png"),p_kk_compare,width = 12, height =8, dpi = 300)
  }
}
}
# 3.
omics_enrichment <- function(pro_diff_path,phos_diff_path,phos_pro_path,outdir ="./",log2FC=1.2,diff_p_adj=0.05,omics_name1="Proteomics",omics_name2="Phosphoproteomics",pvalueCutoff = 0.05,GO_showCategory=6,KEGG_showCategory=15){
  
  predata <- enrichment_predata(pro_diff_path,phos_diff_path,phos_pro_path,outdir,log2FC,diff_p_adj)
  omics1_list <- predata$protein_list
  omics2_list <- predata$phosphoprotein_list
  
 
  # 将 Gene Symbol 转换为 Entrez ID
  omics1_entrez_list <- bitr(omics1_list, 
                      fromType = "SYMBOL",      # 输入ID类型
                      toType = "ENTREZID",      # 输出ID类型
                      OrgDb = org.Hs.eg.db)$ENTREZID     # 物种数据库
  
  omics2_entrez_list <- bitr(omics2_list, 
                      fromType = "SYMBOL",      # 输入ID类型
                      toType = "ENTREZID",      # 输出ID类型
                      OrgDb = org.Hs.eg.db)$ENTREZID     # 物种数据库
  ###############改为不做单组学，只做两组学比较###############
  # # 单组学GO富集分析
  # single_omics_GO_enrichment(omics1_list,omics_name1,outdir,pvalueCutoff,GO_showCategory)
  # single_omics_GO_enrichment(omics2_list,omics_name2,outdir,pvalueCutoff,GO_showCategory)
  # # 单组学KEGG富集分析
  # single_omics_KEGG_enrichment(omics1_entrez_list,omics_name1,outdir,pvalueCutoff,KEGG_showCategory)
  # single_omics_KEGG_enrichment(omics2_entrez_list,omics_name2,outdir,pvalueCutoff,KEGG_showCategory)
  
  # 两组学比较GO、KEGG富集分析
  two_omicses_GO_enrichment(omics1_list,omics2_list,omics_name1,omics_name2,outdir,pvalueCutoff,GO_showCategory)
  two_omicses_KEGG_enrichment(omics1_entrez_list,omics2_entrez_list,omics_name1,omics_name2,outdir,pvalueCutoff,KEGG_showCategory)
}

omics_enrichment_list <- function(omics1_list,omics2_list,outdir ="./",omics_name1="Proteomics",omics_name2="Phosphoproteomics",pvalueCutoff = 0.05,GO_showCategory=6,KEGG_showCategory=15){
  
  # predata <- enrichment_predata(pro_diff_path,phos_diff_path,phos_pro_path,outdir,log2FC,diff_p_adj)
  # omics1_list <- predata$protein_list                    # pro symbol list
  # omics2_list <- predata$phosphoprotein_list            # phos symbol list
  
  
  # 将 Gene Symbol 转换为 Entrez ID
  omics1_entrez_list <- bitr(omics1_list, 
                             fromType = "SYMBOL",      # 输入ID类型
                             toType = "ENTREZID",      # 输出ID类型
                             OrgDb = org.Hs.eg.db)$ENTREZID     # 物种数据库
  
  omics2_entrez_list <- bitr(omics2_list, 
                             fromType = "SYMBOL",      # 输入ID类型
                             toType = "ENTREZID",      # 输出ID类型
                             OrgDb = org.Hs.eg.db)$ENTREZID     # 物种数据库
#   # 单组学GO富集分析
#   if(length(omics1_list)>=1){
#   single_omics_GO_enrichment(omics1_list,omics_name1,outdir,pvalueCutoff,GO_showCategory)
#   }
#   if(length(omics2_list)>=1){
#   single_omics_GO_enrichment(omics2_list,omics_name2,outdir,pvalueCutoff,GO_showCategory)
#   }
#   # 单组学KEGG富集分析
#   if(length(omics1_entrez_list)>=2){
#     # print(omics1_entrez_list)
#   single_omics_KEGG_enrichment(omics1_entrez_list,omics_name1,outdir,pvalueCutoff,KEGG_showCategory)
#   }
#   if(length(omics2_entrez_list)>=2){
#   single_omics_KEGG_enrichment(omics2_entrez_list,omics_name2,outdir,pvalueCutoff,KEGG_showCategory)
# }
  # 两组学比较GO、KEGG富集分析
  if(length(omics1_list)>=2 &length(omics2_list)>=2){
  two_omicses_GO_enrichment(omics1_list,omics2_list,omics_name1,omics_name2,outdir,pvalueCutoff,GO_showCategory)
  }
  if(length(omics1_entrez_list)>=2 &length(omics2_entrez_list)>=2){
  two_omicses_KEGG_enrichment(omics1_entrez_list,omics2_entrez_list,omics_name1,omics_name2,outdir,pvalueCutoff,KEGG_showCategory)
  }
}


#####################run test ######################
#input: symbol id
# pro_diff_path <- "E:/pro_phosphpro/results/2.DE_Pro_Phospro/single_omics_differential_analysis/T_vs_N_Proteomics.tsv"
# phos_diff_path <- "E:/pro_phosphpro/results/2.DE_Pro_Phospro/single_omics_differential_analysis/T_vs_N_Phosphoproteomics.tsv"
# outdir <- "E:/pro_phosphpro/results/3.functional_integrative_analysis/fundamental_annoation_and_enrichment"
# phos_pro_path <- "E:/pro_phosphpro/caseData/phosphoproSite_Protein.tsv"
# log2FC <- 1.5
# diff_p_adj <- 0.05
# omics_name1 = "Proteomics"
# omics_name2 = "Phosphoproteomics"
# GO_showCategory = 6
# KEGG_showCategory = 15
# pvalueCutoff <- 0.05
#####################run test ######################


# opt <- parameter_validation(opt)
# pro_diff_path <- opt$pro_diff
# phos_diff_path  <- opt$phos_diff
# phos_pro_path <- opt$phos_pro
# outdir <- opt$outdir
# log2FC <- opt$log2FC
# diff_p_adj <- opt$diff_p_adj
# omics_name1 <- opt$omics_name1
# omics_name2 <- opt$omics_name2
# GO_showCategory <- opt$GO_showCategory
# KEGG_showCategory <- opt$KEGG_showCategory
# pvalueCutoff <- opt$pvalueCutoff
# pro_list <- c("BANF2")
# pro_list <- c("FGL1","CTRC")
# phos_list <- c("SYT6","HSPA5","SEC16B")
# # pro_list <- c("ALDH1A1","SULF1","MDK","INHBA","CTHRC1","THBS2","ADHFE1","THBS1","TMEM97","C1QTNF6")
# # phos_list <- c("PDHA1","PEPD","ALDH1A1","GSTA2","JPT1","KIF21A","TTYH1","LAYN","COBLL1","MYRIP","HMGA1","SETMAR","TACC1","PLPP3","HSPA5","TSNAX","XRCC1","SNTB2","RAB11FIP3","SEC16B","CDC42EP1","CDC42EP1")
# # 
# omics1_list <- pro_list
# omics2_list <- phos_list
# pvalueCutoff = 0.05
# outdir = "E:/pro_phosphpro/results/2.DE_Pro_Phospro/Stable_DifferentialNetwork/SubNetwork/subnet_9/functional_enrichment"
# omics_enrichment(pro_diff_path,phos_diff_path,phos_pro_path,outdir,log2FC,diff_p_adj,omics_name1,omics_name2,pvalueCutoff,GO_showCategory,KEGG_showCategory)

