# Differential_subnetwork_test
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))

Differential_subnetwork_plot <- function(subnet_dir,subnet_name,group_name = 'T-vs-N', omics1_name = 'Pro',
                                         omics2_name = 'Phos',
                                         edge_color_pos = "#9b6a65",
                                         edge_color_neg = "#5d8992", 
                                         Enhanced_in_N = "#5d8992", 
                                         Enhanced_in_T = "#9b6a65",
                                         Only_in_N = "#0c2b32",
                                         Only_in_T = "#381512",
                                         Conflict_relation = '#808080',
                                         fill_gradientn_color = c("#175663", "#dce6e9", "#90362d")){
nodes_path <- file.path(subnet_dir,paste0("nodes_",subnet_name,"_",group_name,"_",subnet_name,".txt"))
edges_path <- file.path(subnet_dir,paste0("edges_",subnet_name,"_",group_name,"_",subnet_name,".txt"))

group_parts <- strsplit(group_name, "-vs-")[[1]]
group1 <- group_parts[1]
group2 <- group_parts[2]

nodes <- read.table(nodes_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
edges <- read.table(edges_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 创建igraph对象
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# 准备节点属性
V(g)$case_norm_mean <- nodes$case_norm_mean
V(g)$omics_name <- nodes$omics_name

# 准备边属性 - 处理cor_case中的NA值
E(g)$cor_case <- edges$cor_case
E(g)$cor_T_type <- case_when(
  E(g)$cor_case < 0 ~ "< 0",
  E(g)$cor_case > 0 ~ "> 0",
  TRUE ~ "other"
)

# 直接删除 cor_status 为 "other" 的边
g_filtered <- delete_edges(g, E(g)[cor_T_type == "other"])
# g_filtered <- g


# 创建布局矩阵
layout_matrix <- as.matrix(nodes[, c("x", "y")])
############################################
############################################
############################################
# 详细检查形状映射问题
# cat("=== 形状映射问题诊断 ===\n")
# 
# # 检查节点数据中的Class列
# cat("nodes数据框中的Class列:\n")
# print(unique(nodes$Class))
# cat("Class列类型:", class(nodes$Class), "\n")
# 
# # 检查igraph对象中的Class属性
# cat("igraph对象中的Class属性:\n")
# print(unique(V(g_filtered)$Class))
# cat("Class属性类型:", class(V(g_filtered)$Class), "\n")
# 
# # 检查omics名称
# cat("omics1_name:", omics1_name, "\n")
# cat("omics2_name:", omics2_name, "\n")
# 
# # 检查形状映射
# cat("当前的形状映射:\n")
# current_mapping <- c(omics1_name = 21, omics2_name = 24)
current_mapping <- setNames(c(21, 24), c(omics1_name, omics2_name))
# print(current_mapping)
# 
# # 修复形状映射
# # 获取数据中实际存在的Class值
# actual_classes <- unique(V(g_filtered)$Class)
# cat("数据中实际的Class值:", actual_classes, "\n")
############################################
############################################
############################################

# 使用ggraph绘制
T_subnet <-ggraph(g_filtered, layout = layout_matrix) +
  # 绘制边
  geom_edge_link(
    aes(color = cor_T_type),
    width = 0.6,
    alpha = 0.7
  ) +
  # 绘制节点
  geom_node_point(
    aes(fill = case_norm_mean,shape = Class),
    size = 10,
    colour = "transparent",  # 使用透明颜色
    stroke = 0
  ) +
  # # 增强版本的节点标签
  # geom_node_text(
  #   aes(label = name),
  #   repel = TRUE,
  #   size = 3.5,
  #   fontface = "bold",
  #   color = "grey30",       # 中等灰色
  #   family = "sans",        # 字体家族
  #   max.overlaps = 30,      # 增加最大重叠容忍度
  #   box.padding = 2,        # 标签框的内边距
  #   point.padding = 0.8,    # 点周围的内边距
  #   force = 1,              # 排斥力大小
  #   min.segment.length = 0.1, # 最小连接线段长度
  #   segment.color = "grey70", # 连接线颜色
  #   segment.size = 0.3,     # 连接线粗细
  #   segment.alpha = 0.7     # 连接线透明度
  # ) +
  # 设置颜色标尺
  scale_edge_color_manual(
    name = "Correlation",
    # values = c("< 0" = "#487f8b", "> 0" = "#b3635b", "other" = "gray"),
    # values = c("< 0" = "#487f8b", "> 0" = "#b3635b", "other" = "transparent"),
    values = c("< 0" = edge_color_neg , "> 0" = edge_color_pos, "other" = "transparent"),
    labels = c("< 0", "> 0", "Other")
  ) +
  # 三色渐变
  scale_fill_gradientn(
    name = "Case Norm Mean",
    # colours = colorRampPalette(c("#175663", "#dce6e9", "#90362d"))(50)
    colours = colorRampPalette(fill_gradientn_color)(50)
  )+
  scale_shape_manual(
    name = "Omics Type",
    values =current_mapping,  # 使用21和24（带填充的形状）
    guide = guide_legend(override.aes = list(fill = "grey50", size = 6))  # 添加填充色
  ) +
  # 主题设置
  theme_void() +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 60, face = "bold"),
    # 控制图例字体大小
    legend.title = element_text(size = 36, face = "bold"),      # 图例标题字体大小
    legend.text = element_text(size = 36),                     # 图例文字字体大小
  ) +
  # 添加标题
  ggtitle(paste0("Group ",group1," net")) +
  # 图例布局
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
T_subnet_lable <- T_subnet +
  # 增强版本的节点标签
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 14,
    fontface = "bold",
    color = "grey30",       # 中等灰色
    family = "sans",        # 字体家族
    max.overlaps = 30,      # 增加最大重叠容忍度
    box.padding = 2,        # 标签框的内边距
    point.padding = 0.8,    # 点周围的内边距
    force = 10,              # 排斥力大小
    min.segment.length = 0.1, # 最小连接线段长度
    segment.color = "grey70", # 连接线颜色
    segment.size = 0.3,     # 连接线粗细
    segment.alpha = 0.7     # 连接线透明度
  ) 
T_subnet_lable
ggsave(file.path(subnet_dir,paste0(subnet_name,"_",strsplit(group_name,"-vs-")[[1]][1],"_plot.png")),plot = T_subnet_lable, width = 10, height = 8, dpi = 300)

# 准备边属性 - 处理cor_case中的NA值
E(g)$cor_control <- edges$cor_control
E(g)$cor_N_type <- case_when(
  E(g)$cor_control < 0 ~ "< 0",
  E(g)$cor_control > 0 ~ "> 0",
  TRUE ~ "other"
)
# 直接删除 cor_status 为 "other" 的边
g_filtered <- delete_edges(g, E(g)[cor_N_type == "other"])
# g_filtered <- g
N_subnet <-ggraph(g_filtered, layout = layout_matrix) +
  # 绘制边
  geom_edge_link(
    aes(color = cor_N_type),
    width = 0.6,
    alpha = 0.7
  ) +
  # 绘制节点
  geom_node_point(
    aes(fill = control_norm_mean,shape = Class),
    size = 10,
    colour = "transparent",  # 使用透明颜色
    stroke = 0
  ) +
  # # 增强版本的节点标签
  # geom_node_text(
  #   aes(label = name),
  #   repel = TRUE,
  #   size = 3.5,
  #   fontface = "bold",
  #   color = "grey30",       # 中等灰色
  #   family = "sans",        # 字体家族
  #   max.overlaps = 30,      # 增加最大重叠容忍度
  #   box.padding = 2,        # 标签框的内边距
  #   point.padding = 0.8,    # 点周围的内边距
  #   force = 1,              # 排斥力大小
  #   min.segment.length = 0.1, # 最小连接线段长度
  #   segment.color = "grey70", # 连接线颜色
  #   segment.size = 0.3,     # 连接线粗细
  #   segment.alpha = 0.7     # 连接线透明
  # ) +
  # 设置颜色标尺
  scale_edge_color_manual(
    name = "Correlation",
    # values = c("< 0" = "#487f8b", "> 0" = "#b3635b", "other" = "gray"),
    # values = c("< 0" = "#487f8b", "> 0" = "#b3635b", "other" = "transparent"),
    values = c("< 0" = edge_color_neg , "> 0" = edge_color_pos, "other" = "transparent"),
    labels = c("< 0", "> 0", "Other")
  ) +
  # 三色渐变
  scale_fill_gradientn(
    name = "Case Norm Mean",
    # colours = colorRampPalette(c("#175663", "#dce6e9", "#90362d"))(50)
    colours = colorRampPalette(fill_gradientn_color)(50)
    
  )+
  scale_shape_manual(
    name = "Omics Type",
    values = current_mapping,  # 使用21和24（带填充的形状）
    guide = guide_legend(override.aes = list(fill = "grey50", size = 4))  # 添加填充色
  ) +
  # 主题设置
  theme_void() +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 60, face = "bold"),
    # 控制图例字体大小
    legend.title = element_text(size = 36, face = "bold"),      # 图例标题字体大小
    legend.text = element_text(size = 36),                     # 图例文字字体大小
  ) +
  # 添加标题
  ggtitle(paste0("Group ",group2," net")) +
  # 图例布局
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
N_subnet_lable <- N_subnet+
  # 增强版本的节点标签
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 14,
    fontface = "bold",
    color = "grey30",       # 中等灰色
    family = "sans",        # 字体家族
    max.overlaps = 30,      # 增加最大重叠容忍度
    box.padding = 2,        # 标签框的内边距
    point.padding = 0.8,    # 点周围的内边距
    force = 10,              # 排斥力大小
    min.segment.length = 0.1, # 最小连接线段长度
    segment.color = "grey70", # 连接线颜色
    segment.size = 0.3,     # 连接线粗细
    segment.alpha = 0.7     # 连接线透明度
  ) 
ggsave(file.path(subnet_dir,paste0(subnet_name,"_",strsplit(group_name,"-vs-")[[1]][2],"_plot.png")),plot = N_subnet_lable, width = 10, height = 8, dpi = 300)



g_filtered <- delete_edges(g, E(g)[cor_status == "Non-significant"])


# 先创建命名向量
edge_color_vector <- c(Enhanced_in_N, Enhanced_in_T, Only_in_N, Only_in_T, Conflict_relation)

names(edge_color_vector) <- c(
  paste0("Enhanced in ", group2),
  paste0("Enhanced in ", group1), 
  paste0("Only in ", group2),
  paste0("Only in ", group1),
  "Conflict relation"
)

edge_linetype_vector <- c( "dashed","dashed", "dotted", "solid","solid","dotted")
names(edge_linetype_vector) <- c(
  paste0("Enhanced in ", group2),
  paste0("Enhanced in ", group1), 
  "Non-significant",
  paste0("Only in ", group2),
  paste0("Only in ", group1),
  "Conflict relation"
)


# 使用ggraph绘制
diff_subnet <-ggraph(g_filtered, layout = layout_matrix) +
  # 绘制边
  geom_edge_link(
    aes(color = cor_status,
        linetype = cor_status),
    width = 0.6,
    alpha = 0.5
  ) +
  # 绘制节点
  geom_node_point(
    aes(fill = Log2FC,shape = Class),
    size = 10,
    colour = "transparent",  # 使用透明颜色
    stroke = 0
  ) +
  # 增强版本的节点标签
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 16,
    fontface = "bold",
    color = "grey30",       # 中等灰色
    family = "sans",        # 字体家族
    max.overlaps = 30,      # 增加最大重叠容忍度
    box.padding = 2,        # 标签框的内边距
    point.padding = 0.8,    # 点周围的内边距
    force = 10,              # 排斥力大小
    min.segment.length = 0.1, # 最小连接线段长度
    segment.color = "grey70", # 连接线颜色
    segment.size = 0.3,     # 连接线粗细
    segment.alpha = 0.7     # 连接线透明度
  ) +
  # 设置颜色标尺
  scale_edge_color_manual(
    name = "Correlation",
    # values = c("Enhanced in N" = "#5d8992", "Enhanced in T" = "#9b6a65","Only in N" = "#0c2b32","Only in T" = "#381512"),
    # values = c("Enhanced in N" = Enhanced_in_N, "Enhanced in T" = Enhanced_in_T ,"Only in N" = Only_in_N ,"Only in T" = Only_in_T ),
    # values = c(paste0("Enhanced in ",strsplit(group_name,"-vs-")[[1]][2]) = Enhanced_in_N, paste0("Enhanced in ",strsplit(group_name,"-vs-")[[1]][1]) = Enhanced_in_T ,
    #            paste0("Enhanced in ",strsplit(group_name,"-vs-")[[1]][2])"Only in N" = Only_in_N ,"Only in T" = Only_in_T ),
    values = edge_color_vector
  ) +
  # 添加线型标尺（必须的）
  scale_edge_linetype_manual(
    name = "Correlation",
    # values = c("Enhanced in N" = "dashed", "Enhanced in T" = "dashed", "Non-significant" = "dotted","Only in N" = "solid","Only in T" ="solid")
    values = edge_linetype_vector
  ) +
  scale_shape_manual(
    name = "Omics Type",
    values = current_mapping,  # 使用21和24（带填充的形状）
    guide = guide_legend(override.aes = list(fill = "grey50", size = 4))  # 添加填充色
  ) +
  # 三色渐变
  scale_fill_gradientn(
    # name = "Log2FC (T:N)",
    name = paste0("Log2FC (",group1,":",group2,")"),
    # colours = colorRampPalette(c("#175663", "#dce6e9", "#90362d"))(50)
    colours = colorRampPalette(fill_gradientn_color)(50)
  )+
  # 主题设置
  theme_void() +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 60, face = "bold"),
    # 控制图例字体大小
    legend.title = element_text(size = 36, face = "bold"),      # 图例标题字体大小
    legend.text = element_text(size = 36),                     # 图例文字字体大小
  ) +
  # 添加标题
  ggtitle(paste0(group_name," net")) +
  # 图例布局
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
ggsave(file.path(subnet_dir,paste0(subnet_name,"_diff_plot.png")),plot = diff_subnet, width = 10, height = 8, dpi = 300)


# 将T_subnet 和 N_subnet 放在一张图上左右
combined_plot <- T_subnet + N_subnet + diff_subnet+ plot_layout(ncol = 3)
ggsave(file.path(subnet_dir,paste0(subnet_name,"_",group_name,"_plot.png")), plot = combined_plot, width = 30, height = 8, dpi = 300)

}


# subnet_name = "subnet_1"
# group_name = 'T-vs-N'
# subnet_dir <- "E:/pro_phos_test/results/2.Differential_Analysis/2.3Stable_Differential_Network/2.3.2Network_Clustering/2.3.3Sub_Network/subnet_1"
# Differential_subnetwork_plot(subnet_dir,subnet_name)

# ###all test
# 
# Differential_overall_network_plot <- function(outdir,nodesize = 10,omics1_name = 'Pro',
#                                               omics2_name = 'Phos', 
#                                               edge_color_pos = "#9b6a65",
#                                               edge_color_neg = "#5d8992", 
#                                               Enhanced_in_N = "#5d8992", 
#                                               Enhanced_in_T = "#9b6a65",
#                                               Only_in_N = "#0c2b32",
#                                               Only_in_T = "#381512",
#                                               fill_gradientn_color = c("#175663", "#dce6e9", "#90362d")){
#   # setwd(outdir)
#   nodes_path <- file.path(outdir,"nodes_T-vs-N.txt")
#   edges_path <- file.path(outdir,"edges_T-vs-N.txt")
#   
#   nodes <- read.table(nodes_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   edges <- read.table(edges_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   edges <- edges[!(edges$cor_status %in% c("Conflict relation","Non-significant")),]
#   
#   # 创建igraph对象
#   g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
#   
#   # 准备节点属性
#   V(g)$case_norm_mean <- nodes$case_norm_mean
#   V(g)$omics_name <- nodes$omics_name
#   
#   # 准备边属性 - 处理cor_case中的NA值
#   E(g)$cor_case <- edges$cor_case
#   E(g)$cor_T_type <- case_when(
#     E(g)$cor_case < 0 ~ "< 0",
#     E(g)$cor_case > 0 ~ "> 0",
#     TRUE ~ "other"
#   )
#   
#   # 直接删除 cor_status 为 "other" 的边
#   g_filtered <- delete_edges(g, E(g)[cor_T_type == "other"])
#   # g_filtered <- g
#   
#   
#   # 创建布局矩阵
#   layout_matrix <- as.matrix(nodes[, c("x", "y")])
#   
#   # 使用ggraph绘制
#   T_subnet <-ggraph(g_filtered, layout = layout_matrix) +
#     # 绘制边
#     geom_edge_link(
#       aes(color = cor_T_type),
#       width = 0.6,
#       alpha = 0.7
#     ) +
#     # 绘制节点
#     geom_node_point(
#       aes(fill = case_norm_mean,shape = Class),
#       size = nodesize,
#       colour = "transparent",  # 使用透明颜色
#       stroke = 0
#     ) +
#     # # 增强版本的节点标签
#     # geom_node_text(
#     #   aes(label = name),
#     #   repel = TRUE,
#     #   size = 3.5,
#     #   fontface = "bold",
#     #   color = "grey30",       # 中等灰色
#     #   family = "sans",        # 字体家族
#     #   max.overlaps = 30,      # 增加最大重叠容忍度
#     #   box.padding = 2,        # 标签框的内边距
#     #   point.padding = 0.8,    # 点周围的内边距
#     #   force = 1,              # 排斥力大小
#     #   min.segment.length = 0.1, # 最小连接线段长度
#     #   segment.color = "grey70", # 连接线颜色
#     #   segment.size = 0.3,     # 连接线粗细
#     #   segment.alpha = 0.7     # 连接线透明度
#     # ) +
#     # 设置颜色标尺
#     scale_edge_color_manual(
#       name = "Correlation",
#       # values = c("< 0" = "#487f8b", "> 0" = "#b3635b", "other" = "gray"),
#       values = c("< 0" = edge_color_neg , "> 0" = edge_color_pos, "other" = "transparent"),
#       labels = c("< 0", "> 0", "Other")
#     ) +
#   
#     # 三色渐变
#     scale_fill_gradientn(
#       name = "Case Norm Mean",
#       # colours = colorRampPalette(c("#175663", "#dce6e9", "#90362d"))(50)
#       colours = colorRampPalette(fill_gradientn_color)(50)
#     )+
#     scale_shape_manual(
#       name = "Omics Type",
#       values = current_mapping,  # 使用21和24（带填充的形状）
#       guide = guide_legend(override.aes = list(fill = "grey50", size = 6))  # 添加填充色
#     ) +
#     # 主题设置
#     theme_void() +
#     theme(
#       legend.position = "right",
#       plot.background = element_rect(fill = "white", color = NA),
#       plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
#       # 控制图例字体大小
#       legend.title = element_text(size = 16, face = "bold"),      # 图例标题字体大小
#       legend.text = element_text(size = 16),                     # 图例文字字体大小
#     ) +
#     # 添加标题
#     ggtitle("Group T overallnet") +
#     # 图例布局
#     guides(
#       fill = guide_colorbar(barheight = unit(4, "cm")),
#       shape = guide_legend(
#         override.aes = list(
#           fill = "grey80",
#           # color = "black",    # 关键：覆盖透明边框
#           size = 6,
#           stroke = 0.5
#         ))
#     )
#   T_subnet_lable <- T_subnet +
#     # 增强版本的节点标签
#     geom_node_text(
#       aes(label = name),
#       repel = TRUE,
#       size = 5,
#       fontface = "bold",
#       color = "grey30",       # 中等灰色
#       family = "sans",        # 字体家族
#       max.overlaps = 10,      # 增加最大重叠容忍度
#       box.padding = 2,        # 标签框的内边距
#       point.padding = 0.8,    # 点周围的内边距
#       force = 10,              # 排斥力大小
#       min.segment.length = 0.1, # 最小连接线段长度
#       segment.color = "grey70", # 连接线颜色
#       segment.size = 0.3,     # 连接线粗细
#       segment.alpha = 0.7     # 连接线透明度
#     ) 
#   ggsave(file.path(outdir,"T_plot.png"),plot = T_subnet_lable, width = 10, height = 8, dpi = 300)
#   
#   # 准备边属性 - 处理cor_case中的NA值
#   E(g)$cor_control <- edges$cor_control
#   E(g)$cor_N_type <- case_when(
#     E(g)$cor_control < 0 ~ "< 0",
#     E(g)$cor_control > 0 ~ "> 0",
#     TRUE ~ "other"
#   )
#   # 直接删除 cor_status 为 "other" 的边
#   g_filtered <- delete_edges(g, E(g)[cor_N_type == "other"])
#   # g_filtered <- g
#   N_subnet <-ggraph(g_filtered, layout = layout_matrix) +
#     # 绘制边
#     geom_edge_link(
#       aes(color = cor_N_type),
#       width = 0.6,
#       alpha = 0.7
#     ) +
#     # 绘制节点
#     geom_node_point(
#       aes(fill = control_norm_mean,shape = Class),
#       size = nodesize,
#       colour = "transparent",  # 使用透明颜色
#       stroke = 0
#     ) +
#     # # 增强版本的节点标签
#     # geom_node_text(
#     #   aes(label = name),
#     #   repel = TRUE,
#     #   size = 3.5,
#     #   fontface = "bold",
#     #   color = "grey30",       # 中等灰色
#     #   family = "sans",        # 字体家族
#     #   max.overlaps = 30,      # 增加最大重叠容忍度
#     #   box.padding = 2,        # 标签框的内边距
#     #   point.padding = 0.8,    # 点周围的内边距
#     #   force = 1,              # 排斥力大小
#     #   min.segment.length = 0.1, # 最小连接线段长度
#     #   segment.color = "grey70", # 连接线颜色
#     #   segment.size = 0.3,     # 连接线粗细
#     #   segment.alpha = 0.7     # 连接线透明
#     # ) +
#     # 设置颜色标尺
#     scale_edge_color_manual(
#       name = "Correlation",
#       # values = c("< 0" = "#487f8b", "> 0" = "#b3635b", "other" = "gray"),
#       values = c("< 0" = edge_color_neg , "> 0" = edge_color_pos, "other" = "transparent"),
#       labels = c("< 0", "> 0", "Other")
#     ) +
#     # 三色渐变
#     scale_fill_gradientn(
#       name = "Case Norm Mean",
#       # colours = colorRampPalette(c("#175663", "#dce6e9", "#90362d"))(50)
#       colours = colorRampPalette(fill_gradientn_color)(50)
#     )+
#     scale_shape_manual(
#       name = "Omics Type",
#       values = current_mapping,  # 使用21和24（带填充的形状）
#       guide = guide_legend(override.aes = list(fill = "grey50", size = 4))  # 添加填充色
#     ) +
#     # 主题设置
#     theme_void() +
#     theme(
#       legend.position = "right",
#       plot.background = element_rect(fill = "white", color = NA),
#       plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
#       # 控制图例字体大小
#       legend.title = element_text(size = 16, face = "bold"),      # 图例标题字体大小
#       legend.text = element_text(size = 16),                     # 图例文字字体大小
#     ) +
#     # 添加标题
#     ggtitle("Group N overallnet") +
#     # 图例布局
#     guides(
#       fill = guide_colorbar(barheight = unit(4, "cm")),
#       shape = guide_legend(
#         override.aes = list(
#           fill = "grey80",
#           # color = "black",    # 关键：覆盖透明边框
#           size = 6,
#           stroke = 0.5
#         ))
#     )
#   N_subnet_lable <- N_subnet+
#     # 增强版本的节点标签
#     geom_node_text(
#       aes(label = name),
#       repel = TRUE,
#       size = 5,
#       fontface = "bold",
#       color = "grey30",       # 中等灰色
#       family = "sans",        # 字体家族
#       max.overlaps = 10,      # 增加最大重叠容忍度
#       box.padding = 2,        # 标签框的内边距
#       point.padding = 0.8,    # 点周围的内边距
#       force = 10,              # 排斥力大小
#       min.segment.length = 0.1, # 最小连接线段长度
#       segment.color = "grey70", # 连接线颜色
#       segment.size = 0.3,     # 连接线粗细
#       segment.alpha = 0.7     # 连接线透明度
#     ) 
#   ggsave(file.path(outdir,"N_plot.png"),plot = N_subnet_lable, width = 10, height = 8, dpi = 300)
#   
#   
#   
#   g_filtered <- delete_edges(g, E(g)[cor_status == "Non-significant"])
#   # g_filtered <- g
#   
#   # 使用ggraph绘制
#   diff_subnet <-ggraph(g_filtered, layout = layout_matrix) +
#     # 绘制边
#     geom_edge_link(
#       aes(color = cor_status,
#           linetype = cor_status),
#       width = 0.6,
#       alpha = 0.5
#     ) +
#     # 绘制节点
#     geom_node_point(
#       aes(fill = Log2FC,shape = Class),
#       size = nodesize,
#       colour = "transparent",  # 使用透明颜色
#       stroke = 0
#     ) +
#     # 增强版本的节点标签
#     geom_node_text(
#       aes(label = name),
#       repel = TRUE,
#       size = 5,
#       fontface = "bold",
#       color = "grey30",       # 中等灰色
#       family = "sans",        # 字体家族
#       max.overlaps = 30,      # 增加最大重叠容忍度
#       box.padding = 2,        # 标签框的内边距
#       point.padding = 0.8,    # 点周围的内边距
#       force = 10,              # 排斥力大小
#       min.segment.length = 0.1, # 最小连接线段长度
#       segment.color = "grey70", # 连接线颜色
#       segment.size = 0.3,     # 连接线粗细
#       segment.alpha = 0.7     # 连接线透明度
#     ) +
#     # 设置颜色标尺
#     scale_edge_color_manual(
#       name = "Correlation",
#       # values = c("Enhanced in N" = "#5d8992", "Enhanced in T" = "#9b6a65", "Non-significant" = "gray60","Only in N" = "#0c2b32","Only in T" = "#381512"),
#       values = c("Enhanced in N" = Enhanced_in_N, "Enhanced in T" = Enhanced_in_T ,"Only in N" = Only_in_N ,"Only in T" = Only_in_T ),
#     ) +
#     # 添加线型标尺（必须的）
#     scale_edge_linetype_manual(
#       name = "Correlation",
#       values = c("Enhanced in N" = "dashed", "Enhanced in T" = "dashed", "Non-significant" = "dotted","Only in N" = "solid","Only in T" ="solid")
#     ) +
#     scale_shape_manual(
#       name = "Omics Type",
#       values = current_mapping,  # 使用21和24（带填充的形状）
#       guide = guide_legend(override.aes = list(fill = "grey50", size = 4))  # 添加填充色
#     ) +
#     # 三色渐变
#     scale_fill_gradientn(
#       name = "Log2FC (T:N)",
#       # colours = colorRampPalette(c("#175663", "#dce6e9", "#90362d"))(50)
#       colours = colorRampPalette(fill_gradientn_color)(50)
#     )+
#     # 主题设置
#     theme_void() +
#     theme(
#       legend.position = "right",
#       plot.background = element_rect(fill = "white", color = NA),
#       plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
#       # 控制图例字体大小
#       legend.title = element_text(size = 16, face = "bold"),      # 图例标题字体大小
#       legend.text = element_text(size = 16),                     # 图例文字字体大小
#     ) +
#     # 添加标题
#     ggtitle("T VS N overallnet") +
#     # 图例布局
#     guides(
#       fill = guide_colorbar(barheight = unit(4, "cm")),
#       shape = guide_legend(
#         override.aes = list(
#           fill = "grey80",
#           # color = "black",    # 关键：覆盖透明边框
#           size = 6,
#           stroke = 0.5
#         ))
#     )
#   ggsave(file.path(outdir,"diff_plot.png"),plot = diff_subnet, width = 10, height = 8, dpi = 300)
#   
#   
#   # 将T_subnet 和 N_subnet 放在一张图上左右
#   combined_plot <- T_subnet + N_subnet + diff_subnet+ plot_layout(ncol = 3)
#   ggsave(file.path(outdir,"T_vs_N_plot.png"), plot = combined_plot, width = 30, height = 8, dpi = 300)
#   
# }




# outdir <- "E:/pro_phosphpro/results/2.DE_Pro_Phospro/Stable_DifferentialNetwork/OverallNetwork"
# Differential_overall_network_plot(outdir,nodesize = 10)

# subnet_dir <- "D:/datasest2/results/2.Differential_Analysis/2.3Stable_Differential_Network/2.3.3Sub_Network/subnet_1"
# group_name <- "Tumor-vs-Normal"
# subnet_name <- 'subnet_1'
# omics1_name <- "Protein"
# omics2_name <- "PhosProtein"
# 
# Differential_subnetwork_plot(subnet_dir,subnet_name,group_name = group_name, omics1_name = omics1_name,
#                                   omics2_name = omics2_name)
