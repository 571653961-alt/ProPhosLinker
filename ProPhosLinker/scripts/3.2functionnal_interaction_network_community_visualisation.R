# community_detection_plot
# Description: Detects community-related activities in network traffic.

options(warn = -1)
suppressWarnings(library(optparse))   # 加载 optparse 包
option_list <- list(
  #input parameter
  make_option(c("-e", "--edge_file"), type = "character", default = NULL, 
              metavar = "character", help = "Input the edge file path"),
  make_option(c("-c", "--community_info"), type = "character", default = NULL, 
              metavar = "character", help = "Input the community info file path"),
  #output parameter
  make_option(c("-o", "--outdir"), type = "character", default = NULL, metavar = "character",
              help = "Output directory [Default: %default]"),
  make_option(c("-t", "--top_module_num"), type = "integer", default = 6, metavar = "INT",
              help = "top module number [default: %default]"),
  #general parameter
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE, metavar = "FLAG",
              help = "Print detailed output messages [Default: %default]")
)
# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                           description = "the Community Detection Plot of Protein-Phosphorylation network")
opt <- parse_args(opt_parser)

suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(tidygraph)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggforce)))
suppressWarnings(suppressPackageStartupMessages(library(ggpubr)))
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(ggrepel)))



parameter_validation <- function(opt){
  # 手动检查必需参数
  if (is.null(opt$edge_file)) {
    stop("Error: edge file is required. Use -e or --edge_file to specify the edge file.")
  }
  if (is.null(opt$community_info)) {
    stop("Error: community info file is required. Use -c or --community_info to specify the community info file.")
  }
  #path cleaning
  clean_path <- function(path) {
    # Trim leading and trailing whitespace and quotes
    path <- gsub("^['\" ]+|['\" ]+$","",path)
    # Replace backslashes (\) with forward slashes (/)
    path <- gsub("\\\\","/",path)
    return(path)
  }
  check_input_file <- function(file_path,file_name){
    if(!file.exists(file_path)){
      stop("❌ Error: The file '", file_name, "' does not exist at the specified path: ", file_path)
    }
    if(file.size(file_path) == 0){
      stop("❌ Error: The file '", file_name, "' is empty at the specified path: ", file_path)
    }
  }
  check_output_dir <- function(result_dir){
    if(is.null(result_dir)|| (length(result_dir) == 0)){
      result_dir <- getwd()
    }
    if (!dir.exists(result_dir)) {
      # Recursively create directory and all parent directories
      dir.create(result_dir, recursive = TRUE)
      print(paste("✅ Output directory created:", result_dir))
    } else {
      print(paste("✅ Output directory already exists:", result_dir))
    }
    return(result_dir)
  }
  
  opt$edge_file <- clean_path(opt$edge_file)
  opt$community_info <- clean_path(opt$community_info)
  opt$outdir <- clean_path(opt$outdir)
  opt$outdir  <- check_output_dir(opt$outdir)
  check_input_file(opt$edge_file,"the edge file")
  
  if (opt$verbose) { print(opt)}
  return(opt)
}

opt <- parameter_validation(opt)
edge_file_path <- opt$edge_file
community_info <- opt$community_info
outdir <- opt$outdir
# top_module_num <- 6
top_module_num <- opt$top_module_num
setwd(outdir)

# ##### test #######
# # node_file_path <- "E:\\pro_phosphpro\\results\\Network_analysis\\Node_diff_info.tsv"
# edge_file_path <- "E:/pro_phosphpro/results/3.functional_integrative_analysis/Network_analysis/Edge_diff_info.tsv"
# community_info <- "E:/pro_phosphpro/results/3.functional_integrative_analysis/Network_analysis/DEP_community_detection\\infomap_community_detection_modules.tsv"
# setwd("E:/pro_phosphpro/results/3.functional_integrative_analysis/Network_analysis/DEP_community_detection")
# ##### test #######

# 0.Load node and edge data
# node_df <- read.csv(node_file_path, sep = "\t", header = TRUE)
edge_df <- read.csv(edge_file_path, sep = "\t", header = TRUE)
community_info <- read.csv(community_info, sep = "\t", header = TRUE)

# length(unique(node_df$NodeName))

# 1. Create igraph object, and assign module IDs and Numbers to nodes [A whole graph]
# Create igraph object from data frames
g <- graph_from_data_frame(d = edge_df, directed = TRUE)
# 创建节点到模块ID 和 number 的映射表
module_mapping <- community_info %>%
  separate_rows(ModuleItems, sep = ';') %>%
  rename(node_name = ModuleItems) %>%
  mutate(node_name = trimws(node_name))%>%
  #ModuleSize降序排列
  arrange(desc(ModuleSize))

# head(module_mapping)
# Add module IDs and Numbers to nodes in the igraph object
all_nodes <- V(g)$name
module_ids <- rep(NA,length(all_nodes))
module_nums <- rep(NA,length(all_nodes))
module_labels <- rep(NA,length(all_nodes))
for (i in seq_along(all_nodes)) {
  node <- all_nodes[i]
  match_index <- which(module_mapping$node_name == node)
  if (length(match_index) > 0) {
    module_ids[i] <- as.character(module_mapping$ModuleID[match_index[1]])
    module_nums[i] <- module_mapping$ModuleSize[match_index[1]]
    module_labels[i] <- paste0(module_ids[i], " (", module_nums[i],")")
  }
}
V(g)$module_id <- module_ids
V(g)$module_num <- module_nums
V(g)$module_label <- module_labels

# Visualize the entire graph with modules
g_tbl <- as_tbl_graph(g)

ModuleSize_show = 0
unique_modules <- unique(vertex_attr(g_tbl)$module_id[vertex_attr(g_tbl)$module_num >= ModuleSize_show & vertex_attr(g_tbl)$module_id != -1])


# 使用多种颜色方案的组合
if (length(unique_modules) > 0) {
  # 组合多个专业调色板
  colors_from_palettes <- c(
    RColorBrewer::brewer.pal(9, "Set1"),
    RColorBrewer::brewer.pal(8, "Dark2"),
    
    RColorBrewer::brewer.pal(8, "Set2"),
    RColorBrewer::brewer.pal(12, "Set3"),
    viridis::viridis(20)  # 补充更多颜色
  )
  # 取需要的数量
  module_colors <- colors_from_palettes[1:length(unique_modules)]
  color_mapping <- setNames(module_colors, unique_modules)
} else {
  color_mapping <- c()
}

# 为图中的所有节点设置颜色
g_tbl <- g_tbl %>%
  activate(nodes) %>%
  mutate(node_color = ifelse(module_id %in% unique_modules & !is.na(module_id),
                       as.character(module_id),
                       paste0("other (num <", ModuleSize_show, ")")),
         node_color = factor(node_color,
                             levels = c(unique_modules,paste0("other (num <", ModuleSize_show, ")"))),
         module_id = as.factor(module_id)
         )
# 
if(ModuleSize_show >1){


full_color_mapping <- c(color_mapping, "gray70")
names(full_color_mapping)[length(full_color_mapping)] <- paste0("other (num <", ModuleSize_show, ")")
}else{
  full_color_mapping <- color_mapping
}
# 首先创建一个module_id到module_label的映射
module_label_mapping <- setNames(
  vertex_attr(g_tbl)$module_label,
  vertex_attr(g_tbl)$module_id
)
# 去重，确保每个module_id只有一个对应的module_label
module_label_mapping <- module_label_mapping[!duplicated(names(module_label_mapping))]



# 获取模块信息并排序
module_info <- data.frame(
  module_id = vertex_attr(g_tbl)$module_id,
  module_num = vertex_attr(g_tbl)$module_num,
  module_label = vertex_attr(g_tbl)$module_label,
  node_color = vertex_attr(g_tbl)$node_color
) %>% 
  distinct() %>%
  arrange(ifelse(module_id == -1, Inf, -module_num), module_id)

# 重新排序模块因子水平
vertex_attr(g_tbl)$module_label <- factor(
  vertex_attr(g_tbl)$module_label,
  levels = module_info$module_label
)
vertex_attr(g_tbl)$module_id <- factor(
  vertex_attr(g_tbl)$module_id,
  levels = module_info$module_id
)
vertex_attr(g_tbl)$node_color <- factor(
  vertex_attr(g_tbl)$node_color,
  levels = unique(module_info$node_color)
)

########  all modules color  begin ######
# p <- ggraph(g_tbl, layout = "fr") +
#   geom_edge_link(
#     color = "grey50",
#     alpha = 0.4,
#     width = 0.4
#   ) +
#   # 为每个模块添加背景区域（去边框）
#   ggforce::geom_mark_hull(
#     aes(x, y,
#         group = module_id,
#         fill = module_id
#     ),
#     alpha = 0.1,
#     concavity = 6,
#     expand = unit(2, "mm"),
#     radius = unit(3, "mm"),
#     lwd = 0.05,
#     color = 'grey90'
#   ) +
#   geom_node_point(
#     aes(color = node_color, size = centrality_degree()),
#     alpha = 0.9
#   ) +
#   # 设置颜色 scale - 修改 labels 参数
#   scale_color_manual(
#     name = "Module ID (Item Num.)",
#     values = full_color_mapping,
#     drop = FALSE,
#     labels = function(x) {
#       sapply(x, function(id) {
#         if (id %in% names(module_label_mapping)) {
#           module_label_mapping[as.character(id)]
#         } else {
#           id
#         }
#       })
#     }
#   ) +
#   scale_fill_manual(values = full_color_mapping, na.value = "gray70",
#                     guide = "none") +
#   scale_size_continuous(name = "Nodes Centrality",
#                         range = c(2, 6)) +
#   theme_void() +
#   theme(strip.text = element_text(face = "bold"),
#         plot.background = element_rect(fill = "white", color = NA),
#         plot.title = element_text(
#           hjust = 0.5,           # 水平居中
#           face = "bold",         # 加粗
#           size = 18,             # 加大字体
#           margin = margin(t = 15, b = 15)  # 上下边距各15个单位
#         )) +
#   labs(title = "Differential Protein and Phosphorylation Sites Network with Expanded Modules")
# p
# ggsave("network_expanded_modules.png", plot = p, width = 12, height =10, dpi = 300)



# unique_modules_lables <- unique(vertex_attr(g_tbl)$module_label[vertex_attr(g_tbl)$module_num > ModuleSize_show & vertex_attr(g_tbl)$module_id != -1])
# color_mapping_labels <- setNames(module_colors, unique_modules_lables)

# # 选项3：使用facet按模块分面显示
# p3 <- ggraph(g_tbl, layout = "fr") +
#   geom_edge_link(color = "grey30", alpha = 0.3) +
#   geom_node_point(aes(color = module_id,size = centrality_degree())) +
#   facet_nodes(~ module_label, scales = "free") +
#   scale_color_manual(name = "Module ID (Item Num.)",
#                      values = full_color_mapping,
#                      labels = function(x) {
#                        sapply(x, function(id) {
#                          if (id %in% names(module_label_mapping)) {
#                            module_label_mapping[as.character(id)]
#                          } else {
#                            id
#                          }
#                        })
#                      },
#                      na.value = "gray80") +
#   scale_size_continuous(name = "Nodes Centrality",  # 修改图例名称
#                         range = c(2, 6)) +
#   theme_void() +
#   theme(strip.text = element_text(
#     face = "bold",
#     size = 10,  # 调节小图标题大小，可以根据需要调整
#     margin = margin(b = 5)  # 可选：调节标题与图之间的间距
#   ),
#         plot.background = element_rect(fill = "white", color = NA),
#         plot.title = element_text(
#           hjust = 0.5,           # 水平居中
#           face = "bold",         # 加粗
#           size = 30,             # 加大字体
#           margin = margin(t = 15, b = 15)  # 上下边距各15个单位
#         )) +
#   labs(title = "Differential Protein and Phosphorylation Sites Network by Module")
# p3
# # 保存图形
# ggsave("network_by_module.png", plot = p3, width = 18, height =16, dpi = 300)
# ggsave("network_by_module.pdf", plot = p3, width = 18, height =16)
#  all modules color  end #

###################################
# top_module_num <- 6
top_modules <- unique(module_mapping$ModuleID)
#删除top_modules == -1
top_modules <- top_modules[top_modules != -1][1:top_module_num]
top_module_mapping <- module_mapping[module_mapping$ModuleID %in% top_modules, ]
top_module_nodes <- unique(top_module_mapping$node_name)
# 提取g中的top_module_nodes节点作为子图
top_module_g <- induced_subgraph(g_tbl, vids = V(g)$name %in% top_module_nodes)
top_module_g_tbl <- as_tbl_graph(top_module_g)

top_module_full_color_mapping <- full_color_mapping[names(full_color_mapping) %in% top_modules]


# # 1. 在图中计算并添加 centrality_degree 作为节点列
# top_module_g_with_cent <- top_module_g_tbl %>%
#   activate(nodes) %>%
#   mutate(
#     deg_cent = centrality_degree(),
#     is_top6 = rank(-deg_cent, ties.method = "min") <= 6  # 标记前6
#   )
# 
# # 2. 提取前6个节点名称（用于 geom_node_text 的 filter）
# top6_names <- top_module_g_with_cent %>%
#   activate(nodes) %>%
#   as_tibble() %>%
#   arrange(desc(deg_cent)) %>%
#   head(16) %>%
#   pull(name)


g <- top_module_g_tbl

# Step 1: 将图转为 igraph（如果还不是）
if (!inherits(g, "igraph")) {
  ig <- as.igraph(g)
} else {
  ig <- g
}

# Step 2: 获取每个节点的模块 ID 向量
mod_ids <- vertex_attr(ig, "module_id")

# Step 3: 对每个节点，检查其邻居是否来自不同模块
is_bridge <- sapply(V(ig), function(v) {
  nbrs <- neighbors(ig, v, mode = "all")
  if (length(nbrs) == 0) return(FALSE)
  unique_mods <- unique(mod_ids[nbrs])
  length(unique_mods) > 1
})

# Step 4: 将结果加回 tbl_graph
top_module_g_with_bridge <- g %>%
  activate(nodes) %>%
  mutate(is_bridge = is_bridge)

bridge_node_names <- top_module_g_with_bridge %>%
  activate(nodes) %>%
  as_tibble() %>%
  filter(is_bridge) %>%
  pull(name)

# top_module_p <- ggraph(top_module_g_tbl, layout = "kk") +
#   geom_edge_link(
#     color = "grey50",
#     alpha = 0.4,
#     width = 0.4
#   ) +
#   # 为每个模块添加背景区域（去边框）
#   ggforce::geom_mark_hull(
#     aes(x, y,
#         group = module_id,
#         fill = module_id
#     ),
#     alpha = 0.1,
#     concavity = 6,
#     expand = unit(2, "mm"),
#     radius = unit(3, "mm"),
#     lwd = 0.05,
#     color = 'grey90'
#   ) +
#   geom_node_point(
#     aes(color = node_color, size = centrality_degree()),
#     alpha = 0.9
#   ) +
#   # 👇 第二层：仅桥接节点，添加外边框（空心圆）
#   geom_node_point(
#     data = function(layout_data) {
#       # layout_data 是 ggraph 自动生成的节点布局数据（含 x, y, name 等）
#       layout_data %>%
#         dplyr::filter(name %in% bridge_node_names)
#     },
#     aes(x = x, y = y),
#     shape = 1,          # 空心圆（R 中 shape=1 是圆圈无填充）
#     color = "red",      # 边框颜色
#     stroke = 2,         # 边框粗细（关键！控制“外边框”视觉强度）
#     size = 6            # 可略大于底层点，确保覆盖
#   ) +
#   ggrepel::geom_text_repel(
#     aes(x = x, y = y, label = ifelse(name %in% bridge_node_names, name, NA)),
#     na.rm = TRUE,               # 忽略 NA 标签（即非 top6 节点）
#     size = 4,
#     color = "grey30",
#     fontface = "bold",
#     family = "sans",
#     box.padding = 1,
#     point.padding = 0.5,
#     segment.color = "grey50",
#     segment.size = 0.5,
#     max.overlaps = 30,
#     force = 2,
#     show.legend = FALSE
#   )+
#   # 设置颜色 scale - 修改 labels 参数
#   scale_color_manual(
#     name = "Module ID (Item Num.)",
#     values = top_module_full_color_mapping,
#     drop = FALSE,
#     labels = function(x) {
#       sapply(x, function(id) {
#         if (id %in% names(module_label_mapping)) {
#           module_label_mapping[as.character(id)]
#         } else {
#           id
#         }
#       })
#     }
#   ) +
#   scale_fill_manual(values = full_color_mapping, na.value = "gray70",
#                     guide = "none") +
#   scale_size_continuous(name = "Nodes Centrality",
#                         range = c(4, 8)) +
#   theme_void() +
#   theme(strip.text = element_text(face = "bold"),
#         plot.background = element_rect(fill = "white", color = NA),
#         # 控制图例字体大小
#         legend.title = element_text(size = 12, face = "bold"),      # 图例标题字体大小
#         legend.text = element_text(size = 12),                     # 图例文字字体大小
#         plot.title = element_text(
#           hjust = 0.5,           # 水平居中
#           face = "bold",         # 加粗
#           size = 20,             # 加大字体
#           margin = margin(t = 15, b = 15)  # 上下边距各15个单位
#         )) +
#   labs(title = paste0("Differential Protein and Phosphorylation Sites Network with top ",top_module_num," Modules"))

# 先确保 bridge_node_names 是字符向量（即使为空）
bridge_node_names <- bridge_node_names[!is.na(bridge_node_names)]  # 清理 NA（如有）

# 构建基础图
top_module_p <- ggraph(top_module_g_tbl, layout = "kk") +
  geom_edge_link(
    color = "grey50",
    alpha = 0.4,
    width = 0.4
  ) +
  ggforce::geom_mark_hull(
    aes(x, y, group = module_id, fill = module_id),
    alpha = 0.1,
    concavity = 6,
    expand = unit(2, "mm"),
    radius = unit(3, "mm"),
    lwd = 0.05,
    color = 'grey90'
  ) +
  geom_node_point(
    aes(color = node_color, size = centrality_degree()),
    alpha = 0.9
  )

#  条件添加：桥接节点外边框（仅当存在桥接节点时）
# if (length(bridge_node_names) > 0) {
#   top_module_p <- top_module_p +
#     geom_node_point(
#       data = function(layout_data) {
#         layout_data %>%
#           dplyr::filter(name %in% bridge_node_names)
#       },
#       aes(x = x, y = y),
#       shape = 1,
#       color = "grey30",
#       stroke = 2,
#       size = 6
#     )
# }
# 条件添加：桥接节点外边框 + 图例（仅当存在桥接节点时）
if (length(bridge_node_names) > 0) {
  top_module_p <- top_module_p +
    geom_node_point(
      data = function(layout_data) {
        layout_data %>%
          dplyr::filter(name %in% bridge_node_names)
      },
      aes(x = x, y = y, shape = "Bridge Node"),   # 映射 shape 到固定标签
      color = "grey30",
      stroke = 2,
      size = 6
    ) +
    # 手动定义 shape 图例（注意：不能放在 if 外面，否则可能报错）
    scale_shape_manual(
      name = "Node Type",
      values = c("Bridge Node" = 1),   # shape=1 是空心圆
      labels = c("Bridge Node" = "Bridge Node"),
      guide = guide_legend(override.aes = list(color = "grey30", stroke = 2, size = 6))
    )
}

#条件添加：桥接节点标签（仅当存在时）
if (length(bridge_node_names) > 0) {
  top_module_p <- top_module_p +
    ggrepel::geom_text_repel(
      aes(x = x, y = y, label = ifelse(name %in% bridge_node_names, name, NA)),
      na.rm = TRUE,
      size = 4,
      color = "grey30",
      fontface = "bold",
      family = "sans",
      box.padding = 1,
      point.padding = 0.5,
      segment.color = "grey50",
      segment.size = 0.5,
      max.overlaps = 30,
      force = 2,
      show.legend = FALSE
    )
}

# 添加 scales 和 theme（始终需要）
top_module_p <- top_module_p +
  scale_color_manual(
    name = "Module ID (Item Num.)",
    values = top_module_full_color_mapping,
    drop = FALSE,
    labels = function(x) {
      sapply(x, function(id) {
        if (id %in% names(module_label_mapping)) {
          module_label_mapping[as.character(id)]
        } else {
          id
        }
      })
    }
  ) +
  scale_fill_manual(values = full_color_mapping, na.value = "gray70", guide = "none") +
  scale_size_continuous(name = "Nodes Centrality", range = c(4, 8)) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20,
      margin = margin(t = 15, b = 15)
    )
  ) +
  labs(title = paste0("Differential Protein and Phosphorylation Sites Network with top ", top_module_num, " Modules"))

# top_module_p
ggsave("network_by_top_modules.png", plot = top_module_p, width = 12, height =10, dpi = 300)
################
###############
# 1. 首先获取所有模块的名称
module_labels <- top_module_g_tbl %>%
  activate(nodes) %>%
  as_tibble() %>%
  distinct(module_label) %>%
  pull(module_label)
# print(paste("总共有", length(module_labels), "个模块"))
# print(module_labels)
# 2. 创建一个空的列表来存储每个模块的顶部节点
top_nodes_list <- list()
# 3. 循环处理每个模块
for (module in module_labels) {
  # cat("\n正在处理模块:", module, "\n")
  
  # 3.1 提取当前模块的子图
  module_subgraph <- top_module_g_tbl %>%
    activate(nodes) %>%
    filter(module_label == module)
  
  # 检查子图是否有节点
  n_nodes <- module_subgraph %>%
    activate(nodes) %>%
    as_tibble() %>%
    nrow()
  
  # cat("模块", module, "有", n_nodes, "个节点\n")
  
  if (n_nodes > 0) {
    # 3.2 计算子图内的度中心性（只考虑模块内连接）
    module_nodes <- module_subgraph %>%
      activate(nodes) %>%
      mutate(
        module_degree = centrality_degree()  # 这是在子图内计算的度中心性
      ) %>%
      as_tibble() %>%
      arrange(desc(module_degree)) %>%
      select(name, module_label, module_degree)
    
    # 3.3 选择前3个节点（如果节点数少于3，则选择所有）
    n_to_select <- min(3, nrow(module_nodes))
    top_nodes_module <- module_nodes %>%
      slice_head(n = n_to_select) %>%
      mutate(
        rank_in_module = 1:n_to_select  # 添加在模块内的排名
      )
    
    # cat("选择的节点数:", n_to_select, "\n")
    # print(top_nodes_module)
    
    # 3.4 保存到列表
    top_nodes_list[[module]] <- top_nodes_module
  }
}

# 4. 合并所有模块的结果
top_nodes_per_module <- bind_rows(top_nodes_list) %>%
  arrange(module_label, desc(module_degree)) %>%
  pull(name)
# # 5. 查看结果
# cat("\n=== 所有模块的顶部节点汇总 ===\n")
# print(top_nodes_per_module)



top_module_p3 <- ggraph(top_module_g_tbl, layout = "kk") +
  geom_edge_link(color = "grey30", alpha = 0.3) +
  geom_node_point(aes(color = module_id,size = centrality_degree())) +
  facet_nodes(~ module_label, scales = "free") +
  scale_color_manual(name = "Module ID (Item Num.)",
                     values = full_color_mapping,
                     labels = function(x) {
                       sapply(x, function(id) {
                         if (id %in% names(module_label_mapping)) {
                           module_label_mapping[as.character(id)]
                         } else {
                           id
                         }
                       })
                     },
                     na.value = "gray80") +
  scale_size_continuous(name = "Nodes Centrality",  # 修改图例名称
                        range = c(4, 8)) +
  # # 添加高亮节点的标签
  # geom_node_text(
  #   data = top_nodes_per_module,
  #   aes(label = name),
  #   color = "black", 
  #   size = 4,
  #   fontface = "bold"
  # ) +
  ggrepel::geom_text_repel(
    aes(x = x, y = y, label = ifelse(name %in% top_nodes_per_module, name, NA)),
    na.rm = TRUE,
    size = 6,
    color = "grey30",
    fontface = "bold",
    family = "sans",
    box.padding = 1,
    point.padding = 0.5,
    segment.color = "grey50",
    segment.size = 0.5,
    max.overlaps = 30,
    force = 2,
    show.legend = FALSE
  )+
  theme_void() +
  theme( strip.text = element_text(
    face = "bold",
    size = 20,  # 调节小图标题大小，可以根据需要调整
    margin = margin(b = 5)  # 可选：调节标题与图之间的间距
  ),
        plot.background = element_rect(fill = "white", color = NA),
  # 控制图例字体大小
  legend.title = element_text(size = 20, face = "bold"),      # 图例标题字体大小
  legend.text = element_text(size = 20),                     # 图例文字字体大小
        
        plot.title = element_text(
          hjust = 0.5,           # 水平居中
          face = "bold",         # 加粗
          size = 30,             # 加大字体
          margin = margin(t = 15, b = 15)  # 上下边距各15个单位
        )) +
  labs(title = paste0("Differential Protein and Phosphorylation Sites Network by top ",top_module_num," Modules"))
# top_module_p3
# 保存图形
ggsave("network_by_separate_top_modules.png", plot = top_module_p3, width = 18, height =16, dpi = 300)

