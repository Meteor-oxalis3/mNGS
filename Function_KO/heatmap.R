
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)

# 读取 KO 相对丰度数据
ko_data <- read.table("KO_relative_abundance.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# 读取样本分组信息
sample_info <- read.table("sample_data.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 统一列名格式（去掉.1后缀）
colnames(ko_data) <- gsub("\\.1$", "", colnames(ko_data))

# 保留 sample_data 中存在的样本列
ko_data <- ko_data[, sample_info$SampleID]

# 选前10个变异度最大的KO条目
ko_top10 <- ko_data %>%
  mutate(sd = apply(., 1, sd)) %>%
  arrange(desc(sd)) %>%
  head(10) %>%
  select(-sd)

# 转为矩阵
data_matrix <- as.matrix(ko_top10)

# 分组信息匹配为向量
group_vector <- sample_info$Group
names(group_vector) <- sample_info$SampleID

# 设置分组颜色（你指定的）
group_colors <- c("Water" = "#1f78b4", "Definite" = "red", "No Evidence" = "green")

# 构建列注释
col_anno <- HeatmapAnnotation(
  Group = group_vector,
  col = list(Group = group_colors),
  show_annotation_name = TRUE
)

# 热图颜色
heatmap_colors <- colorRamp2(
  seq(min(data_matrix), max(data_matrix), length = 9),
  colorRampPalette(brewer.pal(9, "YlOrRd"))(9)
)

pdf("KO_heatmap.pdf", width=15, height=3)
# 绘制热图
Heatmap(data_matrix,
        name = "Relative Abundance",
        top_annotation = col_anno,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
	show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 8),
        col = heatmap_colors
)
dev.off()

