library(data.table)

# 读取文件
df <- fread("taxonomy_relative_abundance.tsv")

# 备份 taxonomy 列
taxonomy_col <- df$taxonomy

# 提取纯数据部分（不动 taxonomy）
data_matrix <- df[, -1, with = FALSE]

# 设置过滤阈值
min_presence_ratio <- 0.1
min_total_abundance <- 0.01
n_samples <- ncol(data_matrix)

# 计算过滤条件
keep_rows <- rowSums(data_matrix > 0) >= (n_samples * min_presence_ratio) &
             rowSums(data_matrix) >= min_total_abundance

# 应用过滤，同时保留原始 taxonomy 列
df_filtered <- cbind(taxonomy = taxonomy_col[keep_rows], data_matrix[keep_rows])

# 写出结果
fwrite(df_filtered, "filtered_taxonomy_relative_abundance.tsv", sep = "\t")

