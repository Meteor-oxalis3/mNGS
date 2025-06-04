library(data.table)
library(tidyverse)

# 设置目录路径
input_dir <- "01_raw_Kraken2_reports"
file_paths <- list.files(input_dir, pattern = "\\.report$", full.names = TRUE)

# Step 1: 读取所有文件为长格式表格（taxa, count, sample）
data_list <- lapply(file_paths, function(file) {
  sample_name <- tools::file_path_sans_ext(basename(file))
  df <- fread(file, header = FALSE, sep = "\t", fill = TRUE)
  df <- df[!grepl("^#", V1) & V1 != "" & !is.na(V2)]  # 过滤注释和空行
  df$V1 <- gsub(" ", "", df$V1)  # 移除 taxonomy 中的空格
  df$sample <- sample_name
  setnames(df, c("taxonomy", "count", "sample"))
  df$count <- as.numeric(df$count)
  df
})

# Step 2: 合并所有数据
merged_long <- rbindlist(data_list, use.names = TRUE, fill = TRUE)

# Step 3: 转为 wide 矩阵，缺失填 0
merged_wide <- dcast(merged_long, taxonomy ~ sample, value.var = "count", fill = 0.0)

# Step 4: 保存结果
fwrite(merged_wide, "merged_taxonomy_matrix.tsv", sep = "\t", quote = FALSE)

