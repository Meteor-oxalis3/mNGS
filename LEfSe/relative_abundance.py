import pandas as pd

# 读取原始绝对丰度矩阵
df_absolute = pd.read_csv("merged_taxonomy_matrix.tsv", sep='\t', index_col=0)

# 计算相对丰度（每列除以该列总和）
df_relative = df_absolute.div(df_absolute.sum(axis=0), axis=1)

# 可选：保留小数位
df_relative = df_relative.round(6)

# 保存结果
df_relative.to_csv("taxonomy_relative_abundance.tsv", sep='\t')

print("✅ 相对丰度矩阵已保存为 taxonomy_relative_abundance.tsv")

