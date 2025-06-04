import pandas as pd

# 读相对丰度矩阵
df = pd.read_csv("taxonomy_relative_abundance_with_group.tsv", sep='\t', index_col=0,low_memory=False)

# 读性别信息
df_sex = pd.read_csv("sex.tsv", sep='\t', header=None, names=["sample_name", "sex"])

# 构造性别行（顺序按相对丰度矩阵列顺序）
sex_row = [df_sex.set_index('sample_name').loc[sample, 'sex'] if sample in df_sex['sample_name'].values else 'Unknown' for sample in df.columns]

# 把性别行添加到相对丰度矩阵第一行
df_with_sex = pd.concat([pd.DataFrame([sex_row], columns=df.columns, index=["Sex"]), df])

# 保存新矩阵
df_with_sex.to_csv("taxonomy_relative_abundance_with_sex.tsv", sep='\t')

print("✅ 添加性别信息后的矩阵已保存：taxonomy_relative_abundance_with_sex.tsv")