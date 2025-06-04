import pandas as pd

# 读取相对丰度矩阵
df = pd.read_csv("filtered_taxonomy_relative_abundance.tsv", sep='\t', index_col=0)

# 读取分组信息文件
def read_group(filename, group_name):
    with open(filename) as f:
        return {sample.strip(): group_name for sample in f if sample.strip()}

# 构建样本分组字典
group_dict = {}
group_dict.update(read_group("Definite_samples.txt", "Definite"))
group_dict.update(read_group("No_Evidence_samples.txt", "No_Evidence"))
# 如有其他组继续添加，例如：
# group_dict.update(read_group("Suspected.txt", "Suspected"))

# 构建 group 行（与 df 列顺序一致）
group_row = [group_dict.get(sample, "Water") for sample in df.columns]

# 将 group 行添加为第一行
df_with_group = pd.concat(
    [pd.DataFrame([group_row], columns=df.columns, index=["Group"]), df]
)

# 输出结果
df_with_group.to_csv("taxonomy_relative_abundance_with_group.tsv", sep='\t')

print("✅ 添加分组行后的相对丰度矩阵已保存为 taxonomy_relative_abundance_with_group.tsv")

