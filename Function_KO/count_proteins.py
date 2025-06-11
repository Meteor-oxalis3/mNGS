import pandas as pd 
from glob import glob
from collections import defaultdict

# 读取protein2KO映射
protein2ko = pd.read_csv('protein2KO.tsv', sep='\t', header=None, names=['protein','KO'])

# 组织所有 m8 文件
m8_files = glob('diamond_output/*.m8')

# 样本分组（按文件名前缀归组）
sample_files = defaultdict(list)
for f in m8_files:
    sample_name = f.split('/')[-1].split('_')[0]
    sample_files[sample_name].append(f)

count_matrix = pd.DataFrame()

for sample, files in sample_files.items():
    # 合并多个文件
    dfs = [pd.read_csv(f, sep='\t', header=None, usecols=[1], names=['protein']) for f in files]
    df_merged = pd.concat(dfs)
    counts = df_merged['protein'].value_counts()
    counts.name = sample
    count_matrix = pd.concat([count_matrix, counts], axis=1)

count_matrix.fillna(0, inplace=True)

# 连接蛋白到KO
count_matrix['KO'] = protein2ko.set_index('protein').reindex(count_matrix.index)['KO']

# 按KO聚合
ko_counts = count_matrix.groupby('KO').sum()

# 输出矩阵
ko_counts.to_csv('KO_counts_matrix.tsv', sep='\t')

