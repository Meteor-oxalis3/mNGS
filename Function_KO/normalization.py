import pandas as pd

ko_counts = pd.read_csv('KO_counts_matrix.tsv', sep='\t', index_col=0)
rel_abundance = ko_counts.div(ko_counts.sum(axis=0), axis=1) * 1e6  # 转换为每百万丰度（可选）

rel_abundance.to_csv('KO_relative_abundance.tsv', sep='\t')

