library(phyloseq)
library(ggplot2)
library(readr)
library(tibble)

# 读取 OTU 表，第一列为 taxonomy 信息
otu_tax <- read_tsv("input_matrix.tsv")
taxonomy <- otu_tax$taxonomy
otu_table_data <- otu_tax[, -1]
rownames(otu_table_data) <- taxonomy

# 拆分 taxonomy 字符串为各层级
tax_split <- strsplit(taxonomy, ";")
max_ranks <- max(sapply(tax_split, length))
tax_mat <- do.call(rbind, lapply(tax_split, function(x) c(x, rep(NA, max_ranks - length(x)))))
rownames(tax_mat) <- taxonomy
colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:max_ranks]

# 构建 phyloseq 对象
otu <- otu_table(as.matrix(otu_table_data), taxa_are_rows = TRUE)
tax <- tax_table(as.matrix(tax_mat))
physeq <- phyloseq(otu, tax)

# 读取并合并样本信息
sample_data_df <- read_tsv("sample_data.tsv")
sample_data_df <- column_to_rownames(sample_data_df, var = "SampleID")
physeq <- merge_phyloseq(physeq, sample_data(sample_data_df))

# 筛选丰度最高的前 20 个 OTU
otu_mat <- as(otu_table(physeq), "matrix")
if (!taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
top20_taxa <- names(sort(rowSums(otu_mat), decreasing = TRUE))[1:20]
physeq_top20 <- prune_taxa(top20_taxa, physeq)

################################################################################
# (OK) 1. Alpha多样性分析
pdf("alpha_group.pdf", width = 18, height = 5)
plot_richness(physeq, x="Group") + 
  geom_boxplot(aes(fill=Group)) +
  theme_minimal()
dev.off()
################################################################################
# (OK) 2. Beta多样性分析
# 基于Bray-Curtis距离的PCoA
pdf("beta_group.pdf", width = 8, height = 6)
ordu_pcoa <- ordinate(physeq, method="PCoA", distance="bray")
plot_ordination(physeq, ordu_pcoa, color="Group") +
  geom_point(size=4) +
  stat_ellipse(aes(color=Group), type="t", linetype=2) +
  ggtitle("PCoA (Bray-Curtis) - Group") +
  theme_minimal()
dev.off()
###############################################################################
# 3. 柱状图（Barplot）展示分类单元组成
# 展示样本中不同分类单元（如属或门）的相对丰度。
# 转换为相对丰度
# 1. 聚合到 Species 层级并转换为相对丰度
physeq_species <- tax_glom(physeq, taxrank = "Species")
physeq_species_relabund <- transform_sample_counts(physeq_species, function(x) x / sum(x))
# 2. 计算所有物种的总相对丰度（在所有样本中求和）
species_sums <- taxa_sums(physeq_species_relabund)
# 3. 获取前 20 丰度物种的 OTU 名称
top20_species <- names(sort(species_sums, decreasing = TRUE))[1:20]
# 4. 筛选 phyloseq 对象，仅保留这前 20 个物种
physeq_top20_species <- prune_taxa(top20_species, physeq_species_relabund)
# 5. 绘图
pdf("stacked_bar_top20.pdf",width = 30, height = 6)
plot_bar(physeq_top20_species, fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
###############################################################################
# 4. 稀释曲线（Rarefaction curve）(OK!)
# 查看测序深度是否足够覆盖微生物多样性。
library(vegan)
pdf("rarefaction_curve.pdf", width = 10, height = 6)
rarecurve(t(otu_mat), step=100, cex=0.5)
dev.off()

###############################################################################
# 5. 热图（top20）
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 读取文件
df <- read.table("taxonomy_relative_abundance_with_group.tsv",
                 sep = "\t", header = TRUE, check.names = FALSE)

# 提取分组信息
groups <- as.character(unlist(df[1, -1]))  # 去掉第一列是"Group"
sample_names <- colnames(df)[-1]

# 创建分组注释数据框
group_df <- data.frame(Group = factor(groups, levels = c("Water", "Definite", "No_Evidence")))
rownames(group_df) <- sample_names

# 提取表达矩阵
data_matrix <- df[-1, ]  # 去掉第一行是 Group 信息
rownames(data_matrix) <- data_matrix[[1]]  # 第一列是物种名
data_matrix <- data_matrix[, -1]  # 去掉物种列
data_matrix <- as.matrix(data_matrix)
storage.mode(data_matrix) <- "numeric"

# 按 group 排列样本顺序
ordered_cols <- rownames(group_df)[order(group_df$Group)]
data_matrix <- data_matrix[, ordered_cols]
group_df <- group_df[ordered_cols, , drop = FALSE]

# 只保留行名中含有 "s__" 的行并简化行名
data_matrix <- data_matrix[grep("s__", rownames(data_matrix)), ]
rownames(data_matrix) <- sub(".*s__", "", rownames(data_matrix))

# 选择 top 20 丰度物种
topN <- 20
row_means <- rowMeans(data_matrix)
data_matrix <- data_matrix[order(row_means, decreasing = TRUE)[1:topN], ]

# 定义颜色
group_colors <- c("Water" = "#1f78b4", "Definite" = "red", "No_Evidence" = "green")
heatmap_colors <- colorRamp2(seq(min(data_matrix), max(data_matrix), length = 9),
                             colorRampPalette(brewer.pal(9, "YlOrRd"))(9))

# 创建列注释
ha_column <- HeatmapAnnotation(
  Group = group_df$Group,
  col = list(Group = group_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left"
)

# 创建热图
pdf("heatmap_complex.pdf", width =15, height = 5)
Heatmap(data_matrix,
        name = "Abundance",
        col = heatmap_colors,
        top_annotation = ha_column,
        row_names_side = "left",
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 10),
        cluster_columns = TRUE,  # 不聚类列以保持分组顺序
        cluster_rows = FALSE)     # 不聚类行以保持丰度顺序
dev.off()

