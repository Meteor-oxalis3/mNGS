library(DESeq2)
library(EnhancedVolcano)
library(data.table)

# 1. 读取数据
count_data <- fread("KO_counts_matrix.tsv", data.table = FALSE, check.names=FALSE)
rownames(count_data) <- count_data[,1]
count_data <- count_data[,-1]

# 2. 读取样本信息
sample_info <- fread("sample_data.tsv", data.table = FALSE)
rownames(sample_info) <- sample_info$SampleID

# 3. 保留确诊组和无症状组样本
keep_samples <- sample_info$Group %in% c("Definite", "No Evidence")
sample_info <- sample_info[keep_samples, ]
count_data <- count_data[, rownames(sample_info)]

# 4. 创建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ Group)

# 5. 过滤低表达基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 6. 差异分析
dds <- DESeq(dds)

res <- results(dds, contrast = c("Group", "Definite", "No Evidence"))
res <- lfcShrink(dds, coef=2, type="ashr")  # 收缩log2FoldChange

# 7. 保存差异分析结果
write.csv(as.data.frame(res), file="DESeq2_results.csv")

pdf("volcano.pdf")
# 8. 火山图绘制
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'DESeq2 Volcano Plot: Definite vs No Evidence',
                subtitle = 'KO Abundance',
                caption = 'log2FC >1 & padj <0.05',
                colAlpha = 0.8,
                legendLabels=c('NS','Log2FC','p-value','p-value & Log2FC'))
dev.off()
