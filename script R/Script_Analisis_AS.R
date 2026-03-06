#Gene Expression Analysis in patients with Ankylosing Spondylitis (AS)

#Dataset : GSE25101 (Ankylosing Spondylitis vs Normal)
#Platform : Microarray (Illumina HumanHT-12 V3.0 - GLP6947)
#Purpose : Identifying Differentially Expressed Genes (DEG)


######Data collection from GEO
options(timeout = 1000)
options(download.file.method = 'libcurl')
library(GEOquery)
gset <- getGEO('GSE25101', GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]


######Check whether log2 transformation is needed
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}


######Definition of sample group
meta <- pData(gset)
colnames(meta)

ignore.case = TRUE 
group_info <- ifelse(
grepl('Ankylosing|AS|patient', meta$title, ignore.case = TRUE),
'AS',
'Control'
)   

groups <- factor(group_info)
table(groups)

gset$group <- factor(groups)

group_name <- levels(gset$group)
print(group_name)


######Create a design matrix
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

case_group <- 'AS'
control_group <- 'Control'

contrast_formula <- paste(case_group, '-', control_group)
print(paste('Analyzed contrast:', contrast_formula))


######Differential Expression Analysis (Limma)
library(limma)
fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

topTableResults <- topTable(
fit2,
adjust = 'fdr',
sort.by = 'B',
number = 'Inf',
p.value = 0.05
)
head(topTableResults)


######Gene name annotation
library(illuminaHumanv3.db)
probe_ids <- rownames(topTableResults)

gene_annotation <- AnnotationDbi::select(
illuminaHumanv3.db,
keys = probe_ids,
columns = c('SYMBOL', 'GENENAME'),
keytype = 'PROBEID'
)

topTableResults$PROBEID <-
rownames(topTableResults)
topTableResults <- merge(
topTableResults,
gene_annotation,
by = 'PROBEID',
all.x = TRUE
)

head(topTableResults[, c('PROBEID', 'SYMBOL', 'GENENAME')])


######Visualization using Boxplot
library(ggplot2)
library(tidyr)

significant_gene <- c('LSM3', 'TMA7', 'GNG11', 'UQCRBP1', 'COX7B', 'S100A12')

matches <- topTableResults[topTableResults$SYMBOL %in% significant_gene, ]
target_indices <- match(significant_gene, matches$SYMBOL)
target_ids <- matches$PROBEID[target_indices]

df_plot <- as.data.frame(t(ex[target_ids, ]))
colnames(df_plot) <- significant_gene
df_plot$Group <- gset$group

df_long <- pivot_longer(df_plot, cols = -Group, names_to = 'Gene', values_to = 'Expression')

ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
geom_jitter(width = 0.2, alpha = 0.4, size = 1.2) +
facet_wrap(~Gene, scales = 'free_y') + 
scale_fill_manual(values = c('AS' = 'firebrick3', 'Control' = 'navy')) +
theme_minimal() +
labs(title = 'Validation of Candidates Genes for Ankylosing Spondylitis',
subtitle = 'Verified Top Genes from Differential Expression Analysis',
y = 'Log2 Expression Value', x = '') +
theme(legend.position = 'none', 
strip.text = element_text (face = 'bold', size = 11))


######Visualization using Density Plot
library(ggplot2)
expr_long <- data.frame(
Expression = as.vector(ex),
Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, fill = Group, color = Group)) + 
geom_density(alpha = 0.2, linewidth = 1) +
scale_fill_manual(values = c('AS' = 'firebrick3', 'Control' = 'navy')) +
scale_color_manual(values = c('AS' = 'firebrick3', 'Control' = 'navy')) +
theme_minimal() +
labs(
title = 'Global Gene Expression Density : AS vs Control',
subtitle = 'Distribution of 18,000+ Genes in Whole Blood Samples',
x = 'Expression Value (log2)',
y = 'Density'
) +
theme(legend.position = 'top')


######Visualization using UMAP
library(umap)
library(ggplot2)

gene_vars <- apply(ex, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:1000]
umap_input <- t(ex[top_genes, ])

custom_config <- umap.defaults
custom_config$random_state <- 123
umap_result <- umap(umap_input, config = custom_config)
                    
umap_df <- data.frame(
UMAP1 = umap_result$layout[, 1],
UMAP2 = umap_result$layout[, 2],
Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
geom_point(size = 4, alpha = 0.7) +
scale_color_manual(values = c('AS' = 'firebrick3', 'Control' = 'navy')) +
theme_minimal() +
labs(
title = 'UMAP Plot : Global Gene Expression Profile',
subtitle = 'Based on Top 1000 Most Variable Genes',
x = 'UMAP Dimension 1',
y = 'UMAP Dimension 2',
caption = 'Source : Whole Blood Samples of AS Patients'
) +
theme(legend.position = 'bottom')


######Visualization using Volcano Plot
install.packages('ggrepel')
library(ggplot2)
library(ggrepel)
volcano_data <- data.frame(
logFC = topTableResults$logFC,
adj.P.Val = topTableResults$adj.P.Val,
Gene = topTableResults$SYMBOL
)

p_tresh <- 0.05
fc_tresh <- 1

volcano_data$status <- 'NO'
volcano_data$status[volcano_data$logFC >= fc_tresh & volcano_data$adj.P.Val < p_tresh] <- 'UP'
volcano_data$status[volcano_data$logFC <= -fc_tresh & volcano_data$adj.P.Val < p_tresh] <- 'DOWN'

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
geom_point(alpha = 0.6, size = 1.5) +
scale_color_manual(values = c('DOWN' = 'blue', 'NO' = 'grey', 'UP' = 'red')) +
geom_vline(xintercept = c(-fc_tresh, fc_tresh), linetype = 'dashed', color = 'black') +
geom_hline(yintercept = -log10(p_tresh), linetype = 'dashed', color = 'black') +

geom_text_repel(data = subset (volcano_data, status != 'NO'), aes(label = Gene),
size = 3, box.padding = 0.5) + 
theme_minimal () + 
labs(title = 'Differential Expression Genes : Ankylosing Spondylitis', 
subtitle = paste('Threshold: p <', p_tresh, '& |logFC| >', fc_tresh), 
x = 'log2 Fold Change', y = '-log10 Adjusted P-Value')


######Visualization using Heatmap
library(pheatmap)
library(RColorBrewer)
topTableResults <- topTableResults[order(topTableResults$adj.P.Val),]
top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$PROBEID,]

gene_label <- ifelse(
is.na(top50$SYMBOL)|top50$SYMBOL =='',
top50$PROBEID,
top50$SYMBOL
)
rownames(mat_heatmap) <- gene_label

mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0,]

annotation_col <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
mat_heatmap,
scale = 'row',
annotation_col = annotation_col, 
color = colorRampPalette(c('navy', 'white', 'firebrick3'))(50),
border_color = NA,
show_colnames = FALSE,
show_rownames = TRUE,
fontsize_row = 7,
clustering_distance_rows = 'euclidean',
clustering_distance_cols = 'euclidean',
clustering_method = 'complete',
main = 'Top 50 Differentially Expressed Genes of Ankylosing Spondylitis'
)


#######Save the results
write.csv(topTableResults, 'Results_GSE25101_DEG.csv')
message('Analysis complete. The result file has been saved.')

class(groups)
is.factor(groups)
levels(groups)