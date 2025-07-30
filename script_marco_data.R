# load package DESeq2 (all functions)
library(ggplot2)
library(DESeq2)

# read in the sample sheet
# header = TRUE: the first row is the "header", i.e. it contains the column names.
# sep = "\t": the columns/fields are separated with tabs.
marco_samplesheet <- read.table("samplesheet.txt", header=T, sep=" ")

# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(marco_samplesheet) <- marco_samplesheet$SampleName

# display the first 6 rows
head(marco_samplesheet)

# check the number of rows and the number of columns
nrow(marco_samplesheet) 
ncol(marco_samplesheet)

# Reads in a matrix:

# first read in the matrix
marco_matrix <- read.delim("final_unstranded_counts.tsv", header=T, sep="\t", row.names=1)
ncol(marco_matrix)

#check if the ncol and the nrow of marco_matrix and marco_samplesheet corresponds
all(colnames(marco_matrix) == rownames(marco_samplesheet))

# then create the DESeq object
# countData is the matrix containing the counts
# marco_samplesheet is the sample sheet / metadata we created
# design is how we wish to model the data: what we want to measure here is the difference between the treatment times
marco_se_star_matrix <- DESeqDataSetFromMatrix(countData = marco_matrix,
                                         colData = marco_samplesheet,
                                         design = ~ Time)
#filtering process
#We choose to keep only those genes that have more than 10 summed raw counts across the 12 samples

# Number of genes before filtering:
nrow(marco_se_star_matrix)

# Filter
marco_se_star_matrix <- marco_se_star_matrix[rowSums(counts(marco_se_star_matrix)) > 10, ]

# Number of genes left after low-count filtering:
nrow(marco_se_star_matrix)

#from ~78k to ~37k genes

#Fit statistical model
marco_se_star2_matrix <- DESeq(marco_se_star_matrix)

#count normalisation

# Read in the two-column data.frame linking transcript id (column 1) to gene id (column 2)
tx2gene <- read.table("tx2gene.gencode.v29.csv", 
                      sep="\t",
                      header=F)
head(tx2gene)
# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
marco_norm_counts <- log2(counts(marco_se_star2_matrix, normalized = TRUE)+1)

# add the gene symbols
marco_norm_counts_symbols <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(marco_norm_counts), marco_norm_counts), by=1, all=F)

colnames(marco_norm_counts_symbols)[1:2] <- c("gene_ID", "gene_name")

View(marco_norm_counts_symbols)
nrow(marco_norm_counts_symbols)
# write normalized counts to text file
write.table(marco_norm_counts_symbols, "marco_normalized_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")


#sample correlation

# load libraries pheatmap to create the heatmap plot
library(pheatmap)

#vst transformation
vsd <- vst(marco_se_star2_matrix)

# calculate between-sample distance matrix
marco_sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))

# create figure in PNG format
png("sample_distance_heatmap_star.png")
pheatmap(marco_sampleDistMatrix)
# close PNG file after writing figure in it
dev.off() 


#PCA
png("PCA_star.png")
plotPCA(object = vsd,
        intgroup = "Time")
dev.off()


#Differential expression analysis
# check results names: depends on what was modeled. Here it was the "Time"
resultsNames(marco_se_star2_matrix)

# extract results for early vs basal
# contrast: the column from the metadata that is used for the grouping of the samples (Time), then the baseline (basal) and the group compared to the baseline (early) -> results will be as "early vs basal"
de_early_vs_basal <- results(object = marco_se_star2_matrix, 
              name="Time_early_vs_basal")

# extract results for late vs basal
de_late_vs_basal <- results(object = marco_se_star2_matrix, 
                             name="Time_late_vs_basal")

# extract results for reactivated vs basal
de_reactivated_vs_basal <- results(object = marco_se_star2_matrix, 
                             name="Time_reactivated_vs_basal")

#visualisation
head(de_early_vs_basal[order(de_early_vs_basal$pvalue), ])
head(de_late_vs_basal[order(de_late_vs_basal$pvalue), ])
head(de_reactivated_vs_basal[order(de_reactivated_vs_basal$pvalue), ])

#generate more accurate log2 foldchange estimates
lfc_de_early_vs_basal <- lfcShrink(dds = marco_se_star2_matrix,
                       coef="Time_early_vs_basal",
                       type="apeglm")
#visualisation
nrow(lfc_de_early_vs_basal)

#I don't see any difference

# add the more comprehensive gene symbols DE tables
de_early_vs_basal_gene <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_early_vs_basal), de_early_vs_basal), by=1, all=F)
de_late_vs_basal_gene <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_late_vs_basal), de_late_vs_basal), by=1, all=F)
de_reactivated_vs_basal_gene <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_reactivated_vs_basal), de_reactivated_vs_basal), by=1, all=F)

# write differential expression analysis result to a text file
write.table(de_early_vs_basal_gene, "de_early_vs_basal_gene.txt", quote=F, col.names=T, row.names=F, sep="\t")
write.table(de_late_vs_basal_gene, "de_late_vs_basal_gene.txt", quote=F, col.names=T, row.names=F, sep="\t")
write.table(de_reactivated_vs_basal_gene, "de_reactivated_vs_basal_gene.txt", quote=F, col.names=T, row.names=F, sep="\t")

#load the table with the DE genes 
DE_React_Vs_Basal <- read.table("de_reactivated_vs_basal_padj0.05_log2fc0.5_gene.txt")
DE_Late_Vs_Basal <- read.table("de_late_vs_basal_padj0.05_log2fc0.5_gene.txt")
DE_Early_Vs_Basal <- read.table("de_early_vs_basal_padj0.05_log2fc0.5_gene.txt")

nrow(DE_React_Vs_Basal)
nrow(DE_Late_Vs_Basal)
nrow(DE_Early_Vs_Basal)

#DEG plot
#vulcano plot

colnames(de_early_vs_basal_gene)[2] <- "Gene"
colnames(de_late_vs_basal_gene)[2] <- "Gene"
colnames(de_reactivated_vs_basal_gene)[2] <- "Gene"

library(EnhancedVolcano)

Volc_early_vs_basal <- EnhancedVolcano(de_early_vs_basal_gene,
                                       lab = de_early_vs_basal_gene$Gene,
                                       x = 'log2FoldChange',
                                       y = 'padj',
                                       title = 'Early vs Basal',
                                       pCutoff = 0.01,
                                       FCcutoff = 2,
                                       pointSize = 2.0,
                                       labSize = 3.5,
                                       selectLab = NULL,               # No automatic labels
                                       col = c('grey30', 'grey30', 'grey30', 'red2')  # Only the significant ones in red
)

Volc_early_vs_basal
