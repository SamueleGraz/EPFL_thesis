# load package DESeq2
library(ggplot2)
library(DESeq2)

# read in the sample sheet
marco_samplesheet <- read.table("samplesheet.txt", header=T, sep=" ")

# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(marco_samplesheet) <- marco_samplesheet$SampleName

# Reads in a matrix:
marco_matrix <- read.delim("final_unstranded_counts.tsv", header=T, sep="\t", row.names=1)

#check if the ncol and the nrow of marco_matrix and marco_samplesheet corresponds
all(colnames(marco_matrix) == rownames(marco_samplesheet))

# then create the DESeq object
marco_se_star_matrix <- DESeqDataSetFromMatrix(countData = marco_matrix,
                                         colData = marco_samplesheet,
                                         design = ~ Time)

#We choose to keep only those genes that have more than 10 summed raw counts across the 12 samples
marco_se_star_matrix <- marco_se_star_matrix[rowSums(counts(marco_se_star_matrix)) > 10, ]

#from ~78k to ~37k genes

#Fit statistical model
marco_se_star2_matrix <- DESeq(marco_se_star_matrix)

# Read in the two-column data.frame linking transcript id (column 1) to gene id (column 2)
tx2gene <- read.table("tx2gene.gencode.v29.csv", 
                      sep="\t",
                      header=F)

# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
marco_norm_counts <- log2(counts(marco_se_star2_matrix, normalized = TRUE)+1)

# add the gene symbols
marco_norm_counts_symbols <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(marco_norm_counts), marco_norm_counts), by=1, all=F)
colnames(marco_norm_counts_symbols)[1:2] <- c("gene_ID", "gene_name")

# write normalized counts to text file
write.table(marco_norm_counts_symbols, "marco_normalized_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")
