# load package DESeq2 (all functions)
library(ggplot2)
library(DESeq2)
library(readxl)
library(dplyr)

#I filtered the contrasts applying padj lower than 0.01 and log2fc higher than 2
#like what it was done in the marco's paper
#I discover that they compared the contrast between early_basal, early_late and late_reactivated
#So I need to find the contrast between those conditions as well

# read in the sample sheet
# header = TRUE: the first row is the "header", i.e. it contains the column names.
# sep = "\t": the columns/fields are separated with tabs.
marco_samplesheet <- read.table("samplesheet.txt", header=T, sep=" ")

# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(marco_samplesheet) <- marco_samplesheet$SampleName

# Reads in a matrix:
# first read in the matrix
marco_matrix <- read.delim("final_unstranded_counts.tsv", header=T, sep="\t", row.names=1)

#check if the ncol and the nrow of marco_matrix and marco_samplesheet corresponds
all(colnames(marco_matrix) == rownames(marco_samplesheet))

# then create the DESeq object
# countData is the matrix containing the counts
# marco_samplesheet is the sample sheet / metadata we created
# design is how we wish to model the data: what we want to measure here is the difference between the treatment times
marco_se_star_matrix <- DESeqDataSetFromMatrix(countData = marco_matrix,
                                               colData = marco_samplesheet,
                                               design = ~ Time)

# Filter
marco_se_star_matrix <- marco_se_star_matrix[rowSums(counts(marco_se_star_matrix)) > 10, ]

#Fit statistical model
marco_se_star2_matrix <- DESeq(marco_se_star_matrix)

#Differential expression analysis

# check results names: depends on what was modeled. Here it was the "Time"
resultsNames(marco_se_star2_matrix)

#I have to change these contrasts

# early vs late
de_early_vs_late <- results(marco_se_star2_matrix,
                             contrast = c("Time", "early", "late"))

# late vs reactivated
de_late_vs_reactivated <- results(marco_se_star2_matrix,
                                   contrast = c("Time", "late", "reactivated"))

# basal vs early (equivalente al nome già esistente, ma per coerenza)
de_basal_vs_early <- results(marco_se_star2_matrix,
                              contrast = c("Time", "basal", "early"))

# add the more comprehensive gene symbols DE tables
de_early_vs_late <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_early_vs_late), de_early_vs_late), by=1, all=F)
de_late_vs_reactivated <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_late_vs_reactivated), de_late_vs_reactivated), by=1, all=F)
de_basal_vs_early <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_basal_vs_early), de_basal_vs_early), by=1, all=F)

#filter with padj lower than 0.01 and log2fc higher that the module of 2
de_early_vs_late_filtered <- de_early_vs_late[!is.na(de_early_vs_late$padj) & abs(de_early_vs_late$log2FoldChange) >= 2 & de_early_vs_late$padj <= 0.01, ]
de_late_vs_reactivated_filtered <- de_late_vs_reactivated[!is.na(de_late_vs_reactivated$padj) & abs(de_late_vs_reactivated$log2FoldChange) > 2 & de_late_vs_reactivated$padj < 0.01, ]
de_basal_vs_early_filtered <- de_basal_vs_early[!is.na(de_basal_vs_early$padj) & abs(de_basal_vs_early$log2FoldChange) > 2 & de_basal_vs_early$padj < 0.01, ]

#load the table with the original DE analysis 
DE_Basal_Vs_Early_original <- read_excel("nRNA-seq_marco_basal_vs_early_original.xlsx")
#delete useless rows
DE_Basal_Vs_Early_original <- DE_Basal_Vs_Early_original[-c(1), ]

#load the table with the original DE analysis 
DE_Early_vs_Late_original <- read_excel("nRNA-seq_marco_basal_vs_early_original.xlsx", 
                                                     sheet = "Early_vs_Late")
DE_Late_vs_Reactivated_original <- read_excel("nRNA-seq_marco_basal_vs_early_original.xlsx", 
                                                     sheet = "Late_vs_Reactivated")
View(DE_Basal_Vs_Early_original)

#change name of the gene name in "Gene"
colnames(DE_Basal_Vs_Early_original)[1] <- "Gene"
colnames(DE_Early_vs_Late_original)[1] <- "Gene"
colnames(DE_Late_vs_Reactivated_original)[1] <- "Gene"
colnames(de_early_vs_late_filtered)[2] <- "Gene"
colnames(de_late_vs_reactivated_filtered)[2] <- "Gene"
colnames(de_basal_vs_early_filtered)[2] <- "Gene"

nrow(DE_Basal_Vs_Early_original)
nrow(de_basal_vs_early_filtered)
nrow(DE_Early_vs_Late_original)
nrow(de_early_vs_late_filtered)
nrow(DE_Late_vs_Reactivated_original)
nrow(de_late_vs_reactivated_filtered)

#Now we can join accordingly the tables in order to check how many genes are in common

common_basal_vs_early_df <- de_basal_vs_early_filtered %>%
  inner_join(DE_Basal_Vs_Early_original, by = "Gene")

nrow(common_basal_vs_early_df)


common_early_vs_late_df <- de_early_vs_late_filtered %>%
  inner_join(DE_Early_vs_Late_original, by = "Gene")

nrow(common_early_vs_late_df)

common_late_vs_reactivated_df <- de_late_vs_reactivated_filtered %>%
  inner_join(DE_Late_vs_Reactivated_original, by = "Gene")

nrow(common_late_vs_reactivated_df)

#nice Venn plot
library(ggvenn)

Basal_vs_Early_Venn <- ggvenn(list(
  "Basal_vs_Early_Original" = DE_Basal_Vs_Early_original$Gene,
  "Basal_vs_Early_Samu" = de_basal_vs_early_filtered$Gene
), fill_color = c("#FF9999", "#99CCFF"), set_name_size = 4)+
  ggtitle("Basal_vs_Early") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 10))

Basal_vs_Early_Venn

Early_vs_Late_Venn <- ggvenn(list(
  "Early_vs_Late_Original" = DE_Early_vs_Late_original$Gene,
  "Early_vs_Late_Samu" = de_early_vs_late_filtered$Gene
), fill_color = c("#FF9999", "#99CCFF"), set_name_size = 4)+
  ggtitle("Early_vs_Late") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 10))

Early_vs_Late_Venn

Late_vs_Reactivated_Venn <- ggvenn(list(
  "Late_vs_Reactivated_Original" = DE_Late_vs_Reactivated_original$Gene,
  "Late_vs_Reactivated_Samu" = de_late_vs_reactivated_filtered$Gene
), fill_color = c("#FF9999", "#99CCFF"), set_name_size = 3.5)+
  ggtitle("Late vs Reactivated") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 10))

Late_vs_Reactivated_Venn

#put together the images
library(patchwork)

combined <- (Basal_vs_Early_Venn + Early_vs_Late_Venn) / Late_vs_Reactivated_Venn +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(plot.margin = margin(10, 10, 10, 10))

combined



#I risultati sono molto strani, forse ho interpretato male la base su cui è fatto il contrast,
#devo provare a invertire le basi di ogni comparazione fatta
