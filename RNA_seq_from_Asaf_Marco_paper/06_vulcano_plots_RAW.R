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
