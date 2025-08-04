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

# add the more comprehensive gene symbols DE tables
de_early_vs_basal_gene <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_early_vs_basal), de_early_vs_basal), by=1, all=F)
de_late_vs_basal_gene <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_late_vs_basal), de_late_vs_basal), by=1, all=F)
de_reactivated_vs_basal_gene <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(de_reactivated_vs_basal), de_reactivated_vs_basal), by=1, all=F)

# write differential expression analysis result to a text file
write.table(de_early_vs_basal_gene, "de_early_vs_basal_gene.txt", quote=F, col.names=T, row.names=F, sep="\t")
write.table(de_late_vs_basal_gene, "de_late_vs_basal_gene.txt", quote=F, col.names=T, row.names=F, sep="\t")
write.table(de_reactivated_vs_basal_gene, "de_reactivated_vs_basal_gene.txt", quote=F, col.names=T, row.names=F, sep="\t")
