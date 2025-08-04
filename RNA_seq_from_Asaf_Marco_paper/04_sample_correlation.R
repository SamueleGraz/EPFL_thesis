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
