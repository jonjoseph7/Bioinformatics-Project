install.packages("readr",dependencies = TRUE)
warnings()
library(readr)
install.packages("factoextra")
library(factoextra)
library(cluster)

install.packages("pheatmap")
library(pheatmap)


install.packages("alluvial")
library(alluvial)


install.packages("ggplot2") 
install.packages("ggalluvial") 
library(ggplot2)
library(ggalluvial)

install.packages("dendextend")
library(dendextend)

install.packages("circlize")
library(circlize)


install.packages("riverplot")  # Install the riverplot package
library(riverplot)




data <- readr::read_tsv("C:\\Users\\jonjo\\Bioinformatics-Project\\Data\\SRP055573\\SRP055573.tsv")
data <- data.frame(data)
metadata <- readr::read_tsv("C:\\Users\\jonjo\\Bioinformatics-Project\\Data\\SRP055573\\metadata_SRP055573.tsv")
metadata <- data.frame(metadata)

gene_sd <- apply(subset(data, select = -Gene), 1, sd)
gene_sd

data$gene_sd <- gene_sd
sorted_data <- data[order(data$gene_sd, decreasing = TRUE), ]

top5000 <- head(sorted_data, n = 5000)
rownames(top5000) <- top5000$Gene
top5000 <- top5000[,2:69]
top5000 = scale(top5000)

fviz_nbclust(top5000, pam, method ="silhouette")+theme_minimal()
result <-pam(scaleddata, k = 3)
result
result$medoids
result$clusinfo
result$clustering
fviz_cluster(result, 
             palette =c("#FF1812","#D9FF5F","#ffc5cb"), #,"#1050FF"),#,"#ffa500"),
             ellipse.type ="euclid",
             repel =TRUE,
             ggtheme =theme_minimal())

cluster_assignments <- result$clustering

annotations <- data.frame(
  PAM = as.factor(cluster_list)  # Replace 'Clusters' with your cluster method name
)

# Create a color palette for clusters
annotation_colors <- list(
  PAM = colorRampPalette(c("#FF1812", "#D9FF5F", "#ffc5cb"))(length(unique(annotations$Clusters)))
)

# Create the heatmap
pheatmap(
  top5000,
  cluster_cols = FALSE,  # To exclude column clustering
  #annotation_col = annotations,
  #annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust the color scale as needed
  main = "Heatmap",
  legend = TRUE
)
