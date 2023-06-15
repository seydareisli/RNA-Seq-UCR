# R Scripts for RNA-Seq Data Analysis in Late Blight Disease
# Author: Seydanur Reisli (nee. Tikir)
# Contact: seydareisli@gmail.com

#In this code, the ape library (which provides tools for evolutionary and phylogenetic analysis) is used to create a dendrogram, a tree-like diagram that shows the hierarchical clustering of your RNA-seq data.

#First, it reads in a table of RPKM values. RPKM stands for Reads Per Kilobase Million, a normalization method for RNA-seq data that allows for comparison across different samples.

#Next, it calculates the Spearman correlation for these RPKM values, which is a measure of the statistical dependence between the rankings of two variables. The Spearman correlation values are then saved to a file.

#The code then performs hierarchical clustering of the data based on these correlation distances and plots the resulting dendrogram. The plot is saved as a PDF file. The settings for plot.phylo are formatting the plot as a phylogenetic tree (type="p"), with edges colored in blue (edge.col=4), edges with width 2 (edge.width=2), node labels displayed (show.node.label=TRUE), and without extra margins around the plot (no.margin=TRUE).



# Importing the necessary library for this analysis
library(ape)

# Reading the RPKM (Reads Per Kilobase Million) count data file
countDFrpkm <- read.table("./results/tomato_rpkm.xls")

# Subsetting the columns based on the names present in the labels vector
countDFrpkm <- countDFrpkm[, names(labels)]

# Updating the column names with the labels
colnames(countDFrpkm) <- labels

# Computing the spearman correlation for the RPKM data
d <- cor(countDFrpkm, method="spearman")

# Writing the correlation results to a file
write.table(d, "results/correlationMA_MeWB_annotation_June_2012_bowtie2.xls", quote=FALSE, sep="\t", col.names = NA)

# Performing hierarchical clustering based on the correlation distances (1-correlation gives us distance)
hc <- hclust(dist(1-d))

# Creating a PDF file to save the resulting plot
pdf("results/sample_tree_tomato.pdf")

# Plotting the dendrogram as a phylogenetic tree with specific formatting options
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=2, show.node.label=TRUE, no.margin=TRUE)

# Closing the PDF device
dev.off()