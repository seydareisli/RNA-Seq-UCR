# R Scripts for RNA-Seq Data Analysis in Late Blight Disease
# Author: Seydanur Reisli (nee. Tikir)
# Contact: seydareisli@gmail.com

#This section of the script is performing differential gene expression analysis with DESeq2, asking the question "are there any genes that are expressed at significantly different levels between our control and treatment samples?"

#First, the script defines the paths to the BAM files (which contain the aligned sequencing data) and the .gff3 file (which contains the genomic features of interest). Then it defines a simple experimental design with two conditions: control and treatment.

#It then creates a DESeqDataSet object from this information. This is the main object that DESeq2 works with; it contains both the raw counts data and information about the experimental design.

#The DESeq function is then called on the DESeqDataSet object to perform the differential expression analysis, and the results are extracted with the results function. These results, which contain the estimated log2 fold changes and p-values for each gene, are then written to a file for further analysis.



# Import the necessary library
library(DESeq2)

# Define paths to the necessary files
folders_tomato <- list.files(path = "/home/rsun/bigdata/RNA-seq/Carol/result4/", pattern = "Tomato.tophat$")

# "accepted_hits.bam" is a file that typically contains aligned sequencing data
bamFolders <- list.files(path = folders_tomato, pattern = "accepted_hits.bam$")

# Import the .gff3 file, which contains annotations of the genomic features of interest
GFF <- import.gff3("deneme.gff3", asRangedData=FALSE, version="3")

# Define sample information - in this case, it's defining the condition for each sample (control or treatment)
colData <- DataFrame(condition = factor(c("control", "treatment")))

# Create a DESeqDataSet object from the BAM files and the sample information
dds <- DESeqDataSetFromBAM(bamFolders, colData = colData, referenceSequenceFile = "deneme.gff3")

# Perform the differential expression analysis
dds <- DESeq(dds)

# Get the results
res <- results(dds)

# Write the results to a file for further analysis
write.csv(as.data.frame(res), file="DESeq2_results.csv")