# R Scripts for RNA-Seq Data Analysis in Late Blight Disease
# Author: Seydanur Reisli (nee. Tikir)
# Contact: seydareisli@gmail.com


# This section of the script is performing read counting using the GenomicFeatures and Rsamtools packages from Bioconductor. The goal is to find out how many sequencing reads align to each feature of interest (like a gene or exon) in the genome. These counts can then be used for various types of downstream analyses, like identifying differentially expressed genes.



# Import the necessary libraries
library(Rsamtools); library(rtracklayer); library(GenomicFeatures)

# Define the paths to the necessary files
folders_tomato <- list.files(path = "/home/rsun/bigdata/RNA-seq/Carol/result4/", pattern = "Tomato.tophat$")

# "accepted_hits.bam" is a file that typically contains aligned sequencing data
bamFolders <- list.files(path = folders_tomato, pattern = "accepted_hits.bam$")

# Import the .gff3 file, which contains annotations of the genomic features of interest 
GFF <- import.gff3("deneme.gff3", asRangedData=FALSE, version="3")

# Create a list of BamFiles (binary alignment map format files, which store the genomic sequence data)
bamList <- BamFileList(bamFolders)

# Calculate read counts that overlap with the genomic features of interest
# The 'summarizeOverlaps' function does this counting for us
# 'singleEnd = TRUE' indicates that the sequencing data is single-end, as opposed to paired-end
bamData <- assays(summarizeOverlaps(GFF, bamList, singleEnd=TRUE))$counts

# Extract feature IDs and sequence names from the GFF data for further use
IDs <- elementMetadata(GFF)[["ID"]]
seqnames  <-as.character(seqnames(GFF))

# Combine all the data into a table
mytable <- data.frame(seqnames, IDs, bamData)

# Cleanup the bam folder names
my_names <-  gsub("/home/rsun/bigdata/RNA-seq/Carol/result4//", "", bamFolders)
my_names <-  gsub('.tophat/accepted_hits.bam', '' ,my_names)

# Write the final count table to a file, making it easier for further downstream analysis
write.table(mytable, "Tomato_count_table.xls",col.names=c("Seq_Name","Gene_ID",my_names),row.names=F ,sep="\t",quote=FALSE)