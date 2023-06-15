# R Scripts for RNA-Seq Data Analysis in Late Blight Disease
# Author: Seydanur Reisli (nee. Tikir)
# Contact: seydareisli@gmail.com


# This code block handles the same functionality as the previous one (counting_v1.r) but in a slightly different approach.

library(Rsamtools); library(rtracklayer); library(GenomicFeatures);library(Streamer)

folders_tomato <- list.files(path = "/home/rsun/bigdata/RNA-seq/Carol/result4/", pattern = "Tomato.tophat$")
bamFolders <- list.files(path = folders_tomato, pattern = "accepted_hits.bam$")

GFF <- import.gff3("deneme.gff3", asRangedData=FALSE, version="3")
bamList <- BamFileList(bamFolders, yieldSize=25000)
trace(readBamGappedAlignment)
bamData <- assays(summarizeOverlaps(GFF, bamList, singleEnd=TRUE))$counts
untrace(readBamGappedAlignment)

IDs <- elementMetadata(GFF)[["ID"]]
seqnames  <-as.character(seqnames(GFF))
mytable <- data.frame(seqnames, IDs)
mytable<- cbind(mytable,bamData)

my_names <-  gsub("/home/rsun/bigdata/RNA-seq/Carol/result5//", "", bamFolders)
my_names <-  gsub('.tophat/accepted_hits.bam', '' ,my_names)
write.table(table, "Pinfestants_pair1_count_table.xls",col.names=F ,row.names=c("Seq_Name","Gene_ID",mynames),sep="\t",quote=FALSE)

