# R Scripts for RNA-Seq Data Analysis in Late Blight Disease
# Author: Seydanur Reisli (nee. Tikir)
# Contact: seydareisli@gmail.com


# In this part of the code, we are demultiplexing the sequencing data, which is the process of separating multiplexed (i.e., pooled) biological samples that have been tagged with unique barcodes during the sequencing process.



# 1. DEMULTIPLEXING 

library(ShortRead)
adaptors <-c("CTTGTAA","CGATGTA", "TGACCAA", "ACAGTGA", "CAGATCA")
indexread <- "lane1_NoIndex_L001_R2_001.fastq.gz"
fastqread <- "lane1_NoIndex_L001_R3_001.fastq.gz"

indexstream <- FastqStreamer(indexread, 50000)

while(length(index <- yield(indexstream))){
  fastq <- yield(FastqStreamer(fastqread, 50000))
  
  for(ad in adaptors){
    indexID <-which(as.character(sread(index))==ad)
    read <- fast2[indexID]
    newID <- paste(id(read),ad,sep="_")
    newID <- BStringSet(x=newID, start=NA, end=NA, width=NA, use.names=TRUE)
    newobject <- ShortReadQ(sread=sread(read), quality=quality(read), id=newID)
    writeFastq(newobject,paste(gsub(".fastq.gz","",fastqread), "_", ad, ".fastq",sep=""), mode='a')
  }
}

