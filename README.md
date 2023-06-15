# Differential Expression Analysis using DESeq

This repository utilizes the DESeq package in R to perform differential gene expression analysis on RNA-Seq data. It contains materials, scripts, and reflections from my internship at the Bioinformatics Department of the University of California Riverside (UCR) in 2012. The internship focused on differential gene expression analysis using R and Bioconductor. A week-by-week detailed breakdown of activities and learnings from each week of the internship is provided in the repository.

## Table of Contents

1. [Codes](#codes)
2. [Week-by-Week Summary](#week-by-week-summary)
3. [Key Learning Resources](#key-learning-resources)
4. [Acknowledgements](#acknowledgements)

### Codes
A select collection of scripts written during the internship: 

- Demultiplexing Sequencing Data
- Counting Reads and Assigning to Genes
- Differential Expression Analysis with DESeq
- Generating a Heatmap
- Constructing a Sample Tree

### Week-by-Week Summary
Detailed breakdown of activities and learnings from each week of the internship.

### Key Learning Resources

Over the course of the internship, I extensively used the following resources from the UCR Bioinformatics tutorials page prepared by Dr. Thomas Girke and his group at [https://girke.bioinformatics.ucr.edu/manuals/mydoc/home.html](https://girke.bioinformatics.ucr.edu/manuals/mydoc/home.html):

- R Basics Manual
- Bioconductor Manual
- NGS Analysis Workflows with systemPipeR
- NGS Analysis with R/Bioconductor
- Programming in R

#### 1. Demultiplexing Sequencing Data

In the first step, we separate multiplexed biological samples that have been tagged with unique barcodes during the sequencing process. 

Script: `demultiplexing.r`

#### 2. Counting Reads and Assigning to Genes

In the next step, we count the number of reads and assign them to genes. 

Script: `counting.r`

#### 3. Differential Expression Analysis with DESeq

With the count data prepared, we perform differential expression analysis using the DESeq package. 

Script: `deseq_analysis.r`

#### 4. Generating a Heatmap

To visualize the results, we generate a heatmap. 

Script: `heatmap.r`

#### 5. Constructing a Sample Tree

Lastly, we create a sample tree to visualize the relationships between our samples.

Script: `sample_tree.r`

Before running the scripts, make sure to set your working directory to the directory containing the scripts and data. Also, ensure all necessary packages are installed.

Please adjust file paths and parameters as needed for your specific data.

This is a simplified pipeline for differential expression analysis and does not include all steps typically involved in an RNA-Seq analysis pipeline, such as quality control, alignment, and gene annotation steps.

Note: Each section refers to a hypothetical script file (`demultiplexing.r`, `counting.r`, `deseq_analysis.r`, `heatmap.r`, `sample_tree.r`). You may need to adjust these based on how your actual code is organized.

## Acknowledgements

I am deeply grateful to Dr. Thomas Girke for the opportunity to have completed this internship at UC Riverside. The skills and knowledge I acquired during this time have significantly contributed to my development as a computational scientist.

