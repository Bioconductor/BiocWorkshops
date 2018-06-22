
# RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR

## Author: 
* Charity Law (<law@wehi.edu.au>)

## Workshop Description

In this instructor-led live demo, we analyse RNA-sequencing data from the mouse mammary gland, demonstrating use of the popular **edgeR** package to import, organise, filter and normalise the data, followed by the **limma** package with its voom method, linear modelling and empirical Bayes moderation to assess differential expression and graphical representations. This pipeline is further enhanced by the **Glimma** package which enables interactive exploration of the results so that individual samples and genes can be examined by the user. The complete analysis offered by these three packages highlights the ease with which researchers can turn the raw counts from an RNA-sequencing experiment into biological insights using Bioconductor. The complete workflow is available at http://master.bioconductor.org/packages/release/workflows/html/RNAseq123.html .

### Pre-requisites

* Basic knowledge of RNA-sequencing 
* Basic knowledge of R syntax, R object classes and object manipulation

### Workshop Participation

Participants can watch the live demo, or may prefer to follow the demonstration by bringing their laptops along. To follow the analysis on their own laptops, participants need to install the *RNAseq123* workflow by running


```r
source("https://bioconductor.org/biocLite.R")  
biocLite("RNAseq123")
```

in R. The relevant sequencing data should also be download in advance.


```r
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"  
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb")   
utils::untar("GSE63310_RAW.tar", exdir = ".")  
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))  
  R.utils::gunzip(i, overwrite=TRUE)  
```

Due to time restraints, extra help regarding R package installation and coding errors will not be addressed during the workshop. 


### _R_ / _Bioconductor_ packages used

Bioconductor: limma, Glimma, edgeR, Mus.musculus  
CRAN: RColorBrewer, gplots

### Time outline

| Activity                     | Time |
|------------------------------|------|
| Introduction                 | 5mins  |
| Data packaging               | 10mins  |
| Data pre-processing          | 15mins  |
| Differential expression analysis | 30mins   |

## Workshop goals and objectives

The key steps to RNA-seq data analysis are described in this workshop with basic statistical theory of methods used. The goal is to allow beginner-analysts of RNA-seq data to become familiar with each of the steps involved, as well as completing a standard analysis pipeline from start to finish. 

### Learning goals

* learn how to analyse RNA-seq data
* identify methods for pre-processing data
* understand linear models used in differential expression analysis
* examine plots for data exploration and result representation

### Learning objectives

* read in count data and format as a DGEList-object
* annotate Entrez gene identifiers with gene information
* filter out lowly expressed genes
* normalise gene expression values
* unsupervised clustering of samples (standard and interactive plots)
* linear modelling for comparisons of interest
* remove heteroscedascity 
* examine the number of differentially expressed genes
* mean-difference plots (standard and interactive plots)
* heatmaps
