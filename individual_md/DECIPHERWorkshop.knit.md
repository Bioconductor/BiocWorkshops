# Working with Genomic Data in R with the DECIPHER package

Authors:
    Nicholas Cooley^[University of Pittsburgh]
Last Modified: 29 June, 2018
    
<!--
  bookdown::render_book("", "bookdown::gitbook")
-->

## Overview

### Workshop Description

In this workshop we will give an introduction to working with biological
sequence data in R using the Biostrings and DECIPHER packages. We will cover:

* Importing, viewing, and manipulating genomes in R
* Construction of sequence databases for organizing genomes
* Mapping syntenic regions between genomes with the FindSynteny function
* Understanding, accessing, and viewing objects of class Synteny
* Using syntenic information to predict orthologous genes
* Alignment of sequences and genomes with DECIPHER
* Downstream analyses enabled by syntenic mapping

### Pre-requisites

* Familiarity with Biostrings
* Familiarity with DECIPHER Databases [^1]

[^1]: [Wright, E. S. The R Journal 2016, 8 (1), 352–359.](https://journal.r-project.org/archive/2016-1/wright.pdf)

### Workshop Participation

This will be a lab where participants follow along on their computers.

### _R_ / _Bioconductor_ packages used

* Biostrings
* DECIPHER

### Time outline

| Activity                          | Time |
|-----------------------------------|------|
| Packages/Introduction             | 5m   |
| Basic commands for sequences      | 5m   |
| Sequence databases                | 10m  |
| Demonstration of synteny mapping  | 10m  |
| Explanation of function arguments | 10m  |
| Dissecting Synteny objects        | 10m  |
| Visualization of syntenic blocks  | 10m  |
| Alignment of syntenic regions     | 10m  |
| Ortholog prediction from Synteny  | 10m  |
| Constructing phylogenetic trees   | 10m  |

## Workshop goals and objectives

### Learning goals

* Understand a simple workflow for analysis of sequences in R and DECIPHER
* Learn the basic use and appropriate application of functions within DECIPHER

### Learning objectives

* Learn basic commands for working with sequences in R
* Import genomes from online repositories or local files
* Map synteny between genomes
* Analyze a synteny map among multiple genomes
* Develop an understanding of the data structures involved
* Predict orthologs from syntenic maps
* Select core and pan genomes from predicted orthologs
* Construct and interpret phylogenetic trees

### Workshop

In this worshop we will be walking through a comparative genomics pipeline using functions within the R package DECIPHER, and within a package that is currently under construction under the name FindHomology.

Lipid membranes are ubiquitious in living organisms, however the individual lipid molecules that compose membranes have a wide variety of structures. The biosynthesis of common lipid molecules is well understood, conversely the biosynthesis of more exotic lipid molecules is not. The genomes selected today are from the archaeal phylum Thaumarchaeota, which are known for being anammox archaea, meaning they convert ammonia to nitrite. In addition to that, these archaea all produce a series of unusual lipids. Since these lipids are shared across the entire phylum from which these genomes are pulled, it is possible that the biosynthetic genes responsible for lipid construction are shared across the selected genomes, and part of their 'core' genome.

Our goal with this workshop is to show how functions within DECIPHER, in addition to those found within FindHomology can be used to predict and compare the core and pan genome for these selected genomes.

To begin with, we will load the required packages:
DECIPHER and FindHomology will provide tools to complete the phylogenetic analysis, while phytools will provide visualization tools, specifically for the visualization of tangelograms, and stringr is used for a specific string operation.


```r
library(DECIPHER)
```

```
## Loading required package: Biostrings
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colMeans, colnames, colSums, dirname, do.call, duplicated,
##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which, which.max,
##     which.min
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: RSQLite
```

```r
library(devtools)
install_github(repo = "npcooley/FindHomology")
```

```
## Downloading GitHub repo npcooley/FindHomology@master
## from URL https://api.github.com/repos/npcooley/FindHomology/zipball/master
```

```
## Installing FindHomology
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL  \
##   '/tmp/RtmpAXyJvj/devtools64f713e3ee95/npcooley-FindHomology-e9629cf'  \
##   --library='/usr/local/lib/R/R-3.5-bioc-devel' --install-tests
```

```
## 
```

```r
library(FindHomology)
library(stringr)
```

The first step in this process is selecting the genomes we want to work with. A command like the one below, using the NCBI command line utilities, will select FTP directories for the latest, complete genomes for any organism matching the term "Nitrosopumilus" or "Nitrososphaera"from genbank. NCBI command line tools are complex, flexible, and powerful. If you have them installed on your personal system, the esearch command would look like this:


```bash
esearch -db assembly -query 'Nitrosopumilus[organism] OR Nitrososphaera and "complete genome"[filter] AND "latest genbank"[filter]"' | efetch -format docsum | xtract -pattern DocumentSummary -block FtpPath -match "@type:genbank" -element FtpPath | sed 's/$/\//' > ftpdirpaths
```

We can similarly pass this command through R to the terminal if we want to skip this step, or we have a big workflow in R and are unwilling to switch between applications.


```r
system(paste("esearch -db assembly",
             "-query 'Nitrosopumilus[organism] OR Nitrososphaera[Organism] AND \"complete genome\"[filter]",
             "AND \"latest genbank\"[filter]\"\'",
             "|",
             "efetch -format docsum",
             "|",
             "xtract -pattern DocumentSummary -block FtpPath -match",
             "\"@type:genbank\"",
             "-element FtpPath",
             "|",
             "sed ", "'", "s", '/', '$', '/', '\\', '//', "'",
             ">",
             "ftpdirpaths",
             sep = " "))
```

There is currently a bug in using NCBI commands directly through RStudio. If you chose to commandeer this code for your own use, the paste statement, wrapped inside a system call, will work through the regular R console without a problem. The modest workaround below will work in RStudio.[^2]


```r
myTerm <- rstudioapi::terminalCreate(show = FALSE)
rstudioapi::terminalSend(myTerm, paste("esearch -db assembly",
                                       "-query 'Nitrosopumilus[organism] OR Nitrososphaera[Organism] AND \"complete genome\"[filter]",
                                       "AND \"latest genbank\"[filter]\"\'",
                                       "|",
                                       "efetch -format docsum",
                                       "|",
                                       "xtract -pattern DocumentSummary -block FtpPath -match",
                                       "\"@type:genbank\"",
                                       "-element FtpPath",
                                       "|",
                                       "sed ", "'", "s", '/', '$', '/', '\\', '//', "'",
                                       ">",
                                       "ftpdirpaths",
                                       sep = " "))
Sys.sleep(1)
repeat{
    Sys.sleep(0.1)
    if(rstudioapi::terminalBusy(myTerm) == FALSE){
        print("Code Executed")
        break
    }
}
```

Once these file paths are collected, however you choose to do so, we can grab genomic FASTA files for the genomes, and the .gff files for gene calls and annotations. There are several ways to collect gene calls, and annotations, such as Prodigal (for gene calls) and Prokka (for annotations), which by default uses Prodigal to find gene boundaries. Using the NCBI annotations is simple and conducive for this workshop. However if you are working on sequence data that you have generated or collected yourself, NCBI annotations will not be available. 





```r
Genomes <- ftpdirpaths
x <- strsplit(Genomes,
              split = "/")
Annotations <- vector("character",
                      length = length(Genomes))
GeneCalls <- vector("list",
                    length = length(Genomes))
for (i in seq_along(Annotations)) {
  Annotations[i] <- paste(Genomes[i],
                      x[[i]][10L],
                      "_genomic.gff.gz",
                      sep = "")
  z1 <- gzcon(url(Annotations[i]))
  z2 <- textConnection(readLines(z1))
  z3 <- readLines(z2)
  z4 <- strsplit(z3,
                 split = "\t")
  Start <- sapply(z4,
                  function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                     yes = x[4],
                                     no = NA),
                  USE.NAMES = FALSE)
  Stop <- sapply(z4,
                 function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                    yes = x[5],
                                    no = NA),
                 USE.NAMES = FALSE)
  Strand <- sapply(z4,
                   function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                      yes = x[7],
                                      no = NA),
                   USE.NAMES = FALSE)
  Product <- sapply(z4,
                    function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                       yes = x[9],
                                       no = NA),
                    USE.NAMES = FALSE)
  Product <- str_extract(Product,
                         "(?<=product=)(.*)(?=;protein_id)")
  # Subsetting here for positions which have starts, stops, strands, and annotations
  z5 <- as.integer(Start[which(!is.na(Start) &
                                 !is.na(Stop) &
                                 !is.na(Strand) &
                                 !is.na(Product))])
  z6 <- as.integer(Stop[which(!is.na(Start) &
                                !is.na(Stop) &
                                !is.na(Strand) &
                                !is.na(Product))])
  z7 <- Strand[which(!is.na(Start) &
                       !is.na(Stop) &
                       !is.na(Strand) &
                       !is.na(Product))]
  z8 <- Product[which(!is.na(Start) &
                        !is.na(Stop) &
                        !is.na(Strand) &
                        !is.na(Product))]
  GeneCalls[[i]] <- data.frame("Start" = z5,
                               "Stop" = z6,
                               "Strand" = z7,
                               "Annotation" = z8,
                               stringsAsFactors = FALSE)
}
for (i in seq_along(Genomes)) {
  Genomes[i] <- paste(Genomes[i],
                      x[[i]][10L],
                      "_genomic.fna.gz",
                      sep = "")
}
```

Viable genome web addresses can be directly read by BioStrings to provide a DNAStringSet object. A fundamental object for working with sequence data in R.


```r
readDNAStringSet(Genomes[1])
```

```
##   A DNAStringSet instance of length 1
##       width seq                                        names               
## [1] 2833868 TTGTTAAGGTAAAATTAACA...AAACGACCTGTTAACATGA CP002408.1 Candid...
```

## DECIPHER and Databases

DECIPHER works with sequence data by importing it directly into SQL databases. These databases provide a concise method for organization and storage of data. From these databases sequences can be accessed quickly and easily, while accidental changes to sequences during workflow are prevented. DECIPHER will connect to databases through either a Database connection, or a character vector of a filepath to a sqlite file. In this first chunk of code we will use a database connection.


```r
DBPath <- tempfile() # in your own work a temp file is likely inadvisable, for this workshop it is used for simplicity

DBConn <- dbConnect(SQLite(),
                    DBPath)

# Genomes is a file of character vectors of ftp addresses
for (i in seq_along(Genomes)) {
  Seqs2DB(seqs = Genomes[i],
          type = "FASTA",
          dbFile = DBConn,
          identifier = as.character(i),
          tblName = "Seqs")
}
```

```
## Reading FASTA file chunk 1
## 
## 1 total sequence in table Seqs.
## Time difference of 0.61 secs
## 
## Reading FASTA file chunk 1
## 
## Added 1 new sequence to table Seqs.
## 2 total sequences in table Seqs.
## Time difference of 0.83 secs
## 
## Reading FASTA file chunk 1
## 
## Added 1 new sequence to table Seqs.
## 3 total sequences in table Seqs.
## Time difference of 0.6 secs
## 
## Reading FASTA file chunk 1
## 
## Added 1 new sequence to table Seqs.
## 4 total sequences in table Seqs.
## Time difference of 0.67 secs
## 
## Reading FASTA file chunk 1
## 
## Added 1 new sequence to table Seqs.
## 5 total sequences in table Seqs.
## Time difference of 0.69 secs
## 
## Reading FASTA file chunk 1
## 
## Added 1 new sequence to table Seqs.
## 6 total sequences in table Seqs.
## Time difference of 0.53 secs
## 
## Reading FASTA file chunk 1
## 
## Added 1 new sequence to table Seqs.
## 7 total sequences in table Seqs.
## Time difference of 0.67 secs
## 
## Reading FASTA file chunk 1
## 
## Added 1 new sequence to table Seqs.
## 8 total sequences in table Seqs.
## Time difference of 0.65 secs
```

```r
dbDisconnect(DBConn)
# identifying details of the database can be viewed in your systems browser:
# BrowseDB(DBPath)
```

## Comparison of genomes

Comparison of genomes will be accomplished with Synteny. Generically when regions in a genome are the same, they are syntenic matches. Often this is accomplished through matching runs of reciprocal best blast hits. Generically, if in genome 1, genes A-B-C are reciprocal best blast hits to genes Z-Y-X in genome 2, the region encompassing those genes would be considered syntenic.

In DECIPHER this is accomplished using k-mer matching, and by default works in translated space. So in this workflow, we will be determining Synteny between genomes before we consider the gene calls we have already collected.

**Simple Examples of Synteny**


```r
# Here we will use a file path instead of a dbconnection
SyntenyObject <- FindSynteny(dbFile = DBPath,
                             verbose = TRUE)
```

```
## ===========================================================================
## 
## Time difference of 57.2 secs
```

```r
plot(SyntenyObject)
```

<img src="DECIPHERWorkshop_files/figure-html/unnamed-chunk-1-1.png" width="672" />

```r
plot(SyntenyObject,
     "frequency")
```

<img src="DECIPHERWorkshop_files/figure-html/unnamed-chunk-2-1.png" width="672" />

```r
plot(SyntenyObject,
     "neighbor")
```

<img src="DECIPHERWorkshop_files/figure-html/unnamed-chunk-3-1.png" width="672" />

```r
pairs(SyntenyObject[4:5, 4:5])
```

<img src="DECIPHERWorkshop_files/figure-html/unnamed-chunk-4-1.png" width="672" />

```r
pairs(SyntenyObject)
```

<img src="DECIPHERWorkshop_files/figure-html/unnamed-chunk-5-1.png" width="672" />

## Predict Homology between genes that are linked by syntenic hits

The FindHomology package contains 3 functions, the first two of which run below. NucleotideOverlap takes in an object of class Synteny, and a list of GeneCalls. We previously generated the Synteny object  "SyntenyObject", using DECIPHER's function FindSynteny, and a list of Gene calls and annotations, called "GeneCalls". Here we will use these two objects to build what is functionally a list of lists of pairs of genes that are linked by hits.

This is a fun object to scroll through, but is not represented well visually, a simple example would be:

| Genome 1 | Genome 2 | # of Nucleotides |
|----------|----------|------------------|
| Gene A   | Gene Z   |             25L  |
| Gene B   | Gene Y   |             14L  |
| Gene C   | Gene X   |            145L  |
| Gene ... | Gene ... |             ...  |

The reason we name this "MatrixObject" below is because the shape and dimensions of our data are dependent upon our Synteny object. As in, if our Synteny object we built with 5 genomes, that object is a 5 x 5 matrix. NucleotideOverlap accesses the upper triangle of that object to build a 5 x 5 matrix where each position is built from data in the analogous position from the Synteny object. MatrixObject[1, 2] below was created from SyntenyObject[1, 2] and so on and so forth.


```r
MatrixObject <- NucleotideOverlap(SyntenyObject = SyntenyObject,
                                  GeneCalls = GeneCalls,
                                  Verbose = TRUE)
```

```
## ===========================================================================
## Time difference of 1.000226 mins
```




```r
Homologs <- Catalog(MatrixObject,
                    Verbose = TRUE)
```

```
## ===========================================================================
## Time difference of 1.901819 mins
```



```r
hist(sapply(Homologs,
            function(x) nrow(x)))
```

<img src="DECIPHERWorkshop_files/figure-html/find homology objects part 1-1.png" width="672" />

```r
# In this case, this is trivial, but in cases with larger sets of genomes, it may not be
MaxRows <- max(sapply(Homologs,
                      function(x) nrow(x)),
               na.rm = TRUE)
# Genes that are homologous across the whole set of genomes, and are completely transitive, represent the 'core' genome for this set of genomes
CoreSet <- which(sapply(Homologs,
                        function(x) nrow(x)) == MaxRows)
```
Use the core genome to create an alignment

```r
# Select core genes individually from each genome and align them, then concatonate the aligments
CoreGenome <- CoreAligner(Homologs[CoreSet],
                          PATH = DBPath,
                          GeneCalls = GeneCalls,
                          Verbose = TRUE)
```

```
## ===========================================================================
## Time difference of 32.80203 secs
```

```r
# Create a distance matrix, that can be used to create a dendrogram
CoreDist <- DistanceMatrix(myXStringSet = CoreGenome,
                           verbose = FALSE,
                           correction = "Jukes-Cantor")

CoreDend <- IdClusters(myDistMatrix = CoreDist,
                       myXStringSet = CoreGenome,
                       method = "NJ",
                       verbose = FALSE,
                       showPlot = TRUE,
                       type = "dendrogram")
```

<img src="DECIPHERWorkshop_files/figure-html/align core genes-1.png" width="672" />




The presence and absence of genes across the set of genomes can also be used to create a distance matrix

```r
# Create a logical presence/absence matrix for orthologs
LogicalMatrix <- vector("list",
                        length = length(Homologs))
pBar <- txtProgressBar(style = 1L)
for (i in seq_along(Homologs)) {
  LogicalMatrix[[i]] <- vector("logical",
                               length = ncol(Homologs[[i]]))
  for (j in seq_along(LogicalMatrix[[i]])) {
    if (length(unique(Homologs[[i]][, j])) > 1L) {
      LogicalMatrix[[i]][j] <- TRUE
    }
  }
  setTxtProgressBar(pb = pBar,
                    value = i/length(LogicalMatrix))
}
```

```
## ===========================================================================
```

```r
LogicalMatrix <- do.call(cbind, LogicalMatrix)

# Use this logical matrix to create a distance matrix
PanGenome <- dist(LogicalMatrix,
                  method = "binary")

# Create a dendrogram
PanDend <- IdClusters(myDistMatrix = PanGenome,
                      method = "NJ",
                      type = "dendrogram",
                      showPlot = TRUE,
                      verbose = TRUE)
```

```
## ===========================================================================
```

<img src="DECIPHERWorkshop_files/figure-html/build pan genome-1.png" width="672" />

```
## 
## 
## Time difference of 0.02 secs
```

Create a simple tangleogram by comparing the core and pan dendrograms


```r
tf1 <- tempfile()
tf2 <- tempfile()

WriteDendrogram(x =PanDend,
                file = tf1)

WriteDendrogram(x = CoreDend,
                file = tf2)
unlink(tf1)
unlink(tf2)
layout(matrix(1:2, nrow = 1L))
p <- par(mar = c(4, 5, 2, 0))
plot(CoreDend, horiz = TRUE, leaflab = "none")
par(mar = c(4, 4, 2, 1))
plot(PanDend, horiz = TRUE, xlim = c(0, attr(PanDend, "height")), leaflab = "none")
C.Ord <- unlist(CoreDend)
P.Ord <- unlist(PanDend)
segments(-0.05, seq_along(C.Ord), -0.01, match(C.Ord, P.Ord), xpd = TRUE)
```

<img src="DECIPHERWorkshop_files/figure-html/simple tangleogram-1.png" width="672" />

```r
par(p)
```


TODO: More Text, several simpler examples to precede synteny plots, references

[^2]: Thanks to [Serdar Balci](https://twitter.com/serdarbalci) for providing this workaround!







