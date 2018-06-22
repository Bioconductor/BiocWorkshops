
# Introduction to Bioconductor annotation resources

## Instructors

* James W. MacDonald (jmacdon@uw.edu)
* Lori Shepherd (lori.shepherd@roswellpark.org)

## Workshop Description

There are various annotation packages provided by the Bioconductor
project that can be used to incorporate additional information to
results from high-throughput experiments. This can be as simple as
mapping Ensembl IDs to corresponding HUGO gene symbols, to much more
complex queries involving multiple data sources. In this workshop we
will cover the various classes of annotation packages, what they
contain, and how to use them efficiently. 

### Prerequisites

* Basic knowledge of R syntax
* Basic understanding of the various annotation sources (NCBI, EBI/EMBL)

Useful background reading

* The
  [AnnotationDbi](https://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf)
  vignette.
* The
  [biomaRt](https://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html)
  vignette.
* The
  [GenomicFeatures](https://www.bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf)
  vignette.
  

### Workshop Participation

After each type of annotation package is introduced, students will be
given the opportunity to practice making their own queries. 

### _R_ / _Bioconductor_ packages used

* AnnotationDbi
* AnnotationHub
* BSgenome
* biomaRt
* ensembldb
* org.Hs.eg.db
* TxDb.Hsapiens.UCSC.hg19.knownGene
* EnsDb.Hsapiens.v79
* EnsDb.Mmusculus.v79
* Homo.sapiens 
* BSgenome.Hsapiens.UCSC.hg19
* hugene20sttranscriptcluster.db


## Workshop goals and objectives

Annotating data is a complex task. For any high-throughput experiment
the analyst usually starts with a set of identifiers for each thing
that was measured, and in order to make the results useful to
collaborators these identifiers need to be mapped to other identifiers
that are either more familiar to collaborators, or that can be used
for further analyses. As an example, RNA-Seq data may only have Entrez
Gene IDs for each gene measured, and as part of the output you may
want to include the gene symbols, which are more likely to be familiar
to a Biologist.

### Learning goals

* Understand what sort of annotation data are available
* Understand the difference between annotation sources (NCBI and EBI/EMBL)
* Gain familiarity with the various ways to query annotation packages
* Get some practice making queries

### Learning objectives

* Be able to use select and mapIds to map between identifiers
* Be able to extract data from TxDb and EnsDb packages
* Be able to make queries using biomaRt
* Extract and utilize various data from AnnotationHub
