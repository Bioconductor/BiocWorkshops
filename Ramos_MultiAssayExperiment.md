
# Workflow for Multi-omics Analysis with MultiAssayExperiment

## Instructor names and contact information

* Marcel Ramos
* Ludwig Geistlinger
* Levi Waldron

## Workshop Description

This workshop demonstrates data management and analyses of multiple 
assays associated with a single set of biological specimens, 
using the `MultiAssayExperiment` data class and methods. It introduces 
the `RaggedExperiment` data class, which provides efficient and powerful 
operations for representation of copy number and mutation and variant 
data that are represented by different genomic ranges for each specimen. 

### Pre-requisites

List any workshop prerequisites, for example:

* Basic knowledge of R syntax
* Familiarity with the GRanges and SummarizedExperiment classes
* Familiarity with 'omics data types including copy number and gene expression

### workshop Participation

Students will have a chance to build a `MultiAssayExperiment` object 
from scratch, and will also work with more complex objects provided
by the `curatedTCGAData` package.

### R/Bioconductor packages used

* [MultiAssayExperiment]
* [RaggedExperiment]
* [curatedTCGAData]
* [SummarizedExperiment]

### Time outline

1h 45m total

| Activity                            | Time    |
|-------------------------------------|---------|
| Overview of key data classes | 25m |
| Working with RaggedExperiment | 20m |
| Building a MultiAssayExperiment from scratch | 10m |
| Creating and TCGA multi-assay dataset | 10m |
| Subsetting and reshaping multi-assay data | 20m |
| Plotting, correlation, and other analyses | 20m |


## Workshop goals and objectives

### Learning goals

* identify appropriate data structures for different 'omics data types
* gain familiarity with GRangesList and RaggedExperiment

### Learning objectives

* use curatedTCGAData to create custom TCGA MultiAssayExperiment objects
* create a MultiAssayExperiment for TCGA or other multi'omics data
* perform subsetting, reshaping, growing, and extraction of a MultiAssayExperiment
* link MultiAssayExperiment data with packages for differential expression, 
machine learning, and plotting
