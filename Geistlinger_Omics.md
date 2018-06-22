
# Functional enrichment analysis of high-throughput omics data

## Instructor(s) name(s) and contact information

* Ludwig Geistlinger (<Ludwig.Geistlinger@sph.cuny.edu>)
* Levi Waldron

CUNY School of Public Health
55 W 125th St, New York, NY 10027



## Workshop Description

This workshop gives an in-depth overview of existing methods for enrichment 
analysis of gene expression data with regard to functional gene sets, pathways, 
and networks.
The workshop will help participants understand the distinctions between 
assumptions and hypotheses of existing methods as well as the differences in 
objectives and interpretation of results. 
It will provide code and hands-on practice of all necessary steps for differential
 expression analysis, gene set- and network-based enrichment analysis, and 
identification of enriched genomic regions and regulatory elements, along with 
visualization and exploration of results. 
A previous well-attended (>50 participants) version of this workshop was presented at Bioc2017 
(https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/OmicsData/doc/enrichOmics.html).

### Pre-requisites

* Basic knowledge of R syntax
* Familiarity with the SummarizedExperiment class
* Familiarity with the GenomicRanges class

* Familiarity with high-throughput gene expression data as obtained with 
  microarrays and RNA-seq
* Familiarity with the concept of differential expression analysis 
  (with e.g. limma, edgeR, DESeq2)

### Workshop Participation

Execution of example code and hands-on practice

### _R_ / _Bioconductor_ packages used

* EnrichmentBrowser
* regioneR

### Time outline

| Activity                              | Time |
|---------------------------------------|------|
| Background                            | 30m  |
| Differential expression analysis      | 15m  |
| Gene set analysis                     | 30m  |
| Gene network analysis                 | 15m  |
| Genomic region analysis               | 15m  |


## Goals and objectives

Theory
* Gene sets, pathways & regulatory networks
* Resources
* Gene set analysis vs. gene set enrichment analysis
* Underlying null: competitive vs. self-contained
* Generations: ora, fcs & topology-based

Practice:
* Data types: microarray vs. RNA-seq
* Differential expression analysis
* Defining gene sets according to GO and KEGG
* GO/KEGG overrepresentation analysis
* Functional class scoring & permutation testing
* Network-based enrichment analysis
* Genomic region enrichment analysis
