
# Fluent genomic data analysis with plyranges

## Instructor(s) name(s) and contact information

* Stuart Lee (lee.s@wehi.edu.au)
* Michael Lawrence (lawremi@gmail.com)


## Workshop Description

In this workshop, we will give an overview of how to perform low-level
analyses of genomic data using the  grammar of genomic data transformation
defined in the plyranges package. We will cover:

- introduction to GRanges
- overview of the core verbs for arithmetic, restriction, and aggregation of GRanges objects
- performing joins between GRanges objects
- designing pipelines to quickly explore data from AnnotationHub
- reading BAM and other file types as GRanges objects

The workshop will be a computer lab, in which the participants will be able to 
ask questions and interact with the instructors.

### Pre-requisites

This workshop is mostly self-contained however familiarity with the following would be useful:

- plyranges [vignette](https://sa-lee.github.io/plyranges/articles/an-introduction.html)
- the GenomicRanges and IRanges packages
- tidyverse approaches to data analysis

 
### Workshop Participation

Students will work through an Rmarkdown document while the instructors respond to any questions they have.

### _R_ / _Bioconductor_ packages used

- plyranges
- AnnotationHub
- GenomicRanges
- IRanges
- S4Vectors

### Time outline


| Activity                     | Time |
|------------------------------|------|
| Overview of GRanges          | 5m   |
| The plyranges grammar        | 20m  |
| I/O and data pipelines       | 20m  |

## Workshop goals and objectives


### Learning goals


* Understand that GRanges follows tidy data principles
* Apply the plyranges grammar to genomic data analysis


### Learning objectives

* Use AnnotationHub to find and summarise data
* Read files into R as GRanges objects
* Perform coverage analysis 
* Build data pipelines for analysis based on GRanges
