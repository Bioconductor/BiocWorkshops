
---
knit: "bookdown::render_book"
title: "The Bioconductor 2018 Workshop Compilation"
description: "This book contains all the workshops presented at the Bioconductor 2018 Conference"
site: bookdown::bookdown_site
github-repo: Bioconductor/BiocWorkshops
documentclass: book
bibliography: [
  "201_Love_DESeq2.bib",
  "202_Das_SingleCellRNASeq.bib",
  "230_cyto_automation.bib",
  "Lee_Plyranges/plyranges.bib",
  "Waldron_PublicData/Waldron_PublicData.bib"
  ]
---

# Introduction

Editors:
    Levi Waldron^[City University of New York, New York, NY]
    Sean Davis^[Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD]
    Marcel Ramos^[City University of New York and Roswell Park Comprehensive Cancer Center, NY]
    Lori Shepherd^[Roswell Park Comprehensive Cancer Center, Buffalo, NY].
    Martin Morgan^[Roswell Park Comprehensive Cancer Center, Buffalo, NY].
    <br/>
Last modified: 30 July, 2018.

Welcome to Bioc2018. This year's conference includes a wide array of workshops
for audiences ranging from beginner to advance users and developers. Workshop
materials are available as a book in html, pdf, and eBook format at
https://bioconductor.github.io/BiocWorkshops/. Workshops are organized by level
and topic according to numbers, as described below. Every workshop starts with a
syllabus that will help you to decide whether it matches your learning goals.

## The Workshops

This book contains workshops for _R_ / _Bioconductor_
training. The workshops are divided into 3 sections:

- **Learn** (100-series chapters) contains material for beginning
  users of _R_ and _Bioconductor_. However, even experienced _R_ and _Bioconductor_ 
  users may find something new here. 
    - 100: _R_ and _Bioconductor_ for everyone: an introduction
    - 101: Introduction to Bioconductor annotation resources
    - 102: Solving common bioinformatic challenges using GenomicRanges
    - 103: Public data resources and Bioconductor

- **Use** (200-series chapters) contains workshops emphasizing use of
  _Bioconductor_ for common tasks, e.g., RNA-seq differential
  expression, single-cell analysis, gene set enrichment, multi'omics analysis, 
  genome analysis, network analysis, and pharmacogenomics.
    - 200: RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR
    - 201: RNA-seq data analysis with DESeq2
    - 202: Analysis of single-cell RNA-seq data: Dimensionality reduction, clustering, and lineage inference
    - 210: Functional enrichment analysis of high-throughput omics data
    - 220: Workflow for multi-omics analysis with MultiAssayExperiment
    - 230: Cytoscape automation in R using Rcy3
    - 240: Fluent genomic data analysis with plyranges
    - 250: Working with genomic data in R with the DECIPHER package
    - 260: Biomarker discovery from large pharmacogenomics datasets
    
- **Develop** (500-series chapters) contains workshops to help expert
  users hone their skills and contribute their domain-specific
  knowledge to the _Bioconductor_ community. These workshops are presented
  on "Developer Day". 
    - 500: Effectively using the DelayedArray framework to support the analysis of large datasets
    - 510: Maintaining your Bioconductor package

## How to use these workshops

These workshops have a lot of package dependencies, and some use data
that you must have on disk. There are several ways to run the code
from these workshops yourself.

### Install on your own computer

These workshops were developed for Bioconductor 3.8 (development
branch) to allow teaching the most up-to-date methods. Some, but not
all, workshop materials will work on Bioconductor 3.7, and (almost?)
all should work after the release of Bioconductor 3.8 in October
2018. To run the workshops on your own computer, you should install
Bioconductor >= 3.8 (which is the development version before October
2018, and the release version thereafter). The following commands
should then install all needed dependencies, at least for Linux:

```
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("Bioconductor/BiocWorkshops")
```

We have noticed that on Windows and Mac, this may not install the
required annotation and experimental data packages. If you get errors
about missing dependencies, you can install these with additional
calls to the `BiocManager::install()` function.


### Use the exact same AMI as the workshop

Each participant in the Bioc2018 workshops was provided with their own
machine image (called an Amazon Machine Image [AMI]) that contained
up-to-date versions of R, required R packages, all necessary operating
system libraries, and all workshop materials. This image has been
tested by the organizers and more than 100 workshop participants, so
if you use it, *everything will just work.* Using this image is the
most hassle-free way to work through these workshops yourself without
having to worry about setup. 

To use this image yourself (including the AMI image #)... (Lori/Sean)

To offer this image to numerous students in a workshop... (Lori/Sean)

### Make your own Docker image and AMI

If you want to alter the image in some way, you can rebuild it using
[packer](https://packer.io), a toolkit that enables for the AMI
creation process. Clone the GitHub repo:

https://github.com/Bioconductor/BiocWorkshops

Navigate to its *packer* subdirectory, then you should be able to do
`packer build bioc_2018.json` (after packer is installed and AWS keys
sorted out) to create a new AMI based on the Bioc R 3.5.1 devel
AMI. You can see the "provision" section for details of the build
process (Rscripts, basically). Of course, the json file can be altered
for customization, including pegging to specific BiocWorkflow tags, if
interested.


[Bioconductor 2018 Conference]: https://bioc2018.bioconductor.org/
