
---
Authors:
    Levi Waldron^[City University of New York, New York, NY, USA]
    Benjamin Haibe-Kain^[Princess Margaret Cancer Center, Toronto, Canada]
    Sean Davis^[Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD, USA]
bibliography: Waldron_PublicData/Waldron_PublicData.bib
Last modified: 2 July, 2018.
---



# Public Data Resources and Bioconductor

## Workshop Description 

The goal of this workshop is to introduce Bioconductor packages for finding,
accessing, and using large-scale public data resources including the 
Gene Expression Omnibus [GEO](https://www.ncbi.nlm.nih.gov/geo), Sequence
Read Archive [SRA](https://www.ncbi.nlm.nih.gov/sra), the Genomic Data
Commons [GDC](https://portal.gdc.cancer.gov/), and Bioconductor-hosted 
curated data resources for metagenomics, pharmacogenomics, and The Cancer 
Genome Atlas.

### Pre-requisites

* Basic knowledge of R syntax
* Familiarity with the ExpressionSet and SummarizedExperiment classes
* Basic familiarity with 'omics technologies such as microarray and NGS sequencing

Interested students can prepare by reviewing vignettes of the packages listed in "R/Bioconductor packages used" to gain background on aspects of interest to them.

Some more general background on these resources is published in @Kannan2016-yv.

### Workshop Participation 

Each component will include runnable examples of typical usage that students are encouraged to run during demonstration of that component.

### R/Bioconductor packages used

* [GEOquery](http://bioconductor.org/packages/GEOquery/): Access to the NCBI Gene Expression Omnibus (GEO), a public repository of gene expression (primarily microarray) data.
* [GenomicDataCommons](http://bioconductor.org/packages/GenomicDataCommons/): Access to the NIH / NCI Genomic Data Commons RESTful service.
* [SRAdb](http://bioconductor.org/packages/SRAdb/): A compilation of metadata from the NCBI Sequence Read Archive, the largest public repository of sequencing data from the next generation of sequencing platforms, and tools
* [curatedTCGAData](http://bioconductor.org/packages/curatedTCGAData/): Curated data from The Cancer Genome Atlas (TCGA) as MultiAssayExperiment Objects
* [curatedMetagenomicData](http://bioconductor.org/packages/curatedMetagenomicData/): Curated metagenomic data of the human microbiome
* [HMP16SData](http://bioconductor.org/packages/HMP16SData/): Curated metagenomic data of the human microbiome
* [PharmacoGx](https://bioconductor.org/packages/PharmacoGx/): Analysis of large-scale pharmacogenomic data


### Time outline

This is a 1h45m workshop.

| Activity                            | Time    |
|-------------------------------------|---------|
| Overview | 10m |
| GEOquery | 15m |
| GenomicDataCommons | 20m |
| Sequence Read Archive | 20m |
| curatedTCGAData   | 10m |
| curatedMetagenomicData and HMP16SData | 15m |
| PharmacoGx | 20m |

## workshop goals and objectives

Bioconductor provides access to significant amounts of publicly available 
experimental  data. This workshop introduces students to Bioconductor
interfaces to the NCBI's Gene Expression Omnibus, Genomic Data Commons,
Sequence Read Archive and PharmacoDB. It additionally introduces curated resources 
providing The Cancer Genome Atlas, the Human Microbiome Project and other 
microbiome studies, and major pharmacogenomic studies, as native Bioconductor
objects ready for analysis and comparison to in-house datasets.

### Learning goals

* search NCBI resources for publicly available 'omics data
* quickly use data from the TCGA and the Human Microbiome Project

### Learning objectives

* find and download processed microarray and RNA-seq datasets from the Gene Expression Omnibus
* find and download 'omics data from the Genomic Data Commons and Sequence Read Archive
* download and manipulate data from The Cancer Genome Atlas and Human Microbiome Project
* download and explore pharmacogenomics data

## Overview

## GEOquery

## GenomicDataCommons

## Sequence Read Archive

## Accessing The Cancer Genome Atlas (TCGA)

### curatedTCGAData

Load packages:


```r
suppressPackageStartupMessages({
    library(curatedTCGAData)
    library(MultiAssayExperiment)
})
```

Checking available cancer codes and assays in TCGA data:


```r
curatedTCGAData(diseaseCode = "*", assays = "*", dry.run = TRUE)
#> Please see the list below for available cohorts and assays
#> Available Cancer codes:
#>  ACC BLCA BRCA CESC CHOL COAD DLBC ESCA GBM HNSC KICH
#>  KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PCPG
#>  PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM 
#> Available Data Types:
#>  CNACGH CNASeq CNASNP CNVSNP GISTICA GISTICT
#>  Methylation miRNAArray miRNASeqGene mRNAArray
#>  Mutation RNASeq2GeneNorm RNASeqGene RPPAArray
```

Check potential files to be downloaded:


```r
curatedTCGAData(diseaseCode = "COAD", dry.run = TRUE)
#>                                  COAD_CNASeq 
#>                   "COAD_CNASeq-20160128.rda" 
#>                                  COAD_CNASNP 
#>                   "COAD_CNASNP-20160128.rda" 
#>                                  COAD_CNVSNP 
#>                   "COAD_CNVSNP-20160128.rda" 
#>                        COAD_GISTIC_AllByGene 
#>         "COAD_GISTIC_AllByGene-20160128.rda" 
#>                COAD_GISTIC_ThresholdedByGene 
#> "COAD_GISTIC_ThresholdedByGene-20160128.rda" 
#>                            COAD_Methylation1 
#>     "COAD_Methylation_methyl27-20160128.rda" 
#>                            COAD_Methylation2 
#>    "COAD_Methylation_methyl450-20160128.rda" 
#>                            COAD_miRNASeqGene 
#>             "COAD_miRNASeqGene-20160128.rda" 
#>                               COAD_mRNAArray 
#>                "COAD_mRNAArray-20160128.rda" 
#>                                COAD_Mutation 
#>                 "COAD_Mutation-20160128.rda" 
#>                         COAD_RNASeq2GeneNorm 
#>          "COAD_RNASeq2GeneNorm-20160128.rda" 
#>                              COAD_RNASeqGene 
#>               "COAD_RNASeqGene-20160128.rda" 
#>                               COAD_RPPAArray 
#>                "COAD_RPPAArray-20160128.rda"
```

## GBM dataset example


```r
gbm <- curatedTCGAData("COAD", c("RPPA*"), dry.run=FALSE)
gbm
#> A MultiAssayExperiment object of 1 listed
#>  experiment with a user-defined name and respective class. 
#>  Containing an ExperimentList class object of length 1: 
#>  [1] COAD_RPPAArray-20160128: SummarizedExperiment with 208 rows and 360 columns 
#> Features: 
#>  experiments() - obtain the ExperimentList instance 
#>  colData() - the primary/phenotype DataFrame 
#>  sampleMap() - the sample availability DataFrame 
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment 
#>  *Format() - convert into a long or wide DataFrame 
#>  assays() - convert ExperimentList to a SimpleList of matrices
```

**Note**. See the `\BiocPkg("MultiAssayExperiment") vignette for more details on using this object.

### Subtype information

Some cancer datasets contain associated subtype information within the
clinical datasets provided. This subtype information is included in the
metadata of `colData` of the `MultiAssayExperiment` object. To obtain these
variable names, run the `metadata` function on the `colData` of the object
such as:


```r
head(metadata(colData(gbm))[["subtypes"]])
#>        COAD_annotations        COAD_subtype
#> 1            Patient_ID             patient
#> 2                   msi          MSI_status
#> 3  methylation_subtypes methylation_subtype
#> 4         mrna_subtypes  expression_subtype
#> 5 histological_subtypes   histological_type
```



### TCGABiolinks
### GenomicDataCommons

## Microbiome data

Bioconductor provides a couple curated resources of microbiome data. 
Most microbiome data are generated either by targeted amplicon sequencing (usually of variable regions of the 16S ribosomal RNA gene) or by metagenomic shotgun sequencing. 

cheap (multiplex hundreds of samples)
relatively small data
provides genus-level taxonomy and inferred metabolic function for bacteria and archaea


curatedMetagenomicData and HMP16SData

## Pharmacogenomics

