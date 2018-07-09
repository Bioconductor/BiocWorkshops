

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

* [MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment)
* [RaggedExperiment](https://bioconductor.org/packages/RaggedExperiment)
* [curatedTCGAData](https://bioconductor.org/packages/curatedTCGAData)
* [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)
* [TCGAutils](https://bioconductor.org/packages/TCGAutils)
* [UpSetR](https://bioconductor.org/packages/UpSetR)
* [EnsDb.Hsapiens.v86](https://bioconductor.org/packages/EnsDb.Hsapiens.v86)
* [survminer](https://cran.r-project.org/package=survminer)
* [pheatmap](https://cran.r-project.org/package=pheatmap)

### Time outline

1h 45m total

| Activity                            | Time    |
|-------------------------------------|---------|
| Overview of key data classes | 25m |
| Working with RaggedExperiment | 20m |
| Building a MultiAssayExperiment from scratch | 10m |
| TCGA multi-assay dataset | 10m |
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

## Overview of key data classes

### `SummarizedExperiment` / `RaggedExperiment`

![SummarizedExperiment_Schematic](Ramos_MultiAssayExperiment/SummarizedExperiment.svg)

_The_ standard Bioconductor class convenient for storing experimental results
produced by sequencing and microarry experiments. Each object can store
multiple experiments with associated metadata such as features and samples, as
well as phenotype/clinical data. Supersedes the use of `ExpressionSet`.

### `MultiAssayExperiment`

An integrative container for coordinating multi-omics experiment data on a
set of biological specimens.

![MultiAssayExperiment_Schematic](Ramos_MultiAssayExperiment/MultiAssayExperiment.png)

### `RaggedExperiment`

A flexible data representation for copy number, mutation, and other
ragged array schema for genomic location data.

![RaggedExperiment_Schematic](Ramos_MultiAssayExperiment/RaggedExperiment.png)

## Working with RaggedExperiment

### Constructing a `RaggedExperiment` object

We start with a couple of `GRanges` objects, each representing an individual
sample:


```r
library(GenomicRanges)
sample1 <- GRanges(
    c(A = "chr1:1-10:-", B = "chr1:8-14:+", C = "chr2:15-18:+"),
    score = 3:5)
sample2 <- GRanges(
    c(D = "chr1:1-10:-", E = "chr2:11-18:+"),
    score = 1:2)
```

Include column data `colData` to describe the samples:


```r
colDat <- DataFrame(id = 1:2)
```

### Using `GRanges` objects


```r
library(RaggedExperiment)

ragexp <- RaggedExperiment(
    sample1 = sample1,
    sample2 = sample2,
    colData = colDat)

ragexp
#> class: RaggedExperiment 
#> dim: 5 2 
#> assays(1): score
#> rownames(5): A B C D E
#> colnames(2): sample1 sample2
#> colData names(1): id
```

It can also be constructed from a `GRangesList`:


```r
library(GenomicRanges)
example(GRangesList)
#> 
#> GRngsL> ## Construction with GRangesList():
#> GRngsL> gr1 <-
#> GRngsL+   GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
#> GRngsL+           strand = "+", score = 5L, GC = 0.45)
#> 
#> GRngsL> gr2 <-
#> GRngsL+   GRanges(seqnames = c("chr1", "chr1"),
#> GRngsL+           ranges = IRanges(c(7,13), width = 3),
#> GRngsL+           strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
#> 
#> GRngsL> gr3 <-
#> GRngsL+   GRanges(seqnames = c("chr1", "chr2"),
#> GRngsL+           ranges = IRanges(c(1, 4), c(3, 9)),
#> GRngsL+           strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
#> 
#> GRngsL> grl <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
#> 
#> GRngsL> grl
#> GRangesList object of length 3:
#> $gr1 
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |     score        GC
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric>
#>   [1]     chr2       3-6      + |         5      0.45
#> 
#> $gr2 
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames ranges strand | score  GC
#>   [1]     chr1    7-9      + |     3 0.3
#>   [2]     chr1  13-15      - |     4 0.5
#> 
#> $gr3 
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames ranges strand | score  GC
#>   [1]     chr1    1-3      - |     6 0.4
#>   [2]     chr2    4-9      - |     2 0.1
#> 
#> -------
#> seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> GRngsL> ## Summarizing elements:
#> GRngsL> elementNROWS(grl)
#> gr1 gr2 gr3 
#>   1   2   2 
#> 
#> GRngsL> table(seqnames(grl))
#>     chr2 chr1
#> gr1    1    0
#> gr2    0    2
#> gr3    1    1
#> 
#> GRngsL> ## Extracting subsets:
#> GRngsL> grl[seqnames(grl) == "chr1", ]
#> GRangesList object of length 3:
#> $gr1 
#> GRanges object with 0 ranges and 2 metadata columns:
#>    seqnames    ranges strand |     score        GC
#>       <Rle> <IRanges>  <Rle> | <integer> <numeric>
#> 
#> $gr2 
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames ranges strand | score  GC
#>   [1]     chr1    7-9      + |     3 0.3
#>   [2]     chr1  13-15      - |     4 0.5
#> 
#> $gr3 
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames ranges strand | score  GC
#>   [1]     chr1    1-3      - |     6 0.4
#> 
#> -------
#> seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> GRngsL> grl[seqnames(grl) == "chr1" & strand(grl) == "+", ]
#> GRangesList object of length 3:
#> $gr1 
#> GRanges object with 0 ranges and 2 metadata columns:
#>    seqnames    ranges strand |     score        GC
#>       <Rle> <IRanges>  <Rle> | <integer> <numeric>
#> 
#> $gr2 
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames ranges strand | score  GC
#>   [1]     chr1    7-9      + |     3 0.3
#> 
#> $gr3 
#> GRanges object with 0 ranges and 2 metadata columns:
#>      seqnames ranges strand | score GC
#> 
#> -------
#> seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> GRngsL> ## Renaming the underlying sequences:
#> GRngsL> seqlevels(grl)
#> [1] "chr2" "chr1"
#> 
#> GRngsL> seqlevels(grl) <- sub("chr", "Chrom", seqlevels(grl))
#> 
#> GRngsL> grl
#> GRangesList object of length 3:
#> $gr1 
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |     score        GC
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric>
#>   [1]   Chrom2       3-6      + |         5      0.45
#> 
#> $gr2 
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames ranges strand | score  GC
#>   [1]   Chrom1    7-9      + |     3 0.3
#>   [2]   Chrom1  13-15      - |     4 0.5
#> 
#> $gr3 
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames ranges strand | score  GC
#>   [1]   Chrom1    1-3      - |     6 0.4
#>   [2]   Chrom2    4-9      - |     2 0.1
#> 
#> -------
#> seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> GRngsL> ## Coerce to IRangesList (seqnames and strand information is lost):
#> GRngsL> as(grl, "IRangesList")
#> IRangesList of length 3
#> $gr1
#> IRanges object with 1 range and 2 metadata columns:
#>           start       end     width |     score        GC
#>       <integer> <integer> <integer> | <integer> <numeric>
#>   [1]         3         6         4 |         5      0.45
#> 
#> $gr2
#> IRanges object with 2 ranges and 2 metadata columns:
#>           start       end     width |     score        GC
#>       <integer> <integer> <integer> | <integer> <numeric>
#>   [1]         7         9         3 |         3       0.3
#>   [2]        13        15         3 |         4       0.5
#> 
#> $gr3
#> IRanges object with 2 ranges and 2 metadata columns:
#>           start       end     width |     score        GC
#>       <integer> <integer> <integer> | <integer> <numeric>
#>   [1]         1         3         3 |         6       0.4
#>   [2]         4         9         6 |         2       0.1
#> 
#> 
#> GRngsL> ## isDisjoint():
#> GRngsL> isDisjoint(grl)
#>  gr1  gr2  gr3 
#> TRUE TRUE TRUE 
#> 
#> GRngsL> ## disjoin():
#> GRngsL> disjoin(grl)  # metadata columns and order NOT preserved
#> GRangesList object of length 3:
#> $gr1 
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]   Chrom2       3-6      +
#> 
#> $gr2 
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames ranges strand
#>   [1]   Chrom1    7-9      +
#>   [2]   Chrom1  13-15      -
#> 
#> $gr3 
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames ranges strand
#>   [1]   Chrom2    4-9      -
#>   [2]   Chrom1    1-3      -
#> 
#> -------
#> seqinfo: 2 sequences from an unspecified genome; no seqlengths
#> 
#> GRngsL> ## Construction with makeGRangesListFromFeatureFragments():
#> GRngsL> filepath <- system.file("extdata", "feature_frags.txt",
#> GRngsL+                         package="GenomicRanges")
#> 
#> GRngsL> featfrags <- read.table(filepath, header=TRUE, stringsAsFactors=FALSE)
#> 
#> GRngsL> grl2 <- with(featfrags,
#> GRngsL+              makeGRangesListFromFeatureFragments(seqnames=targetName,
#> GRngsL+                                                  fragmentStarts=targetStart,
#> GRngsL+                                                  fragmentWidths=blockSizes,
#> GRngsL+                                                  strand=strand))
#> 
#> GRngsL> names(grl2) <- featfrags$RefSeqID
#> 
#> GRngsL> grl2
#> GRangesList object of length 3:
#> $XM_001065892.1 
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames              ranges strand
#>          <Rle>           <IRanges>  <Rle>
#>   [1]     chr4 124513961-124514087      +
#>   [2]     chr4 124514706-124515691      +
#> 
#> $XM_578205.2 
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames              ranges strand
#>   [1]     chr2 155875533-155876067      -
#>   [2]     chr2 155879894-155880030      -
#>   [3]     chr2 155895543-155895690      -
#> 
#> $NM_012543.2 
#> GRanges object with 4 ranges and 0 metadata columns:
#>       seqnames            ranges strand
#>   [1]     chr1 96173572-96174077      +
#>   [2]     chr1 96174920-96175330      +
#>   [3]     chr1 96176574-96176785      +
#>   [4]     chr1 96177991-96178484      +
#> 
#> -------
#> seqinfo: 3 sequences from an unspecified genome; no seqlengths
rgx <- RaggedExperiment(grl)
```

### *Assay functions

#### sparseAssay

The most straightforward matrix representation of a `RaggedExperiment` will
return a matrix of dimensions equal to the product of the number of ranges and
samples.


```r
dim(ragexp)
#> [1] 5 2
Reduce(`*`, dim(ragexp))
#> [1] 10
sparseAssay(ragexp)
#>   sample1 sample2
#> A       3      NA
#> B       4      NA
#> C       5      NA
#> D      NA       1
#> E      NA       2
length(sparseAssay(ragexp))
#> [1] 10
```

#### compactAssay

Samples with identical ranges are placed in the same row. Non-disjoint ranges
are **not** collapsed.


```r
compactAssay(ragexp)
#>              sample1 sample2
#> chr1:8-14:+        4      NA
#> chr1:1-10:-        3       1
#> chr2:11-18:+      NA       2
#> chr2:15-18:+       5      NA
```

#### disjoinAssay

This function returns a matrix of disjoint ranges across all samples. Elements
of the matrix are summarized by applying the `simplifyDisjoin` functional
argument to assay values of overlapping ranges.


```r
disjoinAssay(ragexp, simplifyDisjoin = mean)
#>              sample1 sample2
#> chr1:8-14:+        4      NA
#> chr1:1-10:-        3       1
#> chr2:11-14:+      NA       2
#> chr2:15-18:+       5       2
```

## qreduceAssay

The `qreduceAssay` function works with a `query` parameter that highlights
a window of ranges for the resulting matrix. The returned matrix will have
dimensions `length(query)` by `ncol(x)`. Elements contain assay values for the
_i_ th query range and the _j_ th sample, summarized according to the
`simplifyReduce` functional argument.

For demonstration purposes, we can have a look at the original `GRangesList`
and the associated scores from which the current `ragexp` object is derived:


```r
unlist(grl, use.names = FALSE)
#> GRanges object with 5 ranges and 2 metadata columns:
#>       seqnames    ranges strand |     score        GC
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric>
#>   [1]   Chrom2       3-6      + |         5      0.45
#>   [2]   Chrom1       7-9      + |         3       0.3
#>   [3]   Chrom1     13-15      - |         4       0.5
#>   [4]   Chrom1       1-3      - |         6       0.4
#>   [5]   Chrom2       4-9      - |         2       0.1
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

This data is represented as `rowRanges` and `assays` in `RaggedExperiment`:


```r
rowRanges(ragexp)
#> GRanges object with 5 ranges and 0 metadata columns:
#>     seqnames    ranges strand
#>        <Rle> <IRanges>  <Rle>
#>   A     chr1      1-10      -
#>   B     chr1      8-14      +
#>   C     chr2     15-18      +
#>   D     chr1      1-10      -
#>   E     chr2     11-18      +
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
assay(ragexp, "score")
#>   sample1 sample2
#> A       3      NA
#> B       4      NA
#> C       5      NA
#> D      NA       1
#> E      NA       2
```

Here we provide the "query" region of interest:


```r
(query <- GRanges(c("chr1:1-14:-", "chr2:11-18:+")))
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      1-14      -
#>   [2]     chr2     11-18      +
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

 The `simplifyReduce` argument in `qreduceAssay` allows the user to summarize
overlapping regions with a custom method for the given "query" region of
interest. We provide one for calculating a weighted average score per
query range, where the weight is proportional to the overlap widths between
overlapping ranges and a query range.

_Note_ that there are three arguments to this function. See the documentation
for additional details.


```r
weightedmean <- function(scores, ranges, qranges)
{
    isects <- pintersect(ranges, qranges)
    sum(scores * width(isects)) / sum(width(isects))
}
```

A call to `qreduceAssay` involves the `RaggedExperiment`, the `GRanges` query
and the `simplifyReduce` functional argument.


```r
qreduceAssay(ragexp, query, simplifyReduce = weightedmean)
#>              sample1 sample2
#> chr1:1-14:-        3       1
#> chr2:11-18:+       5       2
```

**Note**: Coercion is possible to `SummarizedExperiment`.

## Working with MultiAssayExperiment

### The MultiAssayExperiment miniACC demo object

Get started by trying out `MultiAssayExperiment` using a subset of the TCGA
adrenocortical carcinoma (ACC) dataset provided with the package. This dataset
provides five assays on 92 patients, although all five assays were not performed
for every patient:

1. **RNASeq2GeneNorm**: gene mRNA abundance by RNA-seq
2. **gistict**: GISTIC genomic copy number by gene
3. **RPPAArray**: protein abundance by Reverse Phase Protein Array
4. **Mutations**: non-silent somatic mutations by gene
5. **miRNASeqGene**: microRNA abundance by microRNA-seq.


```r
suppressPackageStartupMessages({
    library(MultiAssayExperiment)
    library(S4Vectors)
})
data(miniACC)
miniACC
#> A MultiAssayExperiment object of 5 listed
#>  experiments with user-defined names and respective classes. 
#>  Containing an ExperimentList class object of length 5: 
#>  [1] RNASeq2GeneNorm: SummarizedExperiment with 198 rows and 79 columns 
#>  [2] gistict: SummarizedExperiment with 198 rows and 90 columns 
#>  [3] RPPAArray: SummarizedExperiment with 33 rows and 46 columns 
#>  [4] Mutations: matrix with 97 rows and 90 columns 
#>  [5] miRNASeqGene: SummarizedExperiment with 471 rows and 80 columns 
#> Features: 
#>  experiments() - obtain the ExperimentList instance 
#>  colData() - the primary/phenotype DataFrame 
#>  sampleMap() - the sample availability DataFrame 
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment 
#>  *Format() - convert into a long or wide DataFrame 
#>  assays() - convert ExperimentList to a SimpleList of matrices
```

### colData - information biological units

This slot is a `DataFrame` describing the characteristics of biological units,
for example clinical data for patients. In the prepared datasets from
[The Cancer Genome Atlas][], each row is one patient and each column is a
clinical, pathological, subtype, or other variable. The `$` function provides a
shortcut for accessing or setting `colData` columns.


```r
colData(miniACC)[1:4, 1:4]
#> DataFrame with 4 rows and 4 columns
#>                 patientID years_to_birth vital_status days_to_death
#>               <character>      <integer>    <integer>     <integer>
#> TCGA-OR-A5J1 TCGA-OR-A5J1             58            1          1355
#> TCGA-OR-A5J2 TCGA-OR-A5J2             44            1          1677
#> TCGA-OR-A5J3 TCGA-OR-A5J3             23            0            NA
#> TCGA-OR-A5J4 TCGA-OR-A5J4             23            1           423
table(miniACC$race)
#> 
#>                     asian black or african american 
#>                         2                         1 
#>                     white 
#>                        78
```

*Key points:*

* One row per patient
* Each row maps to zero or more observations in each experiment in the
`ExperimentList`, below.

### ExperimentList - experiment data

A base `list` or `ExperimentList` object containing the experimental datasets
for the set of samples collected. This gets converted into a class
`ExperimentList` during construction.


```r
experiments(miniACC)
#> ExperimentList class object of length 5: 
#>  [1] RNASeq2GeneNorm: SummarizedExperiment with 198 rows and 79 columns 
#>  [2] gistict: SummarizedExperiment with 198 rows and 90 columns 
#>  [3] RPPAArray: SummarizedExperiment with 33 rows and 46 columns 
#>  [4] Mutations: matrix with 97 rows and 90 columns 
#>  [5] miRNASeqGene: SummarizedExperiment with 471 rows and 80 columns
```

*Key points:*

* One matrix-like dataset per list element (although they do not even need to be
matrix-like, see for example the `RaggedExperiment` package)
* One matrix column per assayed specimen. Each matrix column must correspond to
exactly one row of `colData`: in other words, you must know which patient or
cell line the observation came from. However, multiple columns can come from the
same patient, or there can be no data for that patient.
* Matrix rows correspond to variables, e.g. genes or genomic ranges
* `ExperimentList` elements can be genomic range-based (e.g.
`SummarizedExperiment::RangedSummarizedExperiment-class` or
`RaggedExperiment::RaggedExperiment-class`) or ID-based data (e.g.
`SummarizedExperiment::SummarizedExperiment-class`, `Biobase::eSet-class`
`base::matrix-class`, `DelayedArray::DelayedArray-class`, and derived classes)
* Any data class can be included in the `ExperimentList`, as long as it
supports: single-bracket subsetting (`[`), `dimnames`, and `dim`. Most data
classes defined in Bioconductor meet these requirements.

### sampleMap - relationship graph

`sampleMap` is a graph representation of the relationship between biological
units and experimental results. In simple cases where the column names of
`ExperimentList` data matrices match the row names of `colData`, the user won't
need to specify or think about a sample map, it can be created automatically by
the `MultiAssayExperiment` constructor.  `sampleMap` is a simple three-column
`DataFrame`:

1. `assay` column: the name of the assay, and found in the names of
`ExperimentList` list names
2. `primary` column: identifiers of patients or biological units, and found in
the row names of `colData`
3.  `colname` column: identifiers of assay results, and found in the column
names of `ExperimentList` elements
Helper functions are available for creating a map from a list. See `?listToMap`


```r
sampleMap(miniACC)
#> DataFrame with 385 rows and 3 columns
#>               assay      primary                      colname
#>            <factor>  <character>                  <character>
#> 1   RNASeq2GeneNorm TCGA-OR-A5J1 TCGA-OR-A5J1-01A-11R-A29S-07
#> 2   RNASeq2GeneNorm TCGA-OR-A5J2 TCGA-OR-A5J2-01A-11R-A29S-07
#> 3   RNASeq2GeneNorm TCGA-OR-A5J3 TCGA-OR-A5J3-01A-11R-A29S-07
#> 4   RNASeq2GeneNorm TCGA-OR-A5J5 TCGA-OR-A5J5-01A-11R-A29S-07
#> 5   RNASeq2GeneNorm TCGA-OR-A5J6 TCGA-OR-A5J6-01A-31R-A29S-07
#> ...             ...          ...                          ...
#> 381    miRNASeqGene TCGA-PA-A5YG TCGA-PA-A5YG-01A-11R-A29W-13
#> 382    miRNASeqGene TCGA-PK-A5H8 TCGA-PK-A5H8-01A-11R-A29W-13
#> 383    miRNASeqGene TCGA-PK-A5H9 TCGA-PK-A5H9-01A-11R-A29W-13
#> 384    miRNASeqGene TCGA-PK-A5HA TCGA-PK-A5HA-01A-11R-A29W-13
#> 385    miRNASeqGene TCGA-PK-A5HB TCGA-PK-A5HB-01A-11R-A29W-13
```

*Key points:*

* relates experimental observations (`colnames`) to `colData`
* permits experiment-specific sample naming, missing, and replicate observations

<p style="text-align: right;"> [back to top](#overview-of-key-data-classes) </p>

### metadata

Metadata can be used to keep additional information about patients, assays
performed on individuals or on the entire cohort, or features such as genes,
proteins, and genomic ranges. There are many options available for storing
metadata. First, `MultiAssayExperiment` has its own metadata for describing the
entire experiment:


```r
metadata(miniACC)
#> $title
#> [1] "Comprehensive Pan-Genomic Characterization of Adrenocortical Carcinoma"
#> 
#> $PMID
#> [1] "27165744"
#> 
#> $sourceURL
#> [1] "http://s3.amazonaws.com/multiassayexperiments/accMAEO.rds"
#> 
#> $RPPAfeatureDataURL
#> [1] "http://genomeportal.stanford.edu/pan-tcga/show_target_selection_file?filename=Allprotein.txt"
#> 
#> $colDataExtrasURL
#> [1] "http://www.cell.com/cms/attachment/2062093088/2063584534/mmc3.xlsx"
```

Additionally, the `DataFrame` class used by `sampleMap` and `colData`, as well
as the `ExperimentList` class, similarly support metadata. Finally, many
experimental data objects that can be used in the `ExperimentList` support
metadata. These provide flexible options to users and to developers of derived
classes.

## MultiAssayExperiment Subsetting

### Single bracket `[`

In pseudo code below, the subsetting operations work on the rows of the
following indices:
1. _i_ experimental data rows
2. _j_ the primary names or the column names (entered as a `list` or `List`)
3. _k_ assay

```
multiassayexperiment[i = rownames, j = primary or colnames, k = assay]
```

Subsetting operations always return another `MultiAssayExperiment`. For example,
the following will return any rows named "MAPK14" or "IGFBP2", and remove any
assays where no rows match:


```r
miniACC[c("MAPK14", "IGFBP2"), , ]
```

The following will keep only patients of pathological stage iv, and all their
associated assays:

```r
miniACC[, miniACC$pathologic_stage == "stage iv", ]
#> harmonizing input:
#>   removing 311 sampleMap rows with 'colname' not in colnames of experiments
#>   removing 74 colData rownames not in sampleMap 'primary'
```

And the following will keep only the RNA-seq dataset, and only patients for
which this assay is available:

```r
miniACC[, , "RNASeq2GeneNorm"]
#> harmonizing input:
#>   removing 13 colData rownames not in sampleMap 'primary'
```

### Subsetting by genomic ranges

If any ExperimentList objects have features represented by genomic ranges (e.g.
`RangedSummarizedExperiment`, `RaggedExperiment`), then a `GRanges` object in
the first subsetting position will subset these objects as in
`GenomicRanges::findOverlaps()`.


### Double bracket `[[`

The "double bracket" method (`[[`) is a convenience function for extracting
a single element of the `MultiAssayExperiment` `ExperimentList`. It avoids
the use of `experiments(mae)[[1L]]`. For example, both of the following extract
the `ExpressionSet` object containing RNA-seq data:


```r
miniACC[[1L]]  #or equivalently, miniACC[["RNASeq2GeneNorm"]]
#> class: SummarizedExperiment 
#> dim: 198 79 
#> metadata(3): experimentData annotation protocolData
#> assays(1): exprs
#> rownames(198): DIRAS3 MAPK14 ... SQSTM1 KCNJ13
#> rowData names(0):
#> colnames(79): TCGA-OR-A5J1-01A-11R-A29S-07
#>   TCGA-OR-A5J2-01A-11R-A29S-07 ... TCGA-PK-A5HA-01A-11R-A29S-07
#>   TCGA-PK-A5HB-01A-11R-A29S-07
#> colData names(0):
```

## Patients with complete data

`complete.cases()` shows which patients have complete data for all assays:


```r
summary(complete.cases(miniACC))
#>    Mode   FALSE    TRUE 
#> logical      49      43
```

The above logical vector could be used for patient subsetting. More simply,
`intersectColumns()` will select complete cases and rearrange each
`ExperimentList` element so its columns correspond exactly to rows of `colData`
in the same order:


```r
accmatched = intersectColumns(miniACC)
#> harmonizing input:
#>   removing 170 sampleMap rows with 'colname' not in colnames of experiments
#>   removing 49 colData rownames not in sampleMap 'primary'
```

Note, the column names of the assays in `accmatched` are not the same because of
assay-specific identifiers, but they have been automatically re-arranged to
correspond to the same patients. In these TCGA assays, the first three `-`
delimited positions correspond to patient, ie the first patient is
*TCGA-OR-A5J2*:


```r
colnames(accmatched)
#> CharacterList of length 5
#> [["RNASeq2GeneNorm"]] TCGA-OR-A5J2-01A-11R-A29S-07 ...
#> [["gistict"]] TCGA-OR-A5J2-01A-11D-A29H-01 ...
#> [["RPPAArray"]] TCGA-OR-A5J2-01A-21-A39K-20 ...
#> [["Mutations"]] TCGA-OR-A5J2-01A-11D-A29I-10 ...
#> [["miRNASeqGene"]] TCGA-OR-A5J2-01A-11R-A29W-13 ...
```

## Row names that are common across assays

`intersectRows()` keeps only rows that are common to each assay, and aligns them
in identical order. For example, to keep only genes where data are available for
RNA-seq, GISTIC copy number, and somatic mutations:


```r
accmatched2 <- intersectRows(miniACC[, , c("RNASeq2GeneNorm", "gistict", "Mutations")])
rownames(accmatched2)
#> CharacterList of length 3
#> [["RNASeq2GeneNorm"]] DIRAS3 G6PD KDR ERBB3 ... RET CDKN2A MACC1 CHGA
#> [["gistict"]] DIRAS3 G6PD KDR ERBB3 AKT1S1 ... PREX1 RET CDKN2A MACC1 CHGA
#> [["Mutations"]] DIRAS3 G6PD KDR ERBB3 AKT1S1 ... RET CDKN2A MACC1 CHGA
```

<p style="text-align: right;"> [back to top](#overview-of-key-data-classes) </p>

## Extraction

### assay and assays

The `assay` and `assays` methods follow `SummarizedExperiment` convention.
The `assay` (singular) method will extract the first element of the
`ExperimentList` and will return a `matrix`.


```r
class(assay(miniACC))
#> [1] "matrix"
```

The `assays` (plural) method will return a `SimpleList` of the data with each
element being a `matrix`.


```r
assays(miniACC)
#> List of length 5
#> names(5): RNASeq2GeneNorm gistict RPPAArray Mutations miRNASeqGene
```

*Key point:*

* Whereas the `[[` returned an assay as its original class, `assay()` and
`assays()` convert the assay data to matrix form.

<p style="text-align: right;"> [back to top](#overview-of-key-data-classes) </p>

## Summary of slots and accessors

Slot in the `MultiAssayExperiment` can be accessed or set using their accessor
functions:

| Slot | Accessor |
|------|----------|
| `ExperimentList` | `experiments()`|
| `colData` | `colData()` and `$` * |
| `sampleMap` | `sampleMap()` |
| `metadata` | `metadata()` |

__*__ The `$` operator on a `MultiAssayExperiment` returns a single
column of the `colData`.

## Transformation / reshaping

The `longFormat` or `wideFormat` functions will "reshape" and combine
experiments with each other and with `colData` into one `DataFrame`. These
functions provide compatibility with most of the common R/Bioconductor functions
for regression, machine learning, and visualization.

### `longFormat`

In _long_ format a single column provides all assay results, with additional
optional `colData` columns whose values are repeated as necessary. Here *assay*
is the name of the ExperimentList element, *primary* is the patient identifier
(rowname of colData), *rowname* is the assay rowname (in this case genes),
*colname* is the assay-specific identifier (column name), *value* is the numeric
measurement (gene expression, copy number, presence of a non-silent mutation,
etc), and following these are the *vital_status* and *days_to_death* colData
columns that have been added:


```r
longFormat(miniACC[c("TP53", "CTNNB1"), , ],
           colDataCols = c("vital_status", "days_to_death"))
#> DataFrame with 518 rows and 7 columns
#>               assay      primary     rowname                      colname
#>         <character>  <character> <character>                  <character>
#> 1   RNASeq2GeneNorm TCGA-OR-A5J1        TP53 TCGA-OR-A5J1-01A-11R-A29S-07
#> 2   RNASeq2GeneNorm TCGA-OR-A5J1      CTNNB1 TCGA-OR-A5J1-01A-11R-A29S-07
#> 3   RNASeq2GeneNorm TCGA-OR-A5J2        TP53 TCGA-OR-A5J2-01A-11R-A29S-07
#> 4   RNASeq2GeneNorm TCGA-OR-A5J2      CTNNB1 TCGA-OR-A5J2-01A-11R-A29S-07
#> 5   RNASeq2GeneNorm TCGA-OR-A5J3        TP53 TCGA-OR-A5J3-01A-11R-A29S-07
#> ...             ...          ...         ...                          ...
#> 514       Mutations TCGA-PK-A5HA      CTNNB1 TCGA-PK-A5HA-01A-11D-A29I-10
#> 515       Mutations TCGA-PK-A5HB        TP53 TCGA-PK-A5HB-01A-11D-A29I-10
#> 516       Mutations TCGA-PK-A5HB      CTNNB1 TCGA-PK-A5HB-01A-11D-A29I-10
#> 517       Mutations TCGA-PK-A5HC        TP53 TCGA-PK-A5HC-01A-11D-A30A-10
#> 518       Mutations TCGA-PK-A5HC      CTNNB1 TCGA-PK-A5HC-01A-11D-A30A-10
#>          value vital_status days_to_death
#>      <numeric>    <integer>     <integer>
#> 1     563.4006            1          1355
#> 2    5634.4669            1          1355
#> 3     165.4811            1          1677
#> 4   62658.3913            1          1677
#> 5     956.3028            0            NA
#> ...        ...          ...           ...
#> 514          0            0            NA
#> 515          0            0            NA
#> 516          0            0            NA
#> 517          0            0            NA
#> 518          0            0            NA
```

### `wideFormat`

In _wide_ format, each feature from each assay goes in a separate column, with
one row per primary identifier (patient). Here, each variable becomes a new
column:


```r
wideFormat(miniACC[c("TP53", "CTNNB1"), , ],
           colDataCols = c("vital_status", "days_to_death"))
#> DataFrame with 92 rows and 9 columns
#>          primary vital_status days_to_death RNASeq2GeneNorm_CTNNB1
#>      <character>    <integer>     <integer>              <numeric>
#> 1   TCGA-OR-A5J1            1          1355              5634.4669
#> 2   TCGA-OR-A5J2            1          1677             62658.3913
#> 3   TCGA-OR-A5J3            0            NA              6337.4256
#> 4   TCGA-OR-A5J4            1           423                     NA
#> 5   TCGA-OR-A5J5            1           365               5979.055
#> ...          ...          ...           ...                    ...
#> 88  TCGA-PK-A5H9            0            NA              5258.9863
#> 89  TCGA-PK-A5HA            0            NA              8120.1654
#> 90  TCGA-PK-A5HB            0            NA              5257.8148
#> 91  TCGA-PK-A5HC            0            NA                     NA
#> 92  TCGA-P6-A5OG            1           383              6390.0997
#>     RNASeq2GeneNorm_TP53 gistict_CTNNB1 gistict_TP53 Mutations_CTNNB1
#>                <numeric>      <numeric>    <numeric>        <numeric>
#> 1               563.4006              0            0                0
#> 2               165.4811              1            0                1
#> 3               956.3028              0            0                0
#> 4                     NA              0            1                0
#> 5              1169.6359              0            0                0
#> ...                  ...            ...          ...              ...
#> 88              890.8663              0            0                0
#> 89              683.5722              0           -1                0
#> 90              237.3697             -1           -1                0
#> 91                    NA              1            1                0
#> 92              815.3446              1           -1               NA
#>     Mutations_TP53
#>          <numeric>
#> 1                0
#> 2                1
#> 3                0
#> 4                0
#> 5                0
#> ...            ...
#> 88               0
#> 89               0
#> 90               0
#> 91               0
#> 92              NA
```

## MultiAssayExperiment class construction and concatenation

### MultiAssayExperiment constructor function
The `MultiAssayExperiment` constructor function can take three arguments:

1. `experiments` - An `ExperimentList` or `list` of data
2. `colData` - A `DataFrame` describing the patients (or cell lines, or other
biological units)
3. `sampleMap` - A `DataFrame` of `assay`, `primary`, and `colname` identifiers

The miniACC object can be reconstructed as follows:

```r
MultiAssayExperiment(experiments=experiments(miniACC),
                     colData=colData(miniACC),
                     sampleMap=sampleMap(miniACC),
                     metadata=metadata(miniACC))
#> A MultiAssayExperiment object of 5 listed
#>  experiments with user-defined names and respective classes. 
#>  Containing an ExperimentList class object of length 5: 
#>  [1] RNASeq2GeneNorm: SummarizedExperiment with 198 rows and 79 columns 
#>  [2] gistict: SummarizedExperiment with 198 rows and 90 columns 
#>  [3] RPPAArray: SummarizedExperiment with 33 rows and 46 columns 
#>  [4] Mutations: matrix with 97 rows and 90 columns 
#>  [5] miRNASeqGene: SummarizedExperiment with 471 rows and 80 columns 
#> Features: 
#>  experiments() - obtain the ExperimentList instance 
#>  colData() - the primary/phenotype DataFrame 
#>  sampleMap() - the sample availability DataFrame 
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment 
#>  *Format() - convert into a long or wide DataFrame 
#>  assays() - convert ExperimentList to a SimpleList of matrices
```


### `prepMultiAssay` - Constructor function helper

The `prepMultiAssay` function allows the user to diagnose typical problems
when creating a `MultiAssayExperiment` object. See `?prepMultiAssay` for more
details.

### `c` - concatenate to MultiAssayExperiment

The `c` function allows the user to concatenate an additional experiment to an
existing `MultiAssayExperiment`. The optional `sampleMap` argument allows
concatenating an assay whose column names do not match the row names of
`colData`. For convenience, the _mapFrom_ argument allows the user to map from a
particular experiment **provided** that the **order** of the colnames is in the
**same**. A `warning` will be issued to make the user aware of this assumption.
For example, to concatenate a matrix of log2-transformed RNA-seq results:


```r
miniACC2 <- c(miniACC, log2rnaseq = log2(assays(miniACC)$RNASeq2GeneNorm), mapFrom=1L)
#> Warning in .local(x, ...): Assuming column order in the data provided 
#>  matches the order in 'mapFrom' experiment(s) colnames
experiments(miniACC2)
#> ExperimentList class object of length 6: 
#>  [1] RNASeq2GeneNorm: SummarizedExperiment with 198 rows and 79 columns 
#>  [2] gistict: SummarizedExperiment with 198 rows and 90 columns 
#>  [3] RPPAArray: SummarizedExperiment with 33 rows and 46 columns 
#>  [4] Mutations: matrix with 97 rows and 90 columns 
#>  [5] miRNASeqGene: SummarizedExperiment with 471 rows and 80 columns 
#>  [6] log2rnaseq: matrix with 198 rows and 79 columns
```

<p style="text-align: right;"> [back to top](#overview-of-key-data-classes) </p>

### Building a MultiAssayExperiment from scratch

To start from scratch building your own MultiAssayExperiment, see the package
[Coordinating Analysis of Multi-Assay Experiments vignette](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html).
The package [cheat sheet](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment_cheatsheet.pdf) is also helpful.

If anything is unclear, please ask a question at
https://support.bioconductor.org/ or create an issue on the [MultiAssayExperiment issue tracker](https://github.com/waldronlab/MultiAssayExperiment/issues).

## The Cancer Genome Atlas (TCGA) as MultiAssayExperiments

Most unrestricted TCGA data are available as MultiAssayExperiment objects from
the `curatedTCGAData` package. This represents a lot of harmonization!


```r
library(curatedTCGAData)
curatedTCGAData("ACC")
#>                                  ACC_CNASNP 
#>                   "ACC_CNASNP-20160128.rda" 
#>                                  ACC_CNVSNP 
#>                   "ACC_CNVSNP-20160128.rda" 
#>                        ACC_GISTIC_AllByGene 
#>         "ACC_GISTIC_AllByGene-20160128.rda" 
#>                ACC_GISTIC_ThresholdedByGene 
#> "ACC_GISTIC_ThresholdedByGene-20160128.rda" 
#>                             ACC_Methylation 
#>              "ACC_Methylation-20160128.rda" 
#>                            ACC_miRNASeqGene 
#>             "ACC_miRNASeqGene-20160128.rda" 
#>                                ACC_Mutation 
#>                 "ACC_Mutation-20160128.rda" 
#>                         ACC_RNASeq2GeneNorm 
#>          "ACC_RNASeq2GeneNorm-20160128.rda" 
#>                               ACC_RPPAArray 
#>                "ACC_RPPAArray-20160128.rda"
acc <- curatedTCGAData("ACC",
    assays=c("miRNASeqGene", "RPPAArray", "Mutation", "RNASeq2GeneNorm"),
    dry.run=FALSE)
#> snapshotDate(): 2018-06-29
#> snapshotDate(): 2018-06-29
#> see ?curatedTCGAData and browseVignettes('curatedTCGAData') for documentation
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.ExperimentHub/565'
#> snapshotDate(): 2018-06-29
#> see ?curatedTCGAData and browseVignettes('curatedTCGAData') for documentation
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.ExperimentHub/566'
#> snapshotDate(): 2018-06-29
#> see ?curatedTCGAData and browseVignettes('curatedTCGAData') for documentation
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.ExperimentHub/567'
#> snapshotDate(): 2018-06-29
#> see ?curatedTCGAData and browseVignettes('curatedTCGAData') for documentation
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.ExperimentHub/568'
#> snapshotDate(): 2018-06-29
#> see ?curatedTCGAData and browseVignettes('curatedTCGAData') for documentation
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.ExperimentHub/560'
#> snapshotDate(): 2018-06-29
#> see ?curatedTCGAData and browseVignettes('curatedTCGAData') for documentation
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.ExperimentHub/563'
#> snapshotDate(): 2018-06-29
#> see ?curatedTCGAData and browseVignettes('curatedTCGAData') for documentation
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.ExperimentHub/569'
#> harmonizing input:
#>   removing 620 sampleMap rows not in names(experiments)
acc
#> A MultiAssayExperiment object of 4 listed
#>  experiments with user-defined names and respective classes. 
#>  Containing an ExperimentList class object of length 4: 
#>  [1] ACC_miRNASeqGene-20160128: SummarizedExperiment with 1046 rows and 80 columns 
#>  [2] ACC_Mutation-20160128: RaggedExperiment with 20166 rows and 90 columns 
#>  [3] ACC_RNASeq2GeneNorm-20160128: SummarizedExperiment with 20501 rows and 79 columns 
#>  [4] ACC_RPPAArray-20160128: SummarizedExperiment with 192 rows and 46 columns 
#> Features: 
#>  experiments() - obtain the ExperimentList instance 
#>  colData() - the primary/phenotype DataFrame 
#>  sampleMap() - the sample availability DataFrame 
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment 
#>  *Format() - convert into a long or wide DataFrame 
#>  assays() - convert ExperimentList to a SimpleList of matrices
```

These objects contain most unrestricted TCGA assay and clinical / pathological
data, as well as material curated from the supplements of published TCGA primary
papers at the end of the colData columns:


```r
tail(colnames(colData(acc)), 10)
#>  [1] "MethyLevel"       "miRNA.cluster"    "SCNA.cluster"    
#>  [4] "protein.cluster"  "COC"              "OncoSign"        
#>  [7] "purity"           "ploidy"           "genome_doublings"
#> [10] "ADS"
```

## Plotting, correlation, and other analyses

### How many samples have data for each combination of assays?

**Solution**

The built-in `upsetSamples` creates an "upset" Venn diagram to answer this
question:


```r
upsetSamples(miniACC)
#> Loading required namespace: UpSetR
```

<img src="Ramos_MultiAssayExperiment_files/figure-html/unnamed-chunk-36-1.png" width="672" />

In this dataset only 43 samples have all 5 assays, 32 are missing reverse-phase
protein (RPPAArray), 2 are missing Mutations, 1 is missing gistict, 12 have only
mutations and gistict, etc.

## `TCGAutils` package

Aside from the available reshaping functions already included in the
`MultiAssayExperiment` package, the `TCGAutils` package provides a suite
of helper functions for working with TCGA data. 

### What types of samples are in the data?

**Solution**

The `sampleTables` function gives you an overview of samples in each assay:


```r
library(TCGAutils)
sampleTables(acc)
#> $`ACC_miRNASeqGene-20160128`
#> 
#> 01 
#> 80 
#> 
#> $`ACC_Mutation-20160128`
#> 
#> 01 
#> 90 
#> 
#> $`ACC_RNASeq2GeneNorm-20160128`
#> 
#> 01 
#> 79 
#> 
#> $`ACC_RPPAArray-20160128`
#> 
#> 01 
#> 46

head(sampleTypes)
#>   Code                                      Definition Short.Letter.Code
#> 1   01                             Primary Solid Tumor                TP
#> 2   02                           Recurrent Solid Tumor                TR
#> 3   03 Primary Blood Derived Cancer - Peripheral Blood                TB
#> 4   04    Recurrent Blood Derived Cancer - Bone Marrow              TRBM
#> 5   05                        Additional - New Primary               TAP
#> 6   06                                      Metastatic                TM
```

### Is there subtype data available in the `MultiAssayExperiment` obtained from
`curatedTCGAData`?

**Solution**

The `getSubtypeMap` function will show actual variable names found in `colData`
that contain subtype information. This can only be obtained from
`MultiAssayExperiment` objects provided by `curatedTCGAData`. 


```r
getSubtypeMap(acc)
#>          ACC_annotations     ACC_subtype
#> 1             Patient_ID          SAMPLE
#> 2  histological_subtypes       Histology
#> 3          mrna_subtypes         C1A/C1B
#> 4          mrna_subtypes         mRNA_K4
#> 5                   cimp      MethyLevel
#> 6      microrna_subtypes   miRNA cluster
#> 7          scna_subtypes    SCNA cluster
#> 8       protein_subtypes protein cluster
#> 9   integrative_subtypes             COC
#> 10     mutation_subtypes        OncoSign
head(colData(acc)$Histology)
#> [1] "Usual Type" "Usual Type" "Usual Type" "Usual Type" "Usual Type"
#> [6] "Usual Type"
```

A number of other utility functions are available in `TCGAutils` including
TCGA identifier translation from TCGA barcodes to UUIDs and vice versa and
helper functions to create Bioconductor objects from TCGA exon files
(legacy archive), copy number and GISTIC objects from `RTCGAToolbox`.

### Kaplan-meier plot stratified by pathology_N_stage

Create a Kaplan-meier plot, using pathology_N_stage as a stratifying variable.

**Solution**

The colData provides clinical data for things like a Kaplan-Meier plot for
overall survival stratified by nodal stage.


```r
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})
Surv(miniACC$days_to_death, miniACC$vital_status)
#>  [1] 1355  1677    NA+  423   365    NA+  490   579    NA+  922   551 
#> [12] 1750    NA+ 2105    NA+  541    NA+   NA+  490    NA+   NA+  562 
#> [23]   NA+   NA+   NA+   NA+   NA+   NA+  289    NA+   NA+   NA+  552 
#> [34]   NA+   NA+   NA+  994    NA+   NA+  498    NA+   NA+  344    NA+
#> [45]   NA+   NA+   NA+   NA+   NA+   NA+   NA+   NA+   NA+  391   125 
#> [56]   NA+ 1852    NA+   NA+   NA+   NA+   NA+   NA+   NA+ 1204   159 
#> [67] 1197   662   445    NA+ 2385   436  1105    NA+ 1613    NA+   NA+
#> [78] 2405    NA+   NA+   NA+   NA+   NA+  207     0    NA+   NA+   NA+
#> [89]   NA+   NA+   NA+  383
```

And remove any patients missing overall survival information:

```r
miniACCsurv <- miniACC[, complete.cases(miniACC$days_to_death, miniACC$vital_status), ]
#> harmonizing input:
#>   removing 248 sampleMap rows with 'colname' not in colnames of experiments
#>   removing 58 colData rownames not in sampleMap 'primary'
```


```r
fit <- survfit(Surv(days_to_death, vital_status) ~ pathology_N_stage, data = colData(miniACCsurv))
ggsurvplot(fit, data = colData(miniACCsurv), risk.table = TRUE)
```

<img src="Ramos_MultiAssayExperiment_files/figure-html/unnamed-chunk-41-1.png" width="672" />

### Multivariate Cox regression including RNA-seq, copy number, and pathology

Choose the *EZH2* gene for demonstration. This subsetting will drop assays with
no row named EZH2:

```r
wideacc = wideFormat(miniACC["EZH2", , ],
    colDataCols=c("vital_status", "days_to_death", "pathology_N_stage"))
wideacc$y = Surv(wideacc$days_to_death, wideacc$vital_status)
head(wideacc)
#> DataFrame with 6 rows and 7 columns
#>        primary vital_status days_to_death pathology_N_stage
#>    <character>    <integer>     <integer>       <character>
#> 1 TCGA-OR-A5J1            1          1355                n0
#> 2 TCGA-OR-A5J2            1          1677                n0
#> 3 TCGA-OR-A5J3            0            NA                n0
#> 4 TCGA-OR-A5J4            1           423                n1
#> 5 TCGA-OR-A5J5            1           365                n0
#> 6 TCGA-OR-A5J6            0            NA                n0
#>   RNASeq2GeneNorm_EZH2 gistict_EZH2      y
#>              <numeric>    <numeric> <Surv>
#> 1              75.8886            0 1355:1
#> 2             326.5332            1 1677:1
#> 3              190.194            1   NA:0
#> 4                   NA           -2  423:1
#> 5             366.3826            1  365:1
#> 6              30.7371            1   NA:0
```

Perform a multivariate Cox regression with *EZH2* copy number (gistict),
log2-transformed *EZH2* expression (RNASeq2GeneNorm), and nodal status
(pathology_N_stage) as predictors:

```r
coxph(Surv(days_to_death, vital_status) ~ gistict_EZH2 + log2(RNASeq2GeneNorm_EZH2) + pathology_N_stage,
      data=wideacc)
#> Call:
#> coxph(formula = Surv(days_to_death, vital_status) ~ gistict_EZH2 + 
#>     log2(RNASeq2GeneNorm_EZH2) + pathology_N_stage, data = wideacc)
#> 
#>                               coef exp(coef) se(coef)     z       p
#> gistict_EZH2               -0.0372    0.9635   0.2821 -0.13 0.89499
#> log2(RNASeq2GeneNorm_EZH2)  0.9773    2.6573   0.2811  3.48 0.00051
#> pathology_N_stagen1         0.3775    1.4586   0.5699  0.66 0.50774
#> 
#> Likelihood ratio test=16.28  on 3 df, p=0.001
#> n= 26, number of events= 26 
#>    (66 observations deleted due to missingness)
```

We see that *EZH2* expression is significantly associated with overal survival
(p < 0.001), but *EZH2* copy number and nodal status are not. This analysis
could easily be extended to the whole genome for discovery of prognostic
features by repeated univariate regressions over columns, penalized multivariate
regression, etc.

For further detail, see the main MultiAssayExperiment vignette.

<p style="text-align: right;"> [back to top](#overview-of-key-data-classes) </p>

### Correlation between RNA-seq and copy number

**Part 1**

For all genes where there is both recurrent copy number (gistict assay) and
RNA-seq, calculate the correlation between log2(RNAseq + 1) and copy number.
Create a histogram of these correlations. Compare this with the histogram of
correlations between all *unmatched* gene - copy number pairs.

**Solution**

First, narrow down `miniACC` to only the assays needed:


```r
subacc <- miniACC[, , c("RNASeq2GeneNorm", "gistict")]
```

Align the rows and columns, keeping only samples with both assays available:

```r
subacc <- intersectColumns(subacc)
#> harmonizing input:
#>   removing 15 sampleMap rows with 'colname' not in colnames of experiments
#>   removing 15 colData rownames not in sampleMap 'primary'
subacc <- intersectRows(subacc)
```

Create a list of numeric matrices:

```r
subacc.list <- assays(subacc)
```

Log-transform the RNA-seq assay:

```r
subacc.list[[1]] <- log2(subacc.list[[1]] + 1)
```

Transpose both, so genes are in columns:

```r
subacc.list <- lapply(subacc.list, t)
```

Calculate the correlation between columns in the first matrix and columns in the
second matrix:


```r
corres <- cor(subacc.list[[1]], subacc.list[[2]])
```

And finally, create the histograms:

```r
hist(diag(corres))
```

<img src="Ramos_MultiAssayExperiment_files/figure-html/unnamed-chunk-50-1.png" width="672" />

```r
hist(corres[upper.tri(corres)])
```

<img src="Ramos_MultiAssayExperiment_files/figure-html/unnamed-chunk-50-2.png" width="672" />

**Part 2**

For the gene with highest correlation to copy number, make a box plot of log2
expression against copy number.

**Solution**

First, identify the gene with highest correlation between expression and copy
number:

```r
which.max(diag(corres))
#> EIF4E 
#>    91
```

You could now make the plot by taking the EIF4E columns from each element of the
list subacc.list *list* that was extracted from the subacc
*MultiAssayExperiment*, but let's do it by subsetting and extracting from the
*MultiAssayExperiment*:


```r
df <- wideFormat(subacc["EIF4E", , ])
head(df)
#> DataFrame with 6 rows and 3 columns
#>        primary RNASeq2GeneNorm_EIF4E gistict_EIF4E
#>    <character>             <numeric>     <numeric>
#> 1 TCGA-OR-A5J1              460.6148             0
#> 2 TCGA-OR-A5J2              371.2252             0
#> 3 TCGA-OR-A5J3              516.0717             0
#> 4 TCGA-OR-A5J5             1129.3571             1
#> 5 TCGA-OR-A5J6              822.0782             0
#> 6 TCGA-OR-A5J7              344.5648            -1
```


```r
boxplot(RNASeq2GeneNorm_EIF4E ~ gistict_EIF4E,
        data=df, varwidth=TRUE,
        xlab="GISTIC Relative Copy Number Call",
        ylab="RNA-seq counts")
```

<img src="Ramos_MultiAssayExperiment_files/figure-html/unnamed-chunk-53-1.png" width="672" />

<p style="text-align: right;"> [back to top](#overview-of-key-data-classes) </p>

### Identifying correlated principal components

Perform Principal Components Analysis of each of the five assays, using samples
available on each assay, log-transforming RNA-seq data first.  Using the first
10 components, calculate Pearson correlation between all scores and plot these
correlations as a heatmap to identify correlated components across assays.

**Solution**

Here's a function to simplify doing the PCAs:

```r
getLoadings <- function(x, ncomp=10, dolog=FALSE, center=TRUE, scale.=TRUE){
  if(dolog){
    x <- log2(x + 1)
  }
  pc = prcomp(x, center=center, scale.=scale.)
  return(t(pc$rotation[, 1:10]))
}
```

Although it would be possible to do the following with a loop, the different
data types do require different options for PCA (e.g. mutations are a 0/1 matrix
with 1 meaning there is a somatic mutation, and gistict varies between -2 for
homozygous loss and 2 for a genome doubling, so neither make sense to scale and
center). So it is just as well to do the following one by one, concatenating
each PCA results to the MultiAssayExperiment:


```r
miniACC2 <- intersectColumns(miniACC)
#> harmonizing input:
#>   removing 170 sampleMap rows with 'colname' not in colnames of experiments
#>   removing 49 colData rownames not in sampleMap 'primary'
miniACC2 <- c(miniACC2, rnaseqPCA=getLoadings(assays(miniACC2)[[1]], dolog=TRUE), mapFrom=1L)
#> Warning in .local(x, ...): Assuming column order in the data provided 
#>  matches the order in 'mapFrom' experiment(s) colnames
miniACC2 <- c(miniACC2, gistictPCA=getLoadings(assays(miniACC2)[[2]], center=FALSE, scale.=FALSE), mapFrom=2L)
#> Warning in .local(x, ...): Assuming column order in the data provided 
#>  matches the order in 'mapFrom' experiment(s) colnames
miniACC2 <- c(miniACC2, proteinPCA=getLoadings(assays(miniACC2)[[3]]), mapFrom=3L)
#> Warning in .local(x, ...): Assuming column order in the data provided 
#>  matches the order in 'mapFrom' experiment(s) colnames
miniACC2 <- c(miniACC2, mutationsPCA=getLoadings(assays(miniACC2)[[4]], center=FALSE, scale.=FALSE), mapFrom=4L)
#> Warning in .local(x, ...): Assuming column order in the data provided 
#>  matches the order in 'mapFrom' experiment(s) colnames
miniACC2 <- c(miniACC2, miRNAPCA=getLoadings(assays(miniACC2)[[5]]), mapFrom=5L)
#> Warning in .local(x, ...): Assuming column order in the data provided 
#>  matches the order in 'mapFrom' experiment(s) colnames
```

Now subset to keep *only* the PCA results:

```r
miniACC2 <- miniACC2[, , 6:10]
experiments(miniACC2)
#> ExperimentList class object of length 5: 
#>  [1] rnaseqPCA: matrix with 10 rows and 43 columns 
#>  [2] gistictPCA: matrix with 10 rows and 43 columns 
#>  [3] proteinPCA: matrix with 10 rows and 43 columns 
#>  [4] mutationsPCA: matrix with 10 rows and 43 columns 
#>  [5] miRNAPCA: matrix with 10 rows and 43 columns
```

Note, it would be equally easy (and maybe better) to do PCA on all samples
available for each assay, then do intersectColumns at this point instead.

Now, steps for calculating the correlations and plotting a heatmap:
* Use *wideFormat* to paste these together, which has the nice property of
adding assay names to the column names.
* The first column always contains the sample identifier, so remove it.
* Coerce to a matrix
* Calculate the correlations, and take the absolute value (since signs of
principal components are arbitrary)
* Set the diagonals to NA (each variable has a correlation of 1 to itself).

```r
df <- wideFormat(miniACC2)[, -1]
mycors <- cor(as.matrix(df))
mycors <- abs(mycors)
diag(mycors) <- NA
```

To simplify the heatmap, show only components that have at least one correlation
greater than 0.5.

```r
has.high.cor <- apply(mycors, 2, max, na.rm=TRUE) > 0.5
mycors <- mycors[has.high.cor, has.high.cor]
pheatmap::pheatmap(mycors)
```

<img src="Ramos_MultiAssayExperiment_files/figure-html/unnamed-chunk-58-1.png" width="672" />

The highest correlation present is between PC2 of the RNA-seq assay, and PC1 of
the protein assay.

<p style="text-align: right;"> [back to top](#overview-of-key-data-classes) </p>

### Annotating with ranges

Convert all the `ExperimentList` elements in `miniACC`  to
`RangedSummarizedExperiment` objects.  Then use `rowRanges` to annotate these
objects with genomic ranges. For the microRNA assay, annotate instead with the
genomic coordinates of predicted targets.

**Solution**


The following shortcut function takes a list of human gene symbols and uses
`AnnotationFilter` and `EnsDb.Hsapiens.v86` to look up the ranges, and return
these as a GRangesList which can be used to replace the rowRanges of the
SummarizedExperiment objects:


```r
getrr <- function(identifiers, EnsDbFilterFunc="SymbolFilter"){
  suppressPackageStartupMessages({
    library(AnnotationFilter)
    library(EnsDb.Hsapiens.v86)
  })
    FUN <- get(EnsDbFilterFunc)
    edb <- EnsDb.Hsapiens.v86
    afl <- AnnotationFilterList(FUN(identifiers),
                                SeqNameFilter(c(1:21, "X", "Y")),
                                TxBiotypeFilter("protein_coding"))
    gr <- genes(edb, filter=afl)
    grl <- split(gr, factor(identifiers))
    grl <- grl[match(identifiers, names(grl))]
    stopifnot(identical(names(grl), identifiers))
    return(grl)
}
```

For example:

```r
getrr(rownames(miniACC)[[1]])
#> GRangesList object of length 198:
#> $DIRAS3 
#> GRanges object with 1 range and 7 metadata columns:
#>                   seqnames          ranges strand |         gene_id
#>                      <Rle>       <IRanges>  <Rle> |     <character>
#>   ENSG00000116288        1 7954291-7985505      + | ENSG00000116288
#>                     gene_name   gene_biotype seq_coord_system      symbol
#>                   <character>    <character>      <character> <character>
#>   ENSG00000116288       PARK7 protein_coding       chromosome       PARK7
#>                   entrezid     tx_biotype
#>                     <list>    <character>
#>   ENSG00000116288    11315 protein_coding
#> 
#> $MAPK14 
#> GRanges object with 1 range and 7 metadata columns:
#>                   seqnames          ranges strand |         gene_id
#>   ENSG00000116285        1 8004404-8026308      - | ENSG00000116285
#>                   gene_name   gene_biotype seq_coord_system symbol
#>   ENSG00000116285    ERRFI1 protein_coding       chromosome ERRFI1
#>                   entrezid     tx_biotype
#>   ENSG00000116285    54206 protein_coding
#> 
#> $YAP1 
#> GRanges object with 1 range and 7 metadata columns:
#>                   seqnames            ranges strand |         gene_id
#>   ENSG00000198793        1 11106535-11262507      - | ENSG00000198793
#>                   gene_name   gene_biotype seq_coord_system symbol
#>   ENSG00000198793      MTOR protein_coding       chromosome   MTOR
#>                   entrezid     tx_biotype
#>   ENSG00000198793     2475 protein_coding
#> 
#> ...
#> <195 more elements>
#> -------
#> seqinfo: 22 sequences from GRCh38 genome
```

Use this to set the rowRanges of experiments 1-4 with these GRangesList objects

```r
rseACC <- miniACC
withRSE <- c(1:3, 5)
for (i in withRSE){
  rowRanges(rseACC[[i]]) <- getrr(rownames(rseACC[[i]]))
}
```

Note that the class of experiments 1-4 is now `RangedSummarizedExperiment`:

```r
experiments(rseACC)
#> ExperimentList class object of length 5: 
#>  [1] RNASeq2GeneNorm: RangedSummarizedExperiment with 198 rows and 79 columns 
#>  [2] gistict: RangedSummarizedExperiment with 198 rows and 90 columns 
#>  [3] RPPAArray: RangedSummarizedExperiment with 33 rows and 46 columns 
#>  [4] Mutations: matrix with 97 rows and 90 columns 
#>  [5] miRNASeqGene: RangedSummarizedExperiment with 471 rows and 80 columns
```

With ranged objects in the MultiAssayExperiment, you can then do subsetting by
ranges. For example, select all genes on chromosome 1 for the four
*rangedSummarizedExperiment* objects:

```r
rseACC[GRanges(seqnames="1:1-1e9"), , withRSE]
#> A MultiAssayExperiment object of 3 listed
#>  experiments with user-defined names and respective classes. 
#>  Containing an ExperimentList class object of length 3: 
#>  [1] RNASeq2GeneNorm: RangedSummarizedExperiment with 22 rows and 79 columns 
#>  [2] gistict: RangedSummarizedExperiment with 22 rows and 90 columns 
#>  [3] RPPAArray: RangedSummarizedExperiment with 3 rows and 46 columns 
#> Features: 
#>  experiments() - obtain the ExperimentList instance 
#>  colData() - the primary/phenotype DataFrame 
#>  sampleMap() - the sample availability DataFrame 
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment 
#>  *Format() - convert into a long or wide DataFrame 
#>  assays() - convert ExperimentList to a SimpleList of matrices
```

*Note about microRNA*: You can set ranges for the microRNA assay according to
the genomic location of those microRNA, or the locations of their predicted
targets, but we don't do it here. Assigning genomic ranges of microRNA targets
would be an easy way to subset them according to their targets.

<p style="text-align: right;"> [back to top](#overview-of-key-data-classes) </p>
