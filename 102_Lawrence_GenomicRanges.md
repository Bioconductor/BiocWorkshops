
# 102: Solving common bioinformatic challenges using GenomicRanges

## Instructor name and contact information

* Michael Lawrence (michafla@gene.com)

## Workshop Description

We will introduce the fundamental concepts underlying the
GenomicRanges package and related infrastructure. After a structured
introduction, we will follow a realistic workflow, along the way
exploring the central data structures, including GRanges and
SummarizedExperiment, and useful operations in the ranges
algebra. Topics will include data import/export, computing and
summarizing data on genomic features, overlap detection, integration
with reference annotations, scaling strategies, and
visualization. Students can follow along, and there will be plenty of
time for students to ask questions about how to apply the
infrastructure to their particular use case. Michael Lawrence
(Genentech).

### Pre-requisites

* Solid understanding of R
* Basic familiarity with GRanges objects
* Basic familiarity with packages like S4Vectors, IRanges,
  GenomicRanges, rtracklayer, etc.

### Workshop Participation

Describe how students will be expected to participate in the workshop.

### _R_ / _Bioconductor_ packages used

* S4Vectors
* IRanges
* GenomicRanges
* rtracklayer
* GenomicFeatures
* SummarizedExperiment
* GenomicAlignments

### Time outline

| Activity                     | Time |
|------------------------------|------|
| Intro slides                 | 30m  |
| Workflow(s)                  | 1hr  |
| Remaining questions          | 30m  |

## Workshop goals and objectives

### Learning goals

 * Understand how to apply the *Ranges infrastructure to real-world
   problems
 * Gain insight into the design principles of the infrastructure and
   how it was meant to be used

### Learning objectives

* Manipulate GRanges and related objects
* Use the ranges algebra to analyze genomic ranges
* Implement efficient workflows based on the *Ranges infrastructure

## Introduction

### What is the Ranges infrastructure?

The Ranges framework of packages provide data structures and
algorithms for analyzing genomic data. This includes standard genomic
data containers like GRanges and SummarizedExperiment, optimized data
representations like Rle, and fast algorithms for computing overlaps,
finding nearest neighbors, summarizing ranges and metadata, etc.

### Why use the Ranges infrastructure?

Hundreds of Bioconductor packages operate on Ranges data structures,
enabling the construction of complex workflows integrating multiple
packages and data types. The API directly supports data analysis as
well the construction of new genomic software. Code evolves easily
from analysis script to generalized package extending the Bioconductor
ecosystem.

### Who is this workshop for?

If you still think of R as a programming language and want to write
new bioinformatics algorithms and/or build interoperable software on
top of formal genomic data structures, this workshop is for you. For
the tidyverse analog of this workshop, see the plyranges tutorial by
Stuart Lee.

## Setup

To participate in this workshop you'll need to have R >= 3.5 and install
the GenomicRanges, AnnotationHub, and airway Bioconductor 3.7 packages 
(@R-AnnotationHub; @R-airway). You can achieve this
by installing the BiocManager package from CRAN, loading it then running the
install command:


```r
install.packages("BiocManager")
library(BiocManager)
install(c("GenomicRanges", "AnnotationHub", "airway"))
```

## *GRanges*: Genomic Ranges

<div class="figure">
<img src="Lawrence_GenomicRanges/granges.pdf" alt="An illustration of genomic ranges. GRanges represents a set genomic ranges in terms of the sequence name (typically the chromosome), start and end coordinates (as an IRanges object), and strand (either positive, negative, or unstranded). GRanges holds information about its universe of sequences (typically a genome) and an arbitrary set of metadata columns with information particular to the dataset." width="\textwidth" />
<p class="caption">(\#fig:GRanges)An illustration of genomic ranges. GRanges represents a set genomic ranges in terms of the sequence name (typically the chromosome), start and end coordinates (as an IRanges object), and strand (either positive, negative, or unstranded). GRanges holds information about its universe of sequences (typically a genome) and an arbitrary set of metadata columns with information particular to the dataset.</p>
</div>

The central genomic data structure is the *GRanges* class, 
which represents a collection of genomic ranges
that each have a single start and end location on the genome. It can be
used to store the location of genomic features such as binding
sites, read alignments and transcripts. 

## Constructing a *GRanges* object from data.frame

If we have a data.frame containing scores on a set of genomic
ranges, we can call `makeGRangesFromDataFrame()` to promote the
data.frame to a GRanges, thus adding semantics, formal constraints,
and range-specific functionality. For example,


```r
suppressPackageStartupMessages({
 library(BiocStyle)
 library(GenomicRanges)
})
```


```r
df <- data.frame(
    seqnames = rep(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
    start = c(101, 105, 125, 132, 134, 152, 153, 160, 166, 170),
    end = c(104, 120, 133, 132, 155, 154, 159, 166, 171, 190),
    strand = rep(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    score = 1:10,
    GC = seq(1, 0, length=10),
    row.names = head(letters, 10))
gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
```

creates a *GRanges* object with 10 genomic ranges.
The output of the *GRanges* `show()` method separates the
information into a left and right hand region that are separated by
`|` symbols. The genomic coordinates (seqnames, ranges, and strand)
are located on the left-hand side and the metadata columns (annotation)
are located on the right. For this example, the metadata is
comprised of `"score"` and `"GC"` information, but almost
anything can be stored in the metadata portion of a *GRanges*
object.

## Loading a *GRanges* object from a standard file format

We often obtain data on genomic ranges from standard track formats,
like BED, GFF and BigWig. The rtracklayer package parses those files
directly into GRanges objects. The GenomicAlignments package parses
BAM files into GAlignments objects, which behave much like GRanges,
and it is easy to convert a GAlignments to a GRanges. We will see some
examples of loading data from files later in the tutorial.

The `seqnames()`, `ranges()`, and `strand()` accessor functions
extract the components of the genomic coordinates,
 
 
## Basic manipulation of *GRanges* objects


```r
seqnames(gr)
#> factor-Rle of length 10 with 4 runs
#>   Lengths:    1    3    2    4
#>   Values : chr1 chr2 chr1 chr3
#> Levels(3): chr1 chr2 chr3
ranges(gr)
#> IRanges object with 10 ranges and 0 metadata columns:
#>         start       end     width
#>     <integer> <integer> <integer>
#>   a       101       104         4
#>   b       105       120        16
#>   c       125       133         9
#>   d       132       132         1
#>   e       134       155        22
#>   f       152       154         3
#>   g       153       159         7
#>   h       160       166         7
#>   i       166       171         6
#>   j       170       190        21
strand(gr)
#> factor-Rle of length 10 with 5 runs
#>   Lengths: 1 2 2 3 2
#>   Values : - + * + -
#> Levels(3): + - *
```

The `granges()` function extracts genomic ranges without corresponding
metadata,


```r
granges(gr)
#> GRanges object with 10 ranges and 0 metadata columns:
#>     seqnames    ranges strand
#>        <Rle> <IRanges>  <Rle>
#>   a     chr1   101-104      -
#>   b     chr2   105-120      +
#>   c     chr2   125-133      +
#>   d     chr2       132      *
#>   e     chr1   134-155      *
#>   f     chr1   152-154      +
#>   g     chr3   153-159      +
#>   h     chr3   160-166      +
#>   i     chr3   166-171      -
#>   j     chr3   170-190      -
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

The `start()`, `end()`, `width()`, and `range` functions extract basic
interval characteristics,


```r
start(gr)
#>  [1] 101 105 125 132 134 152 153 160 166 170
end(gr)
#>  [1] 104 120 133 132 155 154 159 166 171 190
width(gr)
#>  [1]  4 16  9  1 22  3  7  7  6 21
```

The `mcols()` accessor extracts the metadata as a *DataFrame*,
 

```r
mcols(gr)
#> DataFrame with 10 rows and 2 columns
#>       score                GC
#>   <integer>         <numeric>
#> a         1                 1
#> b         2 0.888888888888889
#> c         3 0.777777777777778
#> d         4 0.666666666666667
#> e         5 0.555555555555556
#> f         6 0.444444444444444
#> g         7 0.333333333333333
#> h         8 0.222222222222222
#> i         9 0.111111111111111
#> j        10                 0
mcols(gr)$score
#>  [1]  1  2  3  4  5  6  7  8  9 10
score(gr)
#>  [1]  1  2  3  4  5  6  7  8  9 10
```

The lengths and other properties of the sequences containing the
ranges can (and should) be stored in the *GRanges* object. Formal
tracking of the sequence universe, typically the genome build, ensures
data integrity and prevents accidental mixing of ranges from
incompatible contexts. Assuming these data are of *Homo sapiens*, we
could add the sequence information like this:
 

```r
seqinfo(gr) <- Seqinfo(genome="hg38")
```
The `Seqinfo()` function automatically loads the sequence information
for the specified `genome=` by querying the UCSC database.
 
And then retrieves as:

```r
seqinfo(gr)
#> Seqinfo object with 455 sequences (1 circular) from hg38 genome:
#>   seqnames         seqlengths isCircular genome
#>   chr1              248956422      FALSE   hg38
#>   chr2              242193529      FALSE   hg38
#>   chr3              198295559      FALSE   hg38
#>   chr4              190214555      FALSE   hg38
#>   chr5              181538259      FALSE   hg38
#>   ...                     ...        ...    ...
#>   chrUn_KI270753v1      62944      FALSE   hg38
#>   chrUn_KI270754v1      40191      FALSE   hg38
#>   chrUn_KI270755v1      36723      FALSE   hg38
#>   chrUn_KI270756v1      79590      FALSE   hg38
#>   chrUn_KI270757v1      71251      FALSE   hg38
```
 
Methods for accessing the `length` and `names` have
also been defined.
 

```r
names(gr)
#>  [1] "a" "b" "c" "d" "e" "f" "g" "h" "i" "j"
length(gr)
#> [1] 10
```

## Subsetting  *GRanges* objects

*GRanges* objects act like vectors of ranges, with the expected
vector-like subsetting operations available
 

```r
gr[2:3]
#> GRanges object with 2 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   c     chr2   125-133      + |         3 0.777777777777778
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```
 
A second argument to the `[` subset operator specifies which metadata
columns to extract from the *GRanges* object. For example,
 

```r
gr[2:3, "GC"]
#> GRanges object with 2 ranges and 1 metadata column:
#>     seqnames    ranges strand |                GC
#>        <Rle> <IRanges>  <Rle> |         <numeric>
#>   b     chr2   105-120      + | 0.888888888888889
#>   c     chr2   125-133      + | 0.777777777777778
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

The `subset()` function provides an easy way to subset based on
attributes of the ranges and columns in the metadata. For example,

```r
subset(gr, strand == "+" & score > 5, select = score)
#> GRanges object with 3 ranges and 1 metadata column:
#>     seqnames    ranges strand |     score
#>        <Rle> <IRanges>  <Rle> | <integer>
#>   f     chr1   152-154      + |         6
#>   g     chr3   153-159      + |         7
#>   h     chr3   160-166      + |         8
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

Elements can also be assigned to the *GRanges* object.  This example
replaces the the second row of a *GRanges* object with the first row
of `gr`.
 

```r
grMod <- gr
grMod[2] <- gr[1]
head(grMod, n=3)
#> GRanges object with 3 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1   101-104      - |         1                 1
#>   b     chr1   101-104      - |         1                 1
#>   c     chr2   125-133      + |         3 0.777777777777778
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

There are methods to repeat, reverse, or select specific portions of
*GRanges* objects.
 

```r
rep(gr[2], times = 3)
#> GRanges object with 3 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
rev(gr)
#> GRanges object with 10 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   j     chr3   170-190      - |        10                 0
#>   i     chr3   166-171      - |         9 0.111111111111111
#>   h     chr3   160-166      + |         8 0.222222222222222
#>   g     chr3   153-159      + |         7 0.333333333333333
#>   f     chr1   152-154      + |         6 0.444444444444444
#>   e     chr1   134-155      * |         5 0.555555555555556
#>   d     chr2       132      * |         4 0.666666666666667
#>   c     chr2   125-133      + |         3 0.777777777777778
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   a     chr1   101-104      - |         1                 1
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
head(gr,n=2)
#> GRanges object with 2 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1   101-104      - |         1                 1
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
tail(gr,n=2)
#> GRanges object with 2 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   i     chr3   166-171      - |         9 0.111111111111111
#>   j     chr3   170-190      - |        10                 0
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
window(gr, start=2,end=4)
#> GRanges object with 3 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   c     chr2   125-133      + |         3 0.777777777777778
#>   d     chr2       132      * |         4 0.666666666666667
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
gr[IRanges(start=c(2,7), end=c(3,9))]
#> GRanges object with 5 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   c     chr2   125-133      + |         3 0.777777777777778
#>   g     chr3   153-159      + |         7 0.333333333333333
#>   h     chr3   160-166      + |         8 0.222222222222222
#>   i     chr3   166-171      - |         9 0.111111111111111
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

## Splitting and combining *GRanges* objects

THe `split()` function divides a *GRanges* into groups, returning a
*GRangesList*, a class that we will describe and demonstrate later.
 

```r
sp <- split(gr, rep(1:2, each=5))
sp
#> GRangesList object of length 2:
#> $1 
#> GRanges object with 5 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1   101-104      - |         1                 1
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   c     chr2   125-133      + |         3 0.777777777777778
#>   d     chr2       132      * |         4 0.666666666666667
#>   e     chr1   134-155      * |         5 0.555555555555556
#> 
#> $2 
#> GRanges object with 5 ranges and 2 metadata columns:
#>     seqnames  ranges strand | score                GC
#>   f     chr1 152-154      + |     6 0.444444444444444
#>   g     chr3 153-159      + |     7 0.333333333333333
#>   h     chr3 160-166      + |     8 0.222222222222222
#>   i     chr3 166-171      - |     9 0.111111111111111
#>   j     chr3 170-190      - |    10                 0
#> 
#> -------
#> seqinfo: 455 sequences (1 circular) from hg38 genome
```

We can split the ranges by metadata columns, like strand,


```r
split(gr, ~ strand)
#> GRangesList object of length 3:
#> $+ 
#> GRanges object with 5 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   c     chr2   125-133      + |         3 0.777777777777778
#>   f     chr1   152-154      + |         6 0.444444444444444
#>   g     chr3   153-159      + |         7 0.333333333333333
#>   h     chr3   160-166      + |         8 0.222222222222222
#> 
#> $- 
#> GRanges object with 3 ranges and 2 metadata columns:
#>     seqnames  ranges strand | score                GC
#>   a     chr1 101-104      - |     1                 1
#>   i     chr3 166-171      - |     9 0.111111111111111
#>   j     chr3 170-190      - |    10                 0
#> 
#> $* 
#> GRanges object with 2 ranges and 2 metadata columns:
#>     seqnames  ranges strand | score                GC
#>   d     chr2     132      * |     4 0.666666666666667
#>   e     chr1 134-155      * |     5 0.555555555555556
#> 
#> -------
#> seqinfo: 455 sequences (1 circular) from hg38 genome
```

The `c()` and `append()` functions combine two (or more in the case of
`c()`) *GRanges* objects.
 

```r
c(sp[[1]], sp[[2]])
#> GRanges object with 10 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1   101-104      - |         1                 1
#>   b     chr2   105-120      + |         2 0.888888888888889
#>   c     chr2   125-133      + |         3 0.777777777777778
#>   d     chr2       132      * |         4 0.666666666666667
#>   e     chr1   134-155      * |         5 0.555555555555556
#>   f     chr1   152-154      + |         6 0.444444444444444
#>   g     chr3   153-159      + |         7 0.333333333333333
#>   h     chr3   160-166      + |         8 0.222222222222222
#>   i     chr3   166-171      - |         9 0.111111111111111
#>   j     chr3   170-190      - |        10                 0
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

The `stack()` function stacks the elements of a *GRangesList* into a
single *GRanges* and adds a column indicating the origin of each
element,

```r
stack(sp, index.var="group")
#> GRanges object with 10 ranges and 3 metadata columns:
#>     seqnames    ranges strand | group     score                GC
#>        <Rle> <IRanges>  <Rle> | <Rle> <integer>         <numeric>
#>   a     chr1   101-104      - |     1         1                 1
#>   b     chr2   105-120      + |     1         2 0.888888888888889
#>   c     chr2   125-133      + |     1         3 0.777777777777778
#>   d     chr2       132      * |     1         4 0.666666666666667
#>   e     chr1   134-155      * |     1         5 0.555555555555556
#>   f     chr1   152-154      + |     2         6 0.444444444444444
#>   g     chr3   153-159      + |     2         7 0.333333333333333
#>   h     chr3   160-166      + |     2         8 0.222222222222222
#>   i     chr3   166-171      - |     2         9 0.111111111111111
#>   j     chr3   170-190      - |     2        10                 0
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

## Aggregating *GRanges* objects

Like other tabular data structures, we can aggregate *GRanges*
objects, for example,


```r
aggregate(gr, score ~ strand, mean)
#> DataFrame with 3 rows and 2 columns
#>     strand            score
#>   <factor>        <numeric>
#> 1        +              5.2
#> 2        - 6.66666666666667
#> 3        *              4.5
```

The `aggregate()` function also supports a syntax similar to
`summarize()` from dplyr,


```r
aggregate(gr, ~ strand, n_score = lengths(score), mean_score = mean(score))
#> DataFrame with 3 rows and 4 columns
#>              grouping   strand   n_score       mean_score
#>   <ManyToOneGrouping> <factor> <integer>        <numeric>
#> 1           2,3,6,...        +         5              5.2
#> 2              1,9,10        -         3 6.66666666666667
#> 3                 4,5        *         2              4.5
```

Note that we need to call `lengths(score)` instead of `length(score)`
because `score` is actually a list-like object in the aggregation
expression.

## Basic interval operations for *GRanges* objects

There are many functions for manipulating *GRanges* objects. The
functions can be classified as *intra-range functions*, *inter-range
functions*, and *between-range functions*.

*Intra-range functions* operate on each element of a
*GRanges* object independent of the other ranges in the
object. For example, the `flank` function can be used to recover
regions flanking the set of ranges represented by the *GRanges*
object. So to get a *GRanges* object containing the ranges that
include the 10 bases upstream according to the direction of
"transcription" (indicated by the strand):
 
 
 ```r
 g <- gr[1:3]
 g <- append(g, gr[10])
 flank(g, 10)
 #> GRanges object with 4 ranges and 2 metadata columns:
 #>     seqnames    ranges strand |     score                GC
 #>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
 #>   a     chr1   105-114      - |         1                 1
 #>   b     chr2    95-104      + |         2 0.888888888888889
 #>   c     chr2   115-124      + |         3 0.777777777777778
 #>   j     chr3   191-200      - |        10                 0
 #>   -------
 #>   seqinfo: 455 sequences (1 circular) from hg38 genome
 ```

And to include the downstream bases:
 

```r
flank(g, 10, start=FALSE)
#> GRanges object with 4 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1    91-100      - |         1                 1
#>   b     chr2   121-130      + |         2 0.888888888888889
#>   c     chr2   134-143      + |         3 0.777777777777778
#>   j     chr3   160-169      - |        10                 0
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

A common use case for `flank()` is generating promoter regions based
on the transcript ranges. There is a convenience function that by
default generates a region starting 2000bp upstream and 200bp
downstream of the TSS,


```r
promoters(g)
#> Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 4 out-of-bound ranges located on sequences
#>   chr1, chr2, and chr3. Note that ranges located on a sequence whose
#>   length is unknown (NA) or on a circular sequence are not
#>   considered out-of-bound (use seqlengths() and isCircular() to get
#>   the lengths and circularity flags of the underlying sequences).
#>   You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> GRanges object with 4 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1  -95-2104      - |         1                 1
#>   b     chr2 -1895-304      + |         2 0.888888888888889
#>   c     chr2 -1875-324      + |         3 0.777777777777778
#>   j     chr3   -9-2190      - |        10                 0
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

To ignore strand/transcription and assume the orientation of left to
right use `unstrand()`,

```r
flank(unstrand(g), 10)
#> GRanges object with 4 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1    91-100      * |         1                 1
#>   b     chr2    95-104      * |         2 0.888888888888889
#>   c     chr2   115-124      * |         3 0.777777777777778
#>   j     chr3   160-169      * |        10                 0
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

Other examples of intra-range functions include `resize()` and
`shift()`. The `shift()` function will move the ranges by a specific number
of base pairs, and the `resize()` function will set a specific width, by
default fixing the "transcription" start (or just the start when
strand is "*"). The `fix=` argument controls whether the "start",
"end" or "center" is held constant.
 

```r
shift(g, 5)
#> GRanges object with 4 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1   106-109      - |         1                 1
#>   b     chr2   110-125      + |         2 0.888888888888889
#>   c     chr2   130-138      + |         3 0.777777777777778
#>   j     chr3   175-195      - |        10                 0
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
resize(g, 30)
#> GRanges object with 4 ranges and 2 metadata columns:
#>     seqnames    ranges strand |     score                GC
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric>
#>   a     chr1    75-104      - |         1                 1
#>   b     chr2   105-134      + |         2 0.888888888888889
#>   c     chr2   125-154      + |         3 0.777777777777778
#>   j     chr3   161-190      - |        10                 0
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```
 
The *[GenomicRanges](http://bioconductor.org/packages/GenomicRanges)* help page `?"intra-range-methods"`
summarizes these methods.

*Inter-range functions* involve comparisons between ranges in a single
*GRanges* object and typically aggregate ranges. For instance, the
`reduce()` function will merge overlapping and adjacent ranges to
produce a minimal set of ranges representing the regions covered by
the original set.
 

```r
reduce(gr)
#> GRanges object with 8 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1   152-154      +
#>   [2]     chr1   101-104      -
#>   [3]     chr1   134-155      *
#>   [4]     chr2   105-120      +
#>   [5]     chr2   125-133      +
#>   [6]     chr2       132      *
#>   [7]     chr3   153-166      +
#>   [8]     chr3   166-190      -
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
reduce(gr, ignore.strand=TRUE)
#> GRanges object with 5 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1   101-104      *
#>   [2]     chr1   134-155      *
#>   [3]     chr2   105-120      *
#>   [4]     chr2   125-133      *
#>   [5]     chr3   153-190      *
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```
 
Rarely, it useful to complement the (reduced) ranges. Note that the
universe is taken as the entire sequence span in all three strands (+,
-, *), which is often surprising when working with unstranded ranges.

```r
gaps(g)
#> GRanges object with 1369 ranges and 0 metadata columns:
#>                  seqnames        ranges strand
#>                     <Rle>     <IRanges>  <Rle>
#>      [1]             chr1   1-248956422      +
#>      [2]             chr1         1-100      -
#>      [3]             chr1 105-248956422      -
#>      [4]             chr1   1-248956422      *
#>      [5]             chr2         1-104      +
#>      ...              ...           ...    ...
#>   [1365] chrUn_KI270756v1       1-79590      -
#>   [1366] chrUn_KI270756v1       1-79590      *
#>   [1367] chrUn_KI270757v1       1-71251      +
#>   [1368] chrUn_KI270757v1       1-71251      -
#>   [1369] chrUn_KI270757v1       1-71251      *
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```
 
The `disjoin` function breaks up the ranges so that they do not
overlap but still cover the same regions:
 

```r
disjoin(g)
#> GRanges object with 4 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1   101-104      -
#>   [2]     chr2   105-120      +
#>   [3]     chr2   125-133      +
#>   [4]     chr3   170-190      -
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```
 
The `coverage` function counts how many ranges overlap each position
in the sequence universe of a *GRanges* object.
 

```r
cov <- coverage(g)
cov[1:3]
#> RleList of length 3
#> $chr1
#> integer-Rle of length 248956422 with 3 runs
#>   Lengths:       100         4 248956318
#>   Values :         0         1         0
#> 
#> $chr2
#> integer-Rle of length 242193529 with 5 runs
#>   Lengths:       104        16         4         9 242193396
#>   Values :         0         1         0         1         0
#> 
#> $chr3
#> integer-Rle of length 198295559 with 3 runs
#>   Lengths:       169        21 198295369
#>   Values :         0         1         0
```
The coverage is stored compactly as an *RleList*, with one *Rle*
vector per sequence. We can convert it to a *GRanges*,

```r
cov_gr <- GRanges(cov)
cov_gr
#> GRanges object with 463 ranges and 1 metadata column:
#>                 seqnames        ranges strand |     score
#>                    <Rle>     <IRanges>  <Rle> | <integer>
#>     [1]             chr1         1-100      * |         0
#>     [2]             chr1       101-104      * |         1
#>     [3]             chr1 105-248956422      * |         0
#>     [4]             chr2         1-104      * |         0
#>     [5]             chr2       105-120      * |         1
#>     ...              ...           ...    ... .       ...
#>   [459] chrUn_KI270753v1       1-62944      * |         0
#>   [460] chrUn_KI270754v1       1-40191      * |         0
#>   [461] chrUn_KI270755v1       1-36723      * |         0
#>   [462] chrUn_KI270756v1       1-79590      * |         0
#>   [463] chrUn_KI270757v1       1-71251      * |         0
#>   -------
#>   seqinfo: 455 sequences from an unspecified genome
```
and even convert the *GRanges* form back to an *RleList* by computing
a weighted coverage,

```r
cov <- coverage(cov_gr, weight="score")
```

The *GRanges* derivative *GPos*, a compact representation of width 1
ranges, is useful for representing coverage, although it cannot yet
represent the coverage for the entire human genome (or any genome with
over ~ 2 billion bp).

```r
GPos(cov[1:3])
#> GPos object with 689445510 positions and 0 metadata columns:
#>               seqnames       pos strand
#>                  <Rle> <integer>  <Rle>
#>           [1]     chr1         1      *
#>           [2]     chr1         2      *
#>           [3]     chr1         3      *
#>           [4]     chr1         4      *
#>           [5]     chr1         5      *
#>           ...      ...       ...    ...
#>   [689445506]     chr3 198295555      *
#>   [689445507]     chr3 198295556      *
#>   [689445508]     chr3 198295557      *
#>   [689445509]     chr3 198295558      *
#>   [689445510]     chr3 198295559      *
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome
```

These inter-range functions all generate entirely new sets of
ranges. The return value is left unannotated, since there is no
obvious way to carry the metadata across the operation. The user is
left to map the metadata to the new ranges. Functions like `reduce()`
and `disjoin()` facilitate this by optionally including in the
returned metadata a one-to-many reverse mapping from the aggregate
ranges to input ranges. For example, to average the score over a
reduction,

```r
rg <- reduce(gr, with.revmap=TRUE)
rg$score <- mean(extractList(gr$score, rg$revmap))
```

See the *[GenomicRanges](http://bioconductor.org/packages/GenomicRanges)* help page
`?"inter-range-methods"` for additional help.

## Interval set operations for *GRanges* objects

*Between-range functions* calculate relationships between different
*GRanges* objects. Of central importance are
`findOverlaps` and related operations; these are discussed
below.  Additional operations treat *GRanges* as mathematical
sets of coordinates; `union(g, g2)` is the union of the
coordinates in `g` and `g2`. Here are examples for
calculating the `union`, the `intersect` and the
asymmetric difference (using `setdiff`).
 

```r
g2 <- head(gr, n=2)
union(g, g2)
#> GRanges object with 4 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1   101-104      -
#>   [2]     chr2   105-120      +
#>   [3]     chr2   125-133      +
#>   [4]     chr3   170-190      -
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
intersect(g, g2)
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1   101-104      -
#>   [2]     chr2   105-120      +
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
setdiff(g, g2)
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr2   125-133      +
#>   [2]     chr3   170-190      -
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

Related functions are available when the structure of the
*GRanges* objects are 'parallel' to one another, i.e., element
1 of object 1 is related to element 1 of object 2, and so on. These
operations all begin with a `p`, which is short for
parallel. The functions then perform element-wise, e.g., the union of
element 1 of object 1 with element 1 of object 2, etc. A requirement
for these operations is that the number of elements in each
*GRanges* object is the same, and that both of the objects have
the same seqnames and strand assignments throughout.
 

```r
g3 <- g[1:2]
ranges(g3[1]) <- IRanges(start=105, end=112)
punion(g2, g3)
#> GRanges object with 2 ranges and 0 metadata columns:
#>     seqnames    ranges strand
#>        <Rle> <IRanges>  <Rle>
#>   a     chr1   101-112      -
#>   b     chr2   105-120      +
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
pintersect(g2, g3)
#> GRanges object with 2 ranges and 3 metadata columns:
#>     seqnames    ranges strand |     score                GC       hit
#>        <Rle> <IRanges>  <Rle> | <integer>         <numeric> <logical>
#>   a     chr1   105-104      - |         1                 1      TRUE
#>   b     chr2   105-120      + |         2 0.888888888888889      TRUE
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
psetdiff(g2, g3)
#> GRanges object with 2 ranges and 0 metadata columns:
#>     seqnames    ranges strand
#>        <Rle> <IRanges>  <Rle>
#>   a     chr1   101-104      -
#>   b     chr2   105-104      +
#>   -------
#>   seqinfo: 455 sequences (1 circular) from hg38 genome
```

For more information on the `GRanges` classes be sure to consult
the manual page.
 

```r
?GRanges
```
 
A relatively comprehensive list of available functions is discovered
with
 

```r
methods(class="GRanges")
```


## Finding overlaps between *GRanges* objects

Interval overlapping is the process of comparing the ranges in two
objects to determine if and when they overlap. As such, it is perhaps
the most common operation performed on *GRanges* objects. 
To this end, the *[GenomicRanges](http://bioconductor.org/packages/GenomicRanges)*
package provides a family of interval overlap functions. The most general
of these functions is `findOverlaps()`, which takes a query and a
subject as inputs and returns a *Hits* object containing
the index pairings for the overlapping elements.

Let us assume that we have three random data.frame objects, each with
annoyingly differing ways of naming the columns defining the ranges,

```r
set.seed(66+105+111+99+49+56)

pos <- sample(1:200, size = 30L)
size <- 10L
end <- size + pos - 1L
chrom <- sample(paste0("chr", 1:3), size = 30L, replace = TRUE)
query_df <- data.frame(chrom = chrom, 
                       start = pos,
                       end = end)
query_dfs <- split(query_df, 1:3)
q1 <- rename(query_dfs[[1L]], start = "pos")
q2 <- rename(query_dfs[[2L]], chrom = "ch", start = "st")
q3 <- rename(query_dfs[[3L]], end = "last")
```
The `makeGRangesFromDataFrame()` function can guess some of these, but
not all of them, so we help it out,

```r
q1 <- makeGRangesFromDataFrame(q1, start.field = "pos")
q2 <- makeGRangesFromDataFrame(q2, seqnames.field = "ch",
                                 start.field = "st")
q3 <- makeGRangesFromDataFrame(q3, end.field = "last")
query <- mstack(q1, q2, q3, .index.var="replicate")
sort(query, by = ~ start)
#> GRanges object with 30 ranges and 1 metadata column:
#>      seqnames    ranges strand | replicate
#>         <Rle> <IRanges>  <Rle> |     <Rle>
#>   25     chr1     11-20      * |         1
#>   22     chr1     16-25      * |         1
#>    2     chr2     21-30      * |         2
#>   30     chr3     50-59      * |         3
#>    6     chr2     51-60      * |         3
#>   ..      ...       ...    ... .       ...
#>   26     chr3   160-169      * |         2
#>   13     chr1   169-178      * |         1
#>   29     chr2   183-192      * |         2
#>   27     chr3   190-199      * |         3
#>   24     chr3   197-206      * |         3
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
```
Above, we use the convenient `mstack()` function, which stacks its
arguments, populating the `.index.var=` column with the origin of each
range (using the argument names or positions).

Perhaps the simplest overlap-based operation is `subsetByOverlaps()`,
which extracts the elements in the query (the first argument) that
overlap at least one element in the subject (the second).


```r
subject <- gr
subsetByOverlaps(query, subject, ignore.strand=TRUE)
#> GRanges object with 9 ranges and 1 metadata column:
#>      seqnames    ranges strand | replicate
#>         <Rle> <IRanges>  <Rle> |     <Rle>
#>   10     chr2   120-129      * |         1
#>   28     chr1   151-160      * |         1
#>    5     chr2   128-137      * |         2
#>   14     chr1   149-158      * |         2
#>   23     chr3   153-162      * |         2
#>   26     chr3   160-169      * |         2
#>    9     chr1    92-101      * |         3
#>   21     chr1   148-157      * |         3
#>   27     chr3   190-199      * |         3
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
```
In every call to an overlap operation, it is necessary to specify
`ignore.strand=TRUE`, except in rare cases when we do not want ranges
on opposite strands to be considered overlapping.

To generally compute on the overlaps, we call `findOverlaps()` to
return a `Hits` object, which is essentially a bipartite graph
matching query ranges to overlapping subject ranges.


```r
hits <- findOverlaps(query, subject, ignore.strand=TRUE)
```

We typically use the hits to perform one of two operations: join and
aggregate. For example, we could inner join the scores from the
subject using the query and subject indexes,

```r
joined <- query[queryHits(hits)]
joined$score <- subject$score[subjectHits(hits)]
```
The above carries over a single metadata column from the
subject. Similar code would carry over other columns and even the
ranges themselves. 

Sometimes, we want to merge the matched query and subject ranges,
typically by finding their intersection,

```r
ranges(joined) <- ranges(pintersect(joined, subject[subjectHits(hits)]))
```

The typical aggregation is counting the number of hits overlapping a
query. In general, aggregation starts by grouping the subject hits by
query hits, which we express as a coercion to a *List*,

```r
hitsByQuery <- as(hits, "List")
```
The result is an *IntegerList*, a type of *AtomicList*. *AtomicList*
objects have many methods for efficient aggregation. In this case, we
just call `lengths()` to get the count:

```r
counts <- lengths(hitsByQuery)
```
Since this a common operation, there are shortcuts,

```r
counts <- countQueryHits(hits)
```
or even shorter and more efficient,

```r
counts <- countOverlaps(query, subject, ignore.strand=TRUE)
unname(counts)
#>  [1] 0 0 0 2 0 0 0 0 0 2 0 2 0 0 2 0 0 2 2 0 0 0 1 0 0 0 2 0 1 0
```

Often, we want to combine joins and aggregations. For example, we may
want to annotate each query with the maximum score among the subject
hits,

```r
query$maxScore <- max(extractList(subject$score, hitsByQuery))
subset(query, maxScore > 0)
#> GRanges object with 9 ranges and 2 metadata columns:
#>      seqnames    ranges strand | replicate  maxScore
#>         <Rle> <IRanges>  <Rle> |     <Rle> <integer>
#>   10     chr2   120-129      * |         1         3
#>   28     chr1   151-160      * |         1         6
#>    5     chr2   128-137      * |         2         4
#>   14     chr1   149-158      * |         2         6
#>   23     chr3   153-162      * |         2         8
#>   26     chr3   160-169      * |         2         9
#>    9     chr1    92-101      * |         3         1
#>   21     chr1   148-157      * |         3         6
#>   27     chr3   190-199      * |         3        10
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

In rare cases, we can more or less arbitrarily select one of the
subject hits. The `select=` argument to `findOverlaps()` automatically
selects an "arbitrary", "first" (in subject order) or "last" subject
range,

```r
hits <- findOverlaps(query, subject, select="first", ignore.strand=TRUE)
hits <- findOverlaps(query, subject, select="arbitrary", ignore.strand=TRUE)
hits
#>  [1] NA NA NA  2 NA NA NA NA NA  5 NA  3 NA NA  5 NA NA  7  8 NA NA NA  1
#> [24] NA NA NA  5 NA 10 NA
```

## Exercises

1. Find the average intensity of the X and Y measurements for each
   each replicate over all positions in the query object
2. Add a new column to the intensities object that is the distance from
   each position to its closest gene (hint `IRanges::distance()`)
3. Find flanking regions downstream of the genes in gr that have width of 8bp
4. Are any of the intensities positions within the flanking region?

## Example: exploring BigWig files from AnnotationHub

In the workflow of ChIP-seq data analysis, we are often interested in
finding peaks from islands of coverage over a chromosome. Here we will
use plyranges to explore ChiP-seq data from the Human Epigenome
Roadmap project @Roadmap-Epigenomics-Consortium2015-pr.

### Extracting data from AnnotationHub

This data is available on Bioconductor's AnnotationHub. First we construct
an AnnotationHub, and then `query()` for all bigWigFiles related to 
the project that correspond to the following conditions:

1. are from methylation marks (H3K4ME in the title)
2. correspond to primary T CD8+ memory cells from peripheral blood
3. correspond to unimputed log10 P-values 

First we construct a hub that contains all references to the EpigenomeRoadMap
data and extract the metadata as a data.frame:


```r
library(AnnotationHub)
ah <- AnnotationHub()
#> snapshotDate(): 2018-06-27
roadmap_hub <- query(ah, "EpigenomeRoadMap") 
metadata <- query(ah, "Metadata")[[1L]]
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.AnnotationHub/47270'
head(metadata)
#>    EID    GROUP   COLOR          MNEMONIC
#> 1 E001      ESC #924965            ESC.I3
#> 2 E002      ESC #924965           ESC.WA7
#> 3 E003      ESC #924965            ESC.H1
#> 4 E004 ES-deriv #4178AE ESDR.H1.BMP4.MESO
#> 5 E005 ES-deriv #4178AE ESDR.H1.BMP4.TROP
#> 6 E006 ES-deriv #4178AE       ESDR.H1.MSC
#>                                     STD_NAME
#> 1                                ES-I3 Cells
#> 2                               ES-WA7 Cells
#> 3                                   H1 Cells
#> 4 H1 BMP4 Derived Mesendoderm Cultured Cells
#> 5 H1 BMP4 Derived Trophoblast Cultured Cells
#> 6          H1 Derived Mesenchymal Stem Cells
#>                                   EDACC_NAME     ANATOMY           TYPE
#> 1                            ES-I3_Cell_Line         ESC PrimaryCulture
#> 2                           ES-WA7_Cell_Line         ESC PrimaryCulture
#> 3                               H1_Cell_Line         ESC PrimaryCulture
#> 4 H1_BMP4_Derived_Mesendoderm_Cultured_Cells ESC_DERIVED     ESCDerived
#> 5 H1_BMP4_Derived_Trophoblast_Cultured_Cells ESC_DERIVED     ESCDerived
#> 6          H1_Derived_Mesenchymal_Stem_Cells ESC_DERIVED     ESCDerived
#>   AGE    SEX SOLID_LIQUID ETHNICITY SINGLEDONOR_COMPOSITE
#> 1  CL Female         <NA>      <NA>                    SD
#> 2  CL Female         <NA>      <NA>                    SD
#> 3  CL   Male         <NA>      <NA>                    SD
#> 4  CL   Male         <NA>      <NA>                    SD
#> 5  CL   Male         <NA>      <NA>                    SD
#> 6  CL   Male         <NA>      <NA>                    SD
```

To find out the name of the sample corresponding to 
primary memory T-cells we can filter the data.frame. We extract the 
sample ID corresponding to our filter.


```r
primary_tcells <- subset(metadata,
                         ANATOMY == "BLOOD" & TYPE == "PrimaryCell" &
                             EDACC_NAME == "CD8_Memory_Primary_Cells")$EID
primary_tcells <- as.character(primary_tcells)
```

Now we can take our roadmap hub and query it based on our other conditions:


```r
methylation_files <-  query(roadmap_hub,
                            c("BigWig", primary_tcells, "H3K4ME[1-3]",
                              "pval.signal"))
methylation_files
#> AnnotationHub with 5 records
#> # snapshotDate(): 2018-06-27 
#> # $dataprovider: BroadInstitute
#> # $species: Homo sapiens
#> # $rdataclass: BigWigFile
#> # additional mcols(): taxonomyid, genome, description,
#> #   coordinate_1_based, maintainer, rdatadateadded, preparerclass,
#> #   tags, rdatapath, sourceurl, sourcetype 
#> # retrieve records with, e.g., 'object[["AH33454"]]' 
#> 
#>             title                                  
#>   AH33454 | E048-H3K4me1.pval.signal.bigwig        
#>   AH33455 | E048-H3K4me3.pval.signal.bigwig        
#>   AH39974 | E048-H3K4me1.imputed.pval.signal.bigwig
#>   AH40101 | E048-H3K4me2.imputed.pval.signal.bigwig
#>   AH40228 | E048-H3K4me3.imputed.pval.signal.bigwig
```

So we'll take the first two entries and download them as BigWigFiles:


```r
bw_files <- lapply(methylation_files[1:2], `[[`, 1L)
#> require("rtracklayer")
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.AnnotationHub/38894'
#> downloading 0 resources
#> loading from cache 
#>     '/home/mramos//.AnnotationHub/38895'
```

We have our desired BigWig files so now we can we can start analyzing them.

### Reading BigWig files

For this analysis, we will call peaks from a score vector over
chromosome 10.

First, we extract the genome information from the first BigWig file and filter
to get the range for chromosome 10. This range will be used as a filter when 
reading the file.


```r
chr10_ranges <- Seqinfo(genome="hg19")["chr10"]
```

Then we read the BigWig file only extracting scores if they overlap chromosome
10.


```r
library(rtracklayer)
chr10_scores <- lapply(bw_files, import, which = chr10_ranges,
                       as = "RleList") 
chr10_scores[[1]]$chr10
#> numeric-Rle of length 135534747 with 5641879 runs
#>   Lengths:               60612                 172 ...                9907
#>   Values :  0.0394200012087822   0.154219999909401 ...                   0
```
Each of element of the list is a run-length encoded vector of the
scores for a particular signal type.

We find the islands by slicing the vectors,

```r
islands <- lapply(chr10_scores, slice, lower=1L)
```
where the islands are represented as *Views* objects, i.e., ranges of
interest over a genomic vector. Then we find the summit within each
island,

```r
summits <- lapply(islands, viewRangeMaxs)
```
using the optimized `viewRangeMaxs()` function. Each element of the
`summits` list is a *RangesList* object, holding the ranges for each
summit. The structure of the *RangesList* keeps track of the
chromosome (10) of the summits (there could be multiple chromosomes in
general). We broaden the summits and reduce them in order to smooth the
peak calls and provide some context,

```r
summits <- lapply(lapply(summits, `+`, 50L), reduce)
```

After this preprocessing, we want to convert the result to a more
familiar and convenient GRanges object containing an *RleList* "score"
column containing the score vector for each summit,

```r
summits_grs <- lapply(summits, GRanges)
score_grs <- mapply(function(scores, summits) {
    summits$score <- scores[summits]
    seqlengths(summits) <- lengths(scores)
    summits
}, chr10_scores, summits_grs)
score_gr <- stack(GenomicRangesList(score_grs), index.var="signal_type")
```
One problem with *RangesList* is that it does not keep track of the
sequence lengths, so we need to add those after forming the *GRanges*.

We could then find summits with the maximum summit height within each
signal type: 

```r
score_gr$score_max <- max(score_gr$score)
chr10_max_score_region <- aggregate(score_gr, score_max ~ signal_type, max)
```

### Exercises

1. Use the `reduce_ranges()` function to find all peaks for each signal type.
2. How could you annotate the scores to find out which genes overlap
each peak found in 1.?
3. Plot a 1000nt window centred around the maximum scores for each signal type 
using the `ggbio` or `Gviz` package.

## Worked example: coverage analysis of BAM files

A common quality control check in a genomics workflow is to perform
coverage analysis over features of interest or over the entire
genome. Here we use the data from the airway package to operate on
read alignment data and compute coverage histograms.

First let's gather all the BAM files available to use in airway (see
`browseVignettes("airway")` for more information about the data and how it 
was prepared):


```r
library(tools)
bams <- list_files_with_exts(system.file("extdata", package = "airway"), "bam")
names(bams) <- sub("_[^_]+$", "", basename(bams))
library(Rsamtools)
#> Loading required package: Biostrings
#> Loading required package: XVector
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
bams <- BamFileList(bams)
```
Casting the vector of filenames to a formal *BamFileList* is critical
for informing the following code about the nature of the files.

To start let's look at a single BAM file (containing only reads from
chr1). We can compute the coverage of the alignments over all contigs
in the BAM as follows:


```r
first_bam <- bams[[1L]]
first_bam_cvg <- coverage(first_bam)
```

The result is a list of *Rle* objects, one per chromosome. Like other
*AtomicList* objects, we call pass our *RleList* to `table()` to
compute the coverage histogram by chromosome,

```r
head(table(first_bam_cvg)[1L,])
#>         0         1         2         3         4         5 
#> 249202844     15607      5247      3055      2030      1280
```

For RNA-seq experiments we are often interested in splitting up
alignments based on whether the alignment has skipped a region from
the reference (that is, there is an "N" in the cigar string,
indicating an intron). We can represent the nested structure using a
*GRangesList* object.

To begin we read the BAM file into a *GAlignments* object using
`readGAlignments()` and extract the ranges, chopping by introns, using
`grglist()`,

```r
library(GenomicAlignments)
#> Loading required package: SummarizedExperiment
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:AnnotationHub':
#> 
#>     cache
#> Loading required package: DelayedArray
#> Loading required package: matrixStats
#> 
#> Attaching package: 'matrixStats'
#> The following objects are masked from 'package:Biobase':
#> 
#>     anyMissing, rowMedians
#> Loading required package: BiocParallel
#> 
#> Attaching package: 'DelayedArray'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
#> The following object is masked from 'package:Biostrings':
#> 
#>     type
#> The following objects are masked from 'package:base':
#> 
#>     aperm, apply
reads <- grglist(readGAlignments(first_bam))
```
Finally, we can find the junction reads:

```r
reads[lengths(reads) >= 2L]
#> GRangesList object of length 3833:
#> [[1]] 
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames            ranges strand
#>          <Rle>         <IRanges>  <Rle>
#>   [1]        1 11072744-11072800      +
#>   [2]        1 11073773-11073778      +
#> 
#> [[2]] 
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames            ranges strand
#>   [1]        1 11072745-11072800      -
#>   [2]        1 11073773-11073779      -
#> 
#> [[3]] 
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames            ranges strand
#>   [1]        1 11072746-11072800      +
#>   [2]        1 11073773-11073780      +
#> 
#> ...
#> <3830 more elements>
#> -------
#> seqinfo: 84 sequences from an unspecified genome
```

We typically want to count how many reads overlap each gene. First, we
get the transcript structures as a *GRangesList* from Ensembl,

```r
library(GenomicFeatures)
#> Loading required package: AnnotationDbi
library(EnsDb.Hsapiens.v75)
#> Loading required package: ensembldb
#> Loading required package: AnnotationFilter
#> 
#> Attaching package: 'ensembldb'
#> The following object is masked from 'package:stats':
#> 
#>     filter
tx <- exonsBy(EnsDb.Hsapiens.v75, "gene")
```

Finally, we count how many reads overlap each transcript,

```r
reads <- keepStandardChromosomes(reads)
counts <- countOverlaps(tx, reads, ignore.strand=TRUE)
head(counts[counts > 0])
#> ENSG00000009724 ENSG00000116649 ENSG00000120942 ENSG00000120948 
#>              89            2006             436            5495 
#> ENSG00000171819 ENSG00000171824 
#>               8            2368
```

To do this over every sample, we use the `summarizeOverlaps()`
convenience function,

```r
airway <- summarizeOverlaps(features=tx, reads=bams,
                            mode="Union", singleEnd=FALSE,
                            ignore.strand=TRUE, fragments=TRUE)
airway
#> class: RangedSummarizedExperiment 
#> dim: 64102 8 
#> metadata(0):
#> assays(1): counts
#> rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
#> rowData names(0):
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(0):
```
The `airway` object is a *SummarizedExperiment* object, the central
Bioconductor data structure for storing results summarized per feature
and sample, along with sample and feature metadata. It is at the point
of summarization that workflows switch focus from the ranges
infrastructure to Bioconductor modeling packages, most of which
consume the *SummarizedExperiment* data structure, so this is an
appropriate point to end this tutorial.

### Exercises 

1. Compute the total depth of coverage across all features.
1. How could you compute the proportion of bases covered over an entire genome?
   (hint: `seqinfo` and `S4Vectors::merge`)
1. How could you compute the strand specific genome wide coverage?
1. Create a workflow for computing the strand specific coverage for
   all BAM files.
1. For each sample plot total breadth of coverage against the number of bases
   covered faceted by each sample name.

## Conclusions

The Bioconductor ranges infrastructure is rich and complex, and it can
be intimidating to new users. However, the effort invested will pay
dividends, especially when generalizing a series of bespoke analyses
into a reusable contribution to Bioconductor.
