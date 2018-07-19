# Working with Genomic Data in R with the DECIPHER package

Authors:
    Nicholas Cooley^[University of Pittsburgh]
Last Modified: 18 July, 2018
    
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
* stringr
* FindHomology

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

#### Workshop

In this worshop we will be walking through a comparative genomics pipeline using functions within the R package DECIPHER, and within a package that is currently under construction under the name FindHomology.

Lipid membranes are ubiquitious in living organisms, however the individual lipid molecules that compose membranes have a wide variety of structures. The biosynthesis of common lipid molecules is well understood, conversely the biosynthesis of more exotic lipid molecules is not. The genomes selected today are from the archaeal phylum Thaumarchaeota, which are known for being anammox archaea, and are important in the nitrogen cycle. In addition to that, these archaea all produce a series of unusual lipids. Since these lipids are shared across the entire phylum from which these genomes are pulled, it is possible that the biosynthetic genes responsible for lipid construction are shared across the selected genomes, and part of their 'core' genome.

Archaea are relatively newly characterized as organisms. First classified in the '70s[^2], for their deviations from Prokaryotes, many Archaea are noteable for being isolated from extreme environments. Despite the colorful beginnings of their classification, not all Archaea are extremophiles.

[^2]: [(1) Woese, C. R.; Fox, G. E. Proc. Natl. Acad. Sci. 1977, 74 (11), 5088–5090.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC432104/)

<div class="figure" style="text-align: center">
<img src="./DECIPHERWorkshop/GrandPrismatic.png" alt="Grand Prismatic Spring in Yellowstone National Park[^3]" width="75%" height="75%" />
<p class="caption">(\#fig:Grand Prismatic)Grand Prismatic Spring in Yellowstone National Park[^3]</p>
</div>

The aforementioned unusual lipds[^4] are striking in their differences from those common in both Eukaryotes and Prokaryotes.

[^3]: [Grand Prismatic Spring](https://commons.wikimedia.org/wiki/File:Grand_prismatic_spring.jpg)

[^4]: [(1) Caforio, A.; Driessen, A. J. M. Biochimica et Biophysica Acta - Molecular and Cell Biology of Lipids. Elsevier November 1, 2017, pp 1325–1339.](https://www.ncbi.nlm.nih.gov/pubmed/28007654)

<div class="figure" style="text-align: center">
<img src="./DECIPHERWorkshop/ArchaealLipids.png" alt="A. A generic lipid with phosphatidyl serine head group B. A generic Archaeal lipid C. A more exotic Archaeal lipid" width="50%" height="50%" />
<p class="caption">(\#fig:Lipids)A. A generic lipid with phosphatidyl serine head group B. A generic Archaeal lipid C. A more exotic Archaeal lipid</p>
</div>

Whereas Lipid construction in eucharyotes and prokaryotes relies on iterative incorporation of malonyl-CoA units, archaea utilize the melvonate pathway to construct their lipids through iterative incorporation of dimethylallyl pyrophosphate. Which is responsible for the iterative methyl groups extant to the lipid chain in **B** above. Archaeal lipids are also notable for their lack of ester linkages between the alkyl chain, and glyerol. Even more interesting, and as yet uncharacterized, Archaea incorporate a range of ring structures into their lipids, and additionally appear to *fuse lipid tails together in the bilayer*, as shown above in **C**. 

#### The above is interesting, but so what?

Our goal with this workshop is to show how functions within DECIPHER, in addition to those found within FindHomology can be used to predict homologs between genomes, as well as predict the core genome of a given set of genomes. As stated above, the chemistry and biology involved in some of the construction of archaeal lipids is unknown. Hypothetically, as this chemistry is universal in these organisms, the genes responsible could all be homologs, and be identified as part of the core genome. Further we can compare the core and pan genomes of this set.

To begin with, we will load the required packages:
DECIPHER and FindHomology will provide tools to complete the phylogenetic analysis, while phytools will provide visualization tools, specifically for the visualization of tangelograms, and stringr is used for a specific string operation.

#### Packages and Sequences


```r
suppressMessages(library(DECIPHER))
library(FindHomology)
library(stringr)
```

#### Genomes, and DECIPHER databases

The first step in this process is selecting the genomes we want to work with. There are several ways to do this, especially when accessing publically available genomes through the NCBI. Entrez Direct[^5] provides powerful unix command line tools, if that's your flavor, but the NCBI Genome List, and other search tools are easily accessible, and have a lower barrier to to entry.

[^5]: [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

For the purposes of this workshop, we will be using a set of preselected genomes that will be accessed via ftp addresses. These are loaded with the FindHomology Package and can be accessed easily.


```r
data("GeneCallAdds")
```

```
## Warning in data("GeneCallAdds"): data set 'GeneCallAdds' not found
```

```r
data("GenomeAdds")
```

```
## Warning in data("GenomeAdds"): data set 'GenomeAdds' not found
```

```r
data("GenomeIDs")
```

```
## Warning in data("GenomeIDs"): data set 'GenomeIDs' not found
```

We have now loaded in a character vector for ftp addresses for genomic fasta files `GenomeAdds`, a character vector of ftp addresses for genomic GFF files `GeneCallAdds` that we will be using to collect gene boundaries, strandedness, and annotations. And lastly a vector of abbreviated names for the genomes we are working with.





















































