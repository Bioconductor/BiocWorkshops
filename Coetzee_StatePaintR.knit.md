# Variant Functional Annotation using StatePaintR, FunciVar & MotifBreakR

## Instructor(s) name(s) and contact information
* Simon G. Coetzee (scoetzee@gmail.com)
* Dennis J. Hazelett (dennis.hazelett@csmc.edu)

## Workshop Description
This workshop will entail lecture and live demo of StateHub/
StatePaintR and funciVar bioconductor packages.

### Pre-requisites
It is recommended that workshop participants have:
* Basic concepts in epigenomics and GWAS
* Intermediate R skills, familiarity with GRanges

Relevant background reading for the workshop.

* [StateHub/StatePaintR][1]
* [BioC2017 Workshop][2]
* [MotifBreakR][3]

### Workshop Participation
Describe how students will be expected to participate in the workshop.

### _R_ / _Bioconductor_ packages used

* StatePaintR
* FunciVar
* MotifBreakR

### Time outline
1 hour Workshop:

| Activity                     | Time |
|------------------------------|------|
| **StatePaintR**              | 30m  |
| * intro and theory           | (10m) |
| * Generate Decision Matrix   | (10m) |
| * live demo                  | (10m) |
| **FunciVar**                 | 15m  |
| * intro                      | (10m) |
| * live demo                  | (5m) |
| **MotifBreakR**              | 15m  |
| * intro & demo               | (10m) |
| * Shiny webapp demo          | (5m) |

## Workshop goals and objectives
This workshop focuses on bioconductor based tools for non-coding 
variant annotation. I will present a high level overview of these concepts and the bioconductor tools our group uses to address them *in silico*. These tools can be used for germline or somatic variants, but as we will see, can also be used for any other type of feature represented as GRanges objects.

### Key Concepts
I will provide brief overview with examples of:

* Genome segmentation and Chromatin State
* Feature enrichment analysis
* Transcription factors and Motif Disruption

### Learning objectives
At the end of this workshop students will understand the basic inputs and outputs of all these analyses, as well as an intuitive understanding of the tools and where to find them:

* [The StateHub website][4]
* StatePaintR genome segmentation software [Bioconductor][5]
* FunciVar for variant annotation and enrichment
* MotifBreakR Transcription Factor motif disruption [Bioconductor][6]

[1]: https://f1000research.com/articles/7-214/v1
[2]: https://www.simoncoetzee.com/bioc2017.html
[3]: https://bioconductor.org/packages/release/bioc/html/motifbreakR.html
[4]: http://statehub.org
[5]: https://github.com/Simon-Coetzee/StatePaintR/blob/master/vignettes/StatePaintR_Analysis.Rmd
[6]: https://www.bioconductor.org/packages/devel/bioc/vignettes/motifbreakR/inst/doc/motifbreakR-vignette.html
