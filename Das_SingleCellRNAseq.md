
# Analysis of single-cell RNA-seq data: Dimensionality reduction, clustering, and lineage inference

## Instructor(s) name(s) and contact information
* Diya Das (@diyadas, <diyadas@berkeley.edu>)
* Davide Risso (@drisso)
* Kelly Street (@kstreet13)

## Workshop Description

This workshop will be presented as a lab session (brief introduction followed by hands-on coding)
that instructs participants in a Bioconductor workflow for the analysis of single-cell RNA-sequencing data, in three parts:
1. dimensionality reduction that accounts for zero inflation, over-dispersion, and batch effects
2. cell clustering that employs a resampling-based approach resulting in robust and stable clusters
3. lineage trajectory analysis that uncovers continuous, branching developmental processes

We will provide worked examples for lab sessions, and a set of stand-alone notes in this repository.

Note to organizers: A previous version of this workshop was well-attended at BioC 2017,
but the tools presented have been significantly updated for
interoperability (most notably, through the use of the `SingleCellExperiment` class), and we have been receiving many requests to provide an
updated workflow. We plan to incorporate feedback from this workshop into a revised version of our F1000 Workflow.

### Pre-requisites

We expect basic knowledge of R syntax. Some familiarity with S4 objects may be helpful, but not required.
More importantly, participants should be familiar with the concept and design of RNA-sequencing experiments. Direct experience with single-cell RNA-seq is not required, and the main challenges of single-cell RNA-seq compared to bulk RNA-seq will be illustrated.

### Workshop Participation

This will be a hands-on workshop, in which each student, using their laptop, will analyze a provided example datasets. The workshop will be a mix of example code that the instructors will show to the students (available through this repository) and short exercises.

### _R_ / _Bioconductor_ packages used

1. _zinbwave_ : https://bioconductor.org/packages/zinbwave
2. _clusterExperiment_: https://bioconductor.org/packages/clusterExperiment
3. _slingshot_: https://github.com/kstreet13/slingshot (will be submitted to Bioconductor shortly)

### Time outline

2 hr workshop:

| Activity                                   | Time |
|--------------------------------------------|------|
| Intro to single-cell RNA-seq analysis      | 15m  |
| zinbwave (dimensionality reduction)        | 30m  |
| clusterExperiment (clustering)             | 30m  |
| slingshot (lineage trajectory analysis)    | 30m  |
| Questions / extensions                     | 15m  |

## Workshop goals and objectives

### Learning goals

* describe the goals of single-cell RNA-seq analysis 
* identify the main steps of a typical single-cell RNA-seq analysis
* evaluate the results of each step in terms of model fit 
* synthesize results of successive steps to interpret biological significance and develop biological models
* apply this workflow to carry out a complete analysis of other single-cell RNA-seq datasets

### Learning objectives

* compute and interpret low-dimensional representations of single-cell data
* identify and remove sources of technical variation from the data
* identify sub-populations of cells (clusters) and evaluate their robustness
* infer lineage trajectories corresponding to differentiating cells
* order cells by developmental "pseudotime"
* identify genes that play an important role in cell differentiation 
