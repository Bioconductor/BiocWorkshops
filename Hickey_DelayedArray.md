
# Effectively using the DelayedArray framework to support the analysis of large  datasets

## Instructor(s) name(s) and contact information
* Peter Hickey (GitHub/Twitter: @PeteHaitch, email: peter.hickey@gmail.com)

## Workshop Description

This workshop will teach the fundamental concepts underlying the DelayedArray framework and related infrastructure. 
It is intended for package developers who want to learn how to use the DelayedArray framework to support the analysis of large datasets, particularly through the use of on-disk data storage.
The first part of the workshop will provide an overview of the DelayedArray 
infrastructure and introduce computing on DelayedArray objects using delayed operations and block-processing.
The second part of the workshop will present strategies for adding support for DelayedArray to an existing package and extending the DelayedArray framework.
Students can expect a mixture of lecture and question-and-answer session to teach the fundamental concepts.
There will be plenty of examples to illustrate common design patterns for writing performant code, although we will not be writing much code during the workshop.

### Pre-requisites

* Solid understanding of R
* Familiarity with common operations on arrays (e.g., `colSums()` and those available in the matrixStats package)
* Familiarity with object oriented programming, particularly S4, will be useful but is not essential
* No familiarity required of technical details of particular data storage backends (e.g., HDF5, sparse matrices)
* No familiarity required of particular biological techniques (e.g., single-cell RNA-seq)

### Workshop Participation

Questions and discussion are encouraged! 
This will be especially important to guide the second half of the worshop which focuses on integrating DelayedArray into an existing or new Bioconductor package.
Students will be expected to be able to follow and reason about example R code. 

### _R_ / _Bioconductor_ packages used

* DelayedArray
* HDF5Array
* SummarizedExperiment
* DelayedMatrixStats
* beachmat

### Time outline

| Activity                                          | Time |
|---------------------------------------------------|------|
| Introductory slides                               | 15m  |
| Part 1: Overview of DelayedArray framework        | 45m  |
| Part 2: Incoporating DelayedArray into a package  | 45m  |
| Questions and discussion                          | 15m  |

## Workshop goals and objectives

### Learning goals

* Identify when it is useful to use a DelayedArray instead of an ordinary array or other array-like data structure.
* Become familiar with the fundamental concepts of delayed operations, block-processing, and realization.
* Learn of existing functions and packages for constructing and computing on DelayedArray objects, avoiding the need to re-invent the wheel.
* Learn common design patterns for writing performant code that operates on a DelayedArray.
* Evaluate whether an existing function that operates on an ordinary array can be readily adapted to work on a DelayedArray.
* Reason about potential bottlenecks in algorithms operating on DelayedArray objects.

### Learning objectives

* Understand the differences between a DelayedArray instance and an instance of a subclass (e.g., HDF5Array, RleArray).
* Know what types of operations 'degrade' an instance of a DelayedArray subclass to a DelayedArray, as well as when and why this matters.
* Construct a DelayedArray:
  * From an in-memory array-like object.
  * From an on-disk data store (e.g., HDF5).
  * From scratch by writing data to a RealizationSink.
* Take a function that operates on rows or columns of a matrix and apply it to a DelayedMatrix.
* Use block-processing on a DelayedArray to compute:
  * A univariate (scalar) summary statistic (e.g., `max()`).
  * A multivariate (vector) summary statistic (e.g., `colSum()` or `rowMean()`).
  * A multivariate (array-like) summary statistic (e.g., `rowRanks()`).
* Design an algorithm that imports data into a DelayedArray.
