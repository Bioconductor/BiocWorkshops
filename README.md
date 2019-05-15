---
knit: "bookdown::render_book"
title: "The Bioconductor 2018 Workshop Compilation"
description: "This book is a central repository for all the workshops submitted to the Bioconductor 2018 Conference"
site: bookdown::bookdown_site
github-repo: Bioconductor/BiocWorkshops
documentclass: book
---

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

The most up-to-date error messages and `.md` files from builds of individual workshops are [here](https://github.com/Bioconductor/BiocWorkshops/tree/individual_builds).

# Introduction

Author:
    Martin Morgan^[Roswell Park Comprehensive Cancer Center, Buffalo, NY].
    <br/>
Last modified: 22 May, 2018.

The most recently compiled version of the book is available at https://bioconductor.github.io/BiocWorkshops/.

## For Everyone

This book contains workshops used in _R_ / _Bioconductor_
training. The workshops are divided into 3 sections:

- **Learn** (100-series chapters) contains material for beginning
  users of _R_ and _Bioconductor_. The _Bioconductor_-related material
  is relevant even for experienced _R_ users who are new to
  _Bioconductor_.

- **Use** (200-series chapters) contains workshops emphasizing use of
  _Bioconductor_ for common tasks, e.g., bulk RNA-seq differential
  expression, ChIP-seq, single-cell analysis, gene set enrichment, and
  annotation.

- **Develop** (500-series chapters) contains workshops to help expert
  users hone their skills and contribute their domain-specific
  knowledge to the _Bioconductor_ community.

## For Workshop Authors

To contribute a new workshop, open a [BiocWorkshops issue][] asking to
be added as a collaborator.

### DESCRIPTION

Update the DESCRIPTION file adding packages utilized in your workshop to
the **Imports** field. Please be kind and don't remove anyone elses packages from
the DESCRIPTION as this is a shared file for all workshops. Also don't add your packages
to the **Depends** field.

### Classifying your workshop

Follow the numbering scheme above to classify your workshop by preceding your title with 
"Bioconductor 1xx:", "Bioconductor 2xx:", or "Bioconductor 5xx:". Final workshop numbers will be 
determined by an editor.

### Bookdown

Write your workshop as a stand-alone markdown document, using your workshop 
syllabus already posted here as a starting point.  We are using [bookdown][] 
('knit-then-merge' strategy) to compile the workshops and produce a book that will be 
available through Continuous Integration at https://bioconductor.github.io/BiocWorkshops/.

You do not need to build the complete book yourself, it is sufficient to ensure that your own
workshop compiles successfully. You can 1) compile draft versions using a standard "knit" 
procedure to compile your chapter, and 2) follow these bookdown procedures to see how your 
chapter will look in the finished product. Install the [bookdown package][] package from CRAN, 
as well as [pandoc][]. Render your chapter with the `preview=` option to `render_book()`, e.g.,

```
Rscript -e "bookdown::render_book(                             \
    'xxx_Your_Workshop.Rmd', 'bookdown::gitbook', preview=TRUE \
)"
```

As this is a shared space for all workshop contributors, in order to
use the above command in the BiocWorkshops directory, the index has to
be built at least once, which can be time consuming depending on how many
workshops have already been submitted.

```
Rscript -e "bookdown::render_book(                             \
    'index.Rmd', 'bookdown::gitbook')"
```

To avoid having to build all workshops but still be able to preview
your individual workshop we recommend creating a soft link to your .Rmd file.
We recommend having the file in the `BiocWorkshop/` and the soft link in
any other directory on your system. By default, this will generate an
html file in `_book/` wherever this command is run.

```
mkdir tmp
cd tmp/
ln -s ../xxx_Your_Workshop.Rmd
Rscript -e "bookdown::render_book(                             \
    'xxx_Your_Workshop.Rmd', 'bookdown::gitbook', preview=TRUE \
)"
```


Push **only** your .Rmd file to the BiocWorkshop repository; the book will be
rebuilt manually or automatically. Eventually the output will be
available for end-users at https://bioconductor.github.io/BiocWorkshops .The
master branch will not contain the built version of the book. Switching to the
[gh-pages branch][] will show built output.  

## Deadlines for Bioc2018

Please be aware of the following deadlines for the [Bioconductor 2018 Conference][] in Toronto

- **Fri June 29:** draft workshop materials submitted to this Bioconductor GitHub bookdown site

- **Fri July 6:** feedback period completes

- **Weds July 18:** workshops must pass checks without errors or warnings (All materials will be checked by Continuous Integration)

- **Thurs / Fri July 26-27:** Bioc2018

[BiocWorkshops issue]: https://github.com/Bioconductor/BiocWorkshops/issues
[bookdown]: https://bookdown.org/yihui/bookdown/
[bookdown package]: https://cran.r-project.org/package=bookdown
[pandoc]: http://pandoc.org/
[gh-pages branch]: https://github.com/Bioconductor/BiocWorkshops/tree/gh-pages
[Bioconductor 2018 Conference]: https://bioc2018.bioconductor.org/
