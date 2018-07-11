#!/bin/bash                                                                                                                                   

Rscript -e "rmarkdown::render('201_Love_DESeq2.Rmd', clean=FALSE, params = list(cache=TRUE))" > 201_Love_DESeq2.log 2> 201_Love_DESeq2.err &
Rscript -e "rmarkdown::render('220_Das_SingleCellRNASeq.Rmd', clean=FALSE, params = list(cache=TRUE))" > 220_Das_SingleCellRNASeq.log 2> 220_Das_SingleCellRNASeq.err &
Rscript -e "rmarkdown::render('Cooley_DECIPHER.Rmd', clean=FALSE, params = list(cache=TRUE))" > Cooley_DECIPHER.log 2> Cooley_DECIPHER.err &
Rscript -e "rmarkdown::render('230_Isserlin_RCy3_intro.Rmd', clean=FALSE, params = list(cache=TRUE))" > 230_Isserlin_RCy3_intro.log 2> 230_Isserlin_RCy3_intro.err &
Rscript -e "rmarkdown::render('500_Effectively_Using_the_DelayedArray_Framework.Rmd', clean=FALSE, params = list(cache=TRUE))" > 500_Effectively_Using_the_DelayedArray_Framework.log 2> 500_Effectively_Using_the_DelayedArray_Framework.err &
Rscript -e "rmarkdown::render('Geistlinger_Omics.Rmd', clean=FALSE, params = list(cache=TRUE))" > Geistlinger_Omics.log 2> Geistlinger_Omics.err &
Rscript -e "rmarkdown::render('Lawrence_GenomicRanges.Rmd', clean=FALSE, params = list(cache=TRUE))" > Lawrence_GenomicRanges.log 2> Lawrence_GenomicRanges.err &
Rscript -e "rmarkdown::render('Law_RNAseq123.Rmd', clean=FALSE, params = list(cache=TRUE))" > Law_RNAseq123.log 2> Law_RNAseq123.err &
Rscript -e "rmarkdown::render('Lee_Plyranges.Rmd', clean=FALSE, params = list(cache=TRUE))" > Lee_Plyranges.log 2> Lee_Plyranges.err &
Rscript -e "rmarkdown::render('MacDonald_Annotation.Rmd', clean=FALSE, params = list(cache=TRUE))" > MacDonald_Annotation.log 2> MacDonald_Annotation.err &
Rscript -e "rmarkdown::render('Ramos_MultiAssayExperiment.Rmd', clean=FALSE, params = list(cache=TRUE))" > Ramos_MultiAssayExperiment.log 2> Ramos_MultiAssayExperiment.err &
Rscript -e "rmarkdown::render('Safikhani_Pharmacogenomics.Rmd', clean=FALSE, params = list(cache=TRUE))" > Safikhani_Pharmacogenomics.log 2> Safikhani_Pharmacogenomics.err &
Rscript -e "rmarkdown::render('Turaga_MaintainBioc.Rmd', clean=FALSE, params = list(cache=TRUE))" > Turaga_MaintainBioc.log 2> Turaga_MaintainBioc.err &
Rscript -e "rmarkdown::render('Waldron_PublicData.Rmd', clean=FALSE, params = list(cache=TRUE))" > Waldron_PublicData.log 2> Waldron_PublicData.err &
