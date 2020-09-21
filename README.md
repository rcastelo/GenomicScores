# GenomicScores: seamless access to genomewide position-specific scores from R and Bioconductor

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/GenomicScores.svg)](https://bioconductor.org/packages/release/bioc/html/GenomicScores.html "How long has been GenomicScores in a release of Bioconductor")
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/release/GenomicScores.svg)](https://bioconductor.org/packages/stats/bioc/GenomicScores.html "Ranking by number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months")
[![Support posts](https://bioconductor.org/shields/posts/GenomicScores.svg)](https://support.bioconductor.org/t/GenomicScores/ "Support site activity on GenomicScores, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")
<img align="right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/GenomicScores/GenomicScores.png" height="200"/>

**Current build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/release/GenomicScores.svg)](https://bioconductor.org/packages/release/bioc/html/GenomicScores.html#archives "Whether GenomicScores release is available on all platforms") 
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/release/GenomicScores.svg)](https://bioconductor.org/packages/release/bioc/html/GenomicScores.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/release/bioc/GenomicScores.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/GenomicScores "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/GenomicScores.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/GenomicScores/ "Bioconductor release build")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/devel/GenomicScores.svg)](https://bioconductor.org/packages/devel/bioc/html/GenomicScores.html#archives "Whether GenomicScores devel is available on all platforms") 
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/devel/GenomicScores.svg)](https://bioconductor.org/packages/devel/bioc/html/GenomicScores.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/devel/bioc/GenomicScores.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/GenomicScores "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/GenomicScores.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/GenomicScores/ "Bioconductor devel build")

The `GenomicScores` package facilitates an efficient storage and seamless access of genomic scores, and their integration into genome analysis workflows on top of R and Bioconductor. Users using these genomic scores should cite their original source included in the metadata of the scores and accessible through the function `citation()`. For citing `GenomicScores` as a software package, please use the following reference:

   Puigdevall, P. and Castelo. R. GenomicScores: seamless access to genomewide position-specific scores from R and Bioconductor. _Bioinformatics_, 18:3208-3210, 2018.

## Installation

This is the __development__ version of the R/Bioconductor package GenomicScores. This version is unstable and should be used only to test new features. If you are looking for the __release__ version of this package please go to its package release landing page at [https://bioconductor.org/packages/GenomicScores](https://bioconductor.org/packages/GenomicScores) and follow the instructions there to install it.

If you were really looking for this development version, then to install it you
need first to install the development version of R that you can find at [https://cran.r-project.org](https://cran.r-project.org) and then type the following instructions from the R shell:

```r
install.packages("BiocManager")
BiocManager::install("GenomicScores", version="devel")
```

Alternatively, you can install it from GitHub using the [devtools](https://github.com/hadley/devtools "devtools") package.

```r
install.packages("devtools")
library(devtools)
install_github("rcastelo/GenomicScores")
```

## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **GenomicScores**
please use the [Bioconductor support site](https://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **GenomicScores**
please use the GitHub issues link at the top-right of this page
([https://github.com/rcastelo/GenomicScores/issues](https://github.com/rcastelo/GenomicScores/issues)).
