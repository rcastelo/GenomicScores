---
title: "About"
author: "Pablo Rodriguez"
date: "11/27/2019"
output: html_document
---


## A Shiny App for the GenomicScore package

Welcome to our shiny gscore app.

The purpose of this app is to facilitate using the [GenomicScores package](https://bioconductor.org/packages/release/bioc/vignettes/GenomicScores/inst/doc/GenomicScores.html) with a simple GUI.

* First, you have to choose which Annotation package (e.g., [phastCons100way.UCSC.hg38](http://bioconductor.org/packages/release/data/annotation/html/phastCons100way.UCSC.hg38.html) ) you want to use. The app will show you the `show()` and `citation()` function for that package.

**NOTE** ****
This app doesn't include any annotation package, it looks up for your libraries and retrieves the ones you have installed in your machine. In order to install them you can refer to their Bioconductor pages or manual.

Example:
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phastCons100way.UCSC.hg38")

```

* Second, you can choose between the populations options from within your chosen package (e.g., `default` or `DP2` ). You can choose more than one population.

* Then, you have to choose between using this app interface to select parameters for the GRange object (`seqname`, `starting` and `ending` position), or you can upload your own [bed](https://genome.ucsc.edu/FAQ/FAQformat#format1) format file with this inputs in order to look up for their genomic scores.

Bear in mind that by using the app parameters, you can choose to see the mean gscore of a range of positions ( **Range** ), or you can see them separately ( **Individual** )

* Finally, you can download the generated tables in bed or csv format.
