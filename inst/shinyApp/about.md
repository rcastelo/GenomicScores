## A Shiny App for the GenomicScore package

Welcome to our shiny gscores app.

The purpose of this app is to facilitate using the [GenomicScores package](https://bioconductor.org/packages/release/bioc/vignettes/GenomicScores/inst/doc/GenomicScores.html) with a simple GUI.

1. Choose which Annotation package (e.g., [phastCons100way.UCSC.hg38](http://bioconductor.org/packages/release/data/annotation/html/phastCons100way.UCSC.hg38.html) ) you want to use. The app will show you the `show()` and `citation()` function for that package.  
**NOTE**  
This app doesn't include any annotation package, but it identifies the ones you have already installed and the ones you don't. In order to install them, you can select them in "Select an Annotation Package" and a modal window will guide you through the installation process. You can also install them as normally, using the `BiocManager` package on your R session.  
Example:
```BiocManager::install("phastCons100way.UCSC.hg38")```

2. Select between the populations options from within your chosen package (e.g., `default` or `DP2` ). You can choose more than one population.

3. Then, choose between using this app interface to select parameters for the GRange object (`seqname`, `starting` and `ending` position), or you can upload your own [bed](https://genome.ucsc.edu/FAQ/FAQformat#format1) format file with this inputs in order to look up for their genomic scores.  
Bear in mind that by using the app parameters, you can choose to see the mean gscores of a range of positions (**Range**), or you can see them separately (**Individual**)

4. Finally, you can download the generated tables in bed or csv format.
