.availableGScoresAH <- function() {
  baseUrl <- "http://functionalgenomics.upf.edu/annotationhub"
  avgs <- character(0)

  if (!.testConnection(baseUrl)) {
    warning(sprintf("No internet connection to %s", baseUrl))
    ## in the absence of internet connectivity, assume all possibly downloaded
    ## GScores resources are the ones cached in inst/extdata/avgs.rds
    avgs <- readRDS(system.file("extdata", "avgs.rds", package="GenomicScores"))
  } else {
    mainDirs <- .getSubDirs(baseUrl)
    mainDirs <- sub("/", "", mainDirs)
    mainDirs <- mainDirs[nchar(mainDirs) > 0]
    if (length(mainDirs) < 1) {
      warning(sprintf("No available genomic scores at %s", baseUrl))
      return(avgs)
    }

    for (d in mainDirs) {
      subDirs <- .getSubDirs(paste(baseUrl, d, sep="/"))
      subDirs <- sub("/", "", subDirs[grep(d, subDirs)])
      avgs <- c(avgs, subDirs)
    }
  }

  ## report only those in the current AnnotationHub database
  ah <- AnnotationHub()
  ah <- query(ah, avgs, pattern.op=`|`)
  mcah <- mcols(ah)

  avgs[!is.na(charmatch(avgs, mcah$title))]
}

availableGScores <- function(installed=FALSE, ah.only=FALSE) {

  ahgspkgs <- readRDS(system.file("extdata", "avgs.rds",
                                  package="GenomicScores"))
  pp <- system.file("scripts", package="GenomicScores")
  mkdatafnames <- list.files(pp, pattern="make-data_*")
  gspkgnames <- sub("make-data_", "", mkdatafnames, fixed=TRUE)
  gspkgnames <- sub(".R", "", gspkgnames, fixed=TRUE)
  ip <- installed.packages()
  res <- data.frame(Name=gspkgnames,
                    Organism=rep(NA, length(gspkgnames)),
                    Category=rep(NA, length(gspkgnames)),
                    Installed=gspkgnames %in% ip,
                    AH=gspkgnames %in% ahgspkgs)
  if (ah.only) {
    res <- res[res$AH, ]
    ahpkgs <- .availableGSscoresAH()
    res <- res[res$Name %in% ahpkgs, ]
  }

  if (installed)
    res <- res[res$Installed, ]

  if (any(res$Installed)) {
    orggrp <- sapply(res$Name[res$Installed],
                     function(pkg) {
                       obj <- getFromNamespace(pkg, pkg)
                       unloadNamespace(pkg)
                       c(organism(obj), gscoresCategory(obj))
                     })
    res[res$Installed, c("Organism", "Category")] <- t(orggrp)
  }


  res
}
