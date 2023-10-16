## adapted from Biobase::testBioCConnection()
.testConnection <- function (urladr) {
  curNetOpt <- getOption("internet.info")
  on.exit(options(internet.info = curNetOpt), add = TRUE)
  options(internet.info = 3)
  http <- as.logical(capabilities(what = "http/ftp"))
  if (http == FALSE)
    return(FALSE)
  fgahURL <- url(urladr)
  options(show.error.messages = FALSE)
  test <- try(readLines(fgahURL)[1])
  options(show.error.messages = TRUE)
  if (inherits(test, "try-error"))
    return(FALSE)
  else
    close(fgahURL)

  return(TRUE)
}

## adapted from AnnotationForge:::.getSubDirs()
.getSubDirs <- function(dname)
{
    getLinks <- function() {
        links <- character(0)
        list(a = function(node, ...) {
                   links <<- c(links, xmlGetAttr(node, "href"))
                   node
                 },
             links = function() links)
    }
    h1 <- getLinks()
    doc <- httr::content(httr::GET(dname))
    htmlTreeParse(doc, handlers = h1)
    res <- h1$links()
    res <- res[!(res %in% c("?C=N;O=D", "?C=M;O=A", "?C=S;O=A", "?C=D;O=A",
                            "/download/current/"))]
    res
}

.availableGScoresAH <- function(use.internet=FALSE) {
  avgs <- readRDS(system.file("extdata", "avgs.rds", package="GenomicScores"))

  if (use.internet) {
    baseUrl <- "https://functionalgenomics.upf.edu/annotationhub"
    if (!.testConnection(baseUrl)) {
      stop(sprintf("No internet connection to %s", baseUrl))
    } else {
      mainDirs <- .getSubDirs(baseUrl)
      mainDirs <- sub("/", "", mainDirs)
      mainDirs <- mainDirs[nchar(mainDirs) > 0]
      if (length(mainDirs) < 1)
        stop(sprintf("No available genomic scores at %s", baseUrl))

      avgs <- character(0)
      for (d in mainDirs) {
        subDirs <- .getSubDirs(paste(baseUrl, d, sep="/"))
        subDirs <- sub("/", "", subDirs[grep(d, subDirs)])
        avgs <- c(avgs, subDirs)
      }
    }
  }

  ahresincache <- rep(FALSE, length(avgs))
  if (file.exists(getAnnotationHubOption("CACHE")) || use.internet) {

    ## fetch information about cached GScore resources
    suppressMessages(ah <- AnnotationHub(localHub=!use.internet))
    ah <- query(ah, avgs, pattern.op=`|`)
    mcah <- mcols(ah)
    avgs <- avgs[!is.na(charmatch(avgs, mcah$title))]
    mt <- regexpr(paste(avgs, collapse="|"), mcah$title)
    stopifnot(all(mt == 1)) ## QC
    mcah$resname <- substr(mcah$title, 1, attr(mt, "match.length"))
    ahidsbyresname <- split(rownames(mcah), mcah$resname)

    ## fetch AH ids in the cache
    bfc <- BiocFileCache(hubCache(ah))
    cachedres <- bfcinfo(bfc)
    cachedres <- sub(" : [0-9]+", "", cachedres$rname)
    ahresincache <- sapply(ahidsbyresname,
                           function(ahids, cachedids) all(ahids %in% cachedids), cachedres)
  }

  data.frame(Name=avgs, Cached=ahresincache[avgs], stringsAsFactors=FALSE)
}

availableGScores <- function(use.internet=FALSE) {

  ahpkgs <- .availableGScoresAH(use.internet)

  pp <- system.file("scripts", package="GenomicScores")
  mkdatafnames <- list.files(pp, pattern="make-data_*")
  gspkgnames <- sub("make-data_", "", mkdatafnames, fixed=TRUE)
  gspkgnames <- sub(".R", "", gspkgnames, fixed=TRUE)
  ip <- rownames(installed.packages())

  ## the BioC core team wants that the newly added AH GenomicScores resources
  ## also have corresponding lightweight annotation packages, but those packages
  ## do not create any GScores object at loading time, so by now we removed them
  ## hardcoding the package names from the list of GenomicScores annotation
  ## packages available through install.
  lightweightpkgs <- c("phastCons30way.UCSC.hg38", "phastCons35way.UCSC.mm39",
                       "phyloP35way.UCSC.mm39", "cadd.v1.6.hg19",
                       "cadd.v1.6.hg38", "AlphaMissense.v2023.hg19",
                       "AlphaMissense.v2023.hg38")
  ip <- setdiff(ip, lightweightpkgs)

  cached <- setdiff(ahpkgs$Name[ahpkgs$Cached], ip)

  ## since package membership of a release does not change within a release,
  ## by default (use.internet=FALSE) we load a pre-saved package list to avoid
  ## investing time checking through the internet
  avgsbmi <- readRDS(system.file("extdata", "avgsbmi.rds", package="GenomicScores"))
  if (use.internet)
    avgsbmi <- BiocManager::available(pattern=paste(gspkgnames, collapse="|"))

  res <- data.frame(Organism=rep(NA, length(gspkgnames)),
                    Category=rep(NA, length(gspkgnames)),
                    Installed=gspkgnames %in% ip,
                    Cached=gspkgnames %in% cached,
                    BiocManagerInstall=gspkgnames %in% avgsbmi,
                    AnnotationHub=(gspkgnames %in% ahpkgs$Name) & (!gspkgnames %in% avgsbmi),
                    row.names=gspkgnames,
                    stringsAsFactors=FALSE)

  ## read frozen GScores resources metadata
  gsrm <- read.csv(gzfile(system.file("extdata", "GScoresResourcesMetadata.csv.gz",
                               package="GenomicScores")), row.names=1)
  stopifnot(all(colnames(gsrm) == c("Organism", "Category"))) ## QC
  mt <- match(rownames(gsrm), rownames(res))
  stopifnot(all(!is.na(mt))) ## QC
  res$Organism[mt] <- gsrm$Organism
  res$Category[mt] <- gsrm$Category
  ## if (any(res$Installed)) {
  ##   orggrp <- sapply(rownames(res)[res$Installed],
  ##                    function(pkg) {
  ##                      obj <- getFromNamespace(pkg, pkg)
  ##                      unloadNamespace(pkg)
  ##                      c(organism(obj), gscoresCategory(obj))
  ##                    })
  ##   res[res$Installed, c("Organism", "Category")] <- t(orggrp)
  ## }

  res
}
