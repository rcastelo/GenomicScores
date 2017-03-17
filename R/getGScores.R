## adapted from Biobase::testBioCConnection()
.testConnection <- function (urladr) {
  curNetOpt <- getOption("internet.info")
  on.exit(options(internet.info = curNetOpt), add = TRUE)
  options(internet.info = 3)
  http <- as.logical(capabilities(what = "http/ftp"))
  if (http == FALSE)
    return(FALSE)
  fgahURL <- url("http://functionalgenomics.upf.edu/annotationhub")
  options(show.error.messages = FALSE)
  test <- try(readLines(fgahURL)[1])
  options(show.error.messages = TRUE)
  if (inherits(test, "try-error"))
    return(FALSE)
  else
    close(fgahURL)

  return(TRUE)
}

availableGScores <- function() {
  baseUrl <- "http://functionalgenomics.upf.edu/annotationhub"
  avgs <- character(0)

  if (!.testConnection(baseUrl)) {
    warning(sprintf("No internet connection to %s", baseUrl))
    return(avgs)
  }

  mainDirs <- AnnotationForge:::.getSubDirs(baseUrl)
  mainDirs <- sub("/", "", mainDirs)
  mainDirs <- mainDirs[nchar(mainDirs) > 0]
  if (length(mainDirs) < 1) {
    warning(sprintf("No available genomic scores at %s", baseUrl))
    return(avgs)
  }

  for (d in mainDirs) {
    subDirs <- AnnotationForge:::.getSubDirs(paste(baseUrl, d, sep="/"))
    subDirs <- sub("/", "", subDirs[grep(d, subDirs)])
    avgs <- c(avgs, subDirs)
  }

  avgs
}

getGScores <- function(x) {
  if (!is.character(x) || length(x) > 1)
    stop("'x' should be a character vector of length 1.")

  ah <- AnnotationHub()
  ah <- query(ah, x)

  if (length(ah) == 0)
    stop("'x' is not available. Please use 'availableGScores()' to know what genomic scores resources are available, or try to update the metadata of the annotation hub.")

  ## use the AnnotationHub metadata to figure out the correspondence
  ## between downloaded filenames and object names
  mdah <- mcols(ah)
  objnames <- mdah$title
  ahids <- rownames(mdah)
  fnames <- cache(ah[names(ah)])
  serializedobjs <- basename(fnames[ahids])
  names(serializedobjs) <- objnames

  ## load the first object to get the metadata
  obj <- ah[[1]]
  mdobj <- metadata(obj)
  gsco <- GScores(provider=mdobj$provider,
                  provider_version=mdobj$provider_version,
                  download_url=mdobj$download_url,
                  download_date=mdobj$download_date,
                  reference_genome=mdobj$reference_genome,
                  data_pkgname=mdobj$data_pkgname,
                  data_dirpath=getAnnotationHubOption("CACHE"),
                  data_serialized_objnames=serializedobjs)
  scorlelist <- get(mdobj$data_pkgname, envir=gsco@.data_cache)
  scorlelist[[mdobj$seqname]] <- obj
  assign(mdobj$data_pkgname, scorlelist, envir=gsco@.data_cache)

  gsco
}
