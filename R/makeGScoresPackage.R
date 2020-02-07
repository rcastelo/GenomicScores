## the following functions .normMaintainer, .getMaintainer
## .mergeMaintainer and .normAuthors are taken from
## GenomicFeatures/R/makeTxDbPackage.R
.normMaintainer <- function(maintainer) {
    maintainer <- as.person(maintainer)
    if (length(maintainer) > 1L) {
        stop("more than one 'maintainer' provided")
    }
    maintainer
}

.getMaintainer <- function(authors) {
    m <- vapply(authors, function(a) "cre" %in% a$role, logical(1L))
    if (sum(m) != 1L) {
        stop("there must be one 'maintainer'")
    }
    maintainer <- authors[m]
    maintainer$role <- list(NULL)
    maintainer$comment <- list(NULL)
    maintainer
}

.mergeMaintainer <- function(authors, maintainer) {
    maintainer <- .normMaintainer(maintainer)
    maintainer$role <- list(union(maintainer$role, "cre"))
    m <- unlist(authors$given) == maintainer$given &
        unlist(authors$family) == maintainer$family
    if (any(m)) {
        authors$role[m] <- list(union(unlist(authors$role[m]), "cre"))
        if (!is.null(maintainer$email)) {
            authors$email[m] <- maintainer$email
        }
    } else {
        authors <- c(authors, maintainer)
    }
    maintainer <- .getMaintainer(authors)
    if (is.null(maintainer$email)) {
        stop("the 'maintainer' must have an email address")
    }
    authors
}

.normAuthor <- function(authors, maintainer) {
  authors <- as.person(authors)
  if (!missing(maintainer)) {
    authors <- .mergeMaintainer(authors, maintainer)
  }
  authors
}

## adapted from makeTxDbPackage() in GenomicFeatures/R/makeTxDbPackage.R
makeGScoresPackage <- function(gsco, version, maintainer, author,
                               destDir=".", license="Artistic-2.0") {
  pkgname <- name(gsco)
  authors <- .normAuthor(author, maintainer)

  template_path <- system.file("gscores-template", package="GenomicScores")

  symvals <- list(PKGTITLE=sprintf("%s genomic scores for %s (%s)", type(gsco),
                                   organism(gsco), providerVersion(genomeDescription(gsco))),
                  PKGDESCRIPTION=sprintf("Store %s genomic scores.", type(gsco)),
                  PKGVERSION=version,
                  AUTHOR=paste(authors, collapse=", "),
                  MAINTAINER=as.character(.getMaintainer(authors)),
                  GSVERSION=packageDescription("GenomicScores")$Version,
                  LIC=license,
                  ORGANISM=organism(gsco),
                  SPECIES=organism(gsco),
                  GENOMEVERSION=providerVersion(genomeDescription(gsco)),
                  PROVIDER=provider(gsco),
                  PROVIDERVERSION=providerVersion(gsco),
                  ORGANISMBIOCVIEW=gsub(" ", "_", organism(gsco))
                 )

  res <- createPackage(pkgname=pkgname, destinationDir=destDir,
                       originDir=template_path, symbolValues=symvals)

  if (length(citation(gsco)) > 0)
    writeLines(capture.output(print(citation(gsco), style="R")),
               file.path(destDir, pkgname, "inst", "CITATION"))

  data_dirpath <- gsco@data_dirpath
  refgenomeGD <- metadata(get(pkgname, envir=gsco@.data_cache)[[1]])$reference_genome
  saveRDS(refgenomeGD,
          file=file.path(destDir, pkgname, "inst", "extdata", "refgenomeGD.rds"))

  data_files <- list.files(pattern=pkgname, path=data_dirpath, full.names=TRUE)
  source_data_files <- gsco@data_serialized_objnames
  target_data_files <- names(gsco@data_serialized_objnames)
  if (length(source_data_files) == 0 || length(target_data_files) == 0)
    stop("Cannot find data files")

  for (i in seq_along(source_data_files)) {
    copied <- file.copy(from=file.path(data_dirpath, source_data_files[i]),
                        to=file.path(destDir, pkgname, "inst", "extdata", target_data_files[i]))
    if (!copied)
      stop(sprintf("Cannot write in the %s directory",
                   file.path(destDir, pkgname, "inst", "extdata")))
  }

  invisible(res$pkgdir)
}
