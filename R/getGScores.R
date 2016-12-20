.availableGScores <- c("phastCons100way.UCSC.hg19",
                       "phastCons100way.UCSC.hg38",
                       "phastCons7way.UCSC.hg38",
                       "fitCons.UCSC.hg19")

availableGScores <- function() .availableGScores

getGScores <- function(x) {
  if (!is.character(x) || length(x) > 1)
    stop("'x' should be a character vector of length 1.")

  if (!x %in% availableGScores())
    stop("'x' is not available. Please use 'availableGScores()' to know what genomic scores sources are avilable")

  ## query the AnnotationHub
  ## FIXME: some error manipulation is needed when
  ##        the resource is not available for whatever reason
  ah <- AnnotationHub()
  ah <- query(ah, x)

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
  scorlelist <- get(gsco@data_pkgname, envir=gsco@.data_cache)
  scorlelist[[mdobj$seqname]] <- obj
  assign(mdobj$data_pkgname, scorlelist, envir=gsco@.data_cache)

  gsco
}
