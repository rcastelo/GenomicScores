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
  fnames <- cache(ah[ahids])

  ## in BioC 3.8 names from 'cache()' became 'AHXXXX : YYYY'
  ## so we need to override those names.
  names(fnames) <- ahids
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
