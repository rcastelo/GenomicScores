getGScores <- function(x, accept.license=FALSE) {
  if (!is.character(x) || length(x) > 1)
    stop("'x' should be a character vector of length 1.")

  if (!is.logical(accept.license) || length(accept.license) > 1)
    stop("'accept.license' should be a logical vector of length 1.")

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
  lic <- licurl <- character(0)
  licreqconsent <- FALSE
  if (!is.null(mdobj$license) && is.character(mdobj$license))
    lic <- mdobj$license
  if (!is.null(mdobj$license_url) && is.character(mdobj$license_url))
    licurl <- mdobj$license_url
  if (!is.null(mdobj$license_reqconsent) && is.logical(mdobj$license_reqconsent))
    licreqconsent <- mdobj$license_reqconsent
  createobj <- TRUE
  if (nchar(lic) > 0 && licreqconsent) {
    licensestr <- lic
    if (nchar(licurl) > 0)
      licensestr <- sprintf("%s (see %s)", lic, licurl)
    ans <- "n"
    if (!accept.license && interactive())
      ans <- .getAnswer(sprintf("These data is shared under the license %s, do you accept it? [y/n]: ", licensestr),
                        allowed=c("y", "Y", "n", "N"))
    else if (accept.license) {
      ans <- "y"
      message(sprintf("Using these data you are accepting the license %s\n", licensestr))
    }
    if (ans != "y")
      stop("Data will not be made available as a 'GScores' object unless you agree to the terms of its license.")
  }
  gsco <- NULL
  if (createobj) {
    gsco <- GScores(provider=mdobj$provider,
                    provider_version=mdobj$provider_version,
                    download_url=mdobj$download_url,
                    download_date=mdobj$download_date,
                    reference_genome=mdobj$reference_genome,
                    data_pkgname=mdobj$data_pkgname,
                    data_dirpath=getAnnotationHubOption("CACHE"),
                    data_serialized_objnames=serializedobjs,
                    license=lic, license_url=licurl,
                    license_reqconsent=licreqconsent)
    scorlelist <- get(mdobj$data_pkgname, envir=gsco@.data_cache)
    scorlelist[[mdobj$seqname]] <- obj
    assign(mdobj$data_pkgname, scorlelist, envir=gsco@.data_cache)
  }

  gsco
}
