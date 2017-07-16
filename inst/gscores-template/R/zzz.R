.onLoad <- function(libname, pkgname) {
  extdata_dirpath <- system.file("extdata", package=pkgname,
                                 lib.loc=libname, mustWork=TRUE)

  ## load GenomeDescription object frozen from BSgenome.Hsapiens.UCSC.hg19
  rg <- readRDS(file.path(extdata_dirpath, "refgenomeGD.rds"))

  ## fetch the metadata from the smallest file to minimize loading time
  finfo <- file.info(list.files(pattern=pkgname, path=extdata_dirpath,
                                full.names=TRUE))
  obj <- readRDS(file.path(extdata_dirpath, basename(rownames(finfo)[which.min(finfo$size)])))
  mdobj <- metadata(obj)
  stopifnot(identical(pkgname, mdobj$data_pkgname)) ## QC

  serializedobjs <- rownames(finfo)
  names(serializedobjs) <- serializedobjs

  ## make and export GScores object.
  gsco <- GScores(provider=mdobj$provider,
                  provider_version=mdobj$provider_version,
                  download_url=mdobj$download_url,
                  download_date=mdobj$download_date,
                  reference_genome=mdobj$reference_genome,
                  data_pkgname=pkgname,
                  data_dirpath=extdata_dirpath,
                  data_serialized_objnames=serializedobjs)
  scorlelist <- get(pkgname, envir=gsco@.data_cache)
  scorlelist[[mdobj$seqname]] <- obj
  assign(pkgname, scorlelist, envir=gsco@.data_cache)

  ns <- asNamespace(pkgname)
  assign(pkgname, gsco, envir=ns)

  namespaceExport(ns, pkgname)
}
