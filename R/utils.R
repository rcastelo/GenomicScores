## this is, at the moment, an internal utility function to
## process RDS files storing MAF genomic scores with Rle vector
## and put them together by population in a single HDF5 file
## arguments: prefix - prefix of the RDS files to be processed
##            inputpath - path where the RDS files
##            outputpath - path where to store the resulting HDF5 files
.makeMafH5 <- function(prefix, inputpath, outputpath) {
  
  ## get populations
  data_pops <- list.files(path=inputpath, pattern=prefix)
  data_pops <- data_pops[grep("nonsnv", data_pops, invert=TRUE)]
  data_pops <- gsub(paste0(prefix, "."), "", data_pops)
  data_pops <- gsub("\\..+$", "", data_pops)
  data_pops <- sort(unique(data_pops))
  cat(sprintf("identified populations: %s\n", paste(data_pops, collapse=", ")))

  if (!dir.exists(inputpath))
    stop(sprintf("input path %s does not exist.", inputpath))

  ## retrieve common metadata, i.e., discarding chromosome-specific metadata
  first_file <- list.files(path=inputpath, pattern=prefix, full.names=TRUE)[1]
  if (!file.exists(first_file))
    stop(sprintf("cannot access any file starting with prefix %s in %s.", prefix, inputpath))

  first_chrom <- readRDS(first_file)
  chrspecificmetadata <- c("seqname", "ecdf", "max_abs_error", "maskREF")
  mt <- match(chrspecificmetadata, names(metadata(first_chrom)))
  stopifnot(all(!is.na(mt))) ## QC
  pkg_metadata <- metadata(first_chrom)[-mt]
  if (!is.null(pkg_metadata$data_pkgname))
    pkg_metadata$data_pkgname <- gsub("MafDb", "MafH5", pkg_metadata$data_pkgname)
  
  ## check if nsites.rds file exists, then save it in metadata
  nsitesfiles <- list.files(path=inputpath, pattern="nsites.rds$", full.names=TRUE)
  pkg_metadata[["nsites"]] <- 0L
  for (fname in nsitesfiles) {
    nsites <- readRDS(fname)
    pkg_metadata[["nsites"]] <- pkg_metadata[["nsites"]] + nsites
  }
  
  if (!dir.exists(outputpath))
    dir.create(outputpath)

  outputprefix <- gsub("MafDb", "MafH5", prefix)
  
  ## save metadata as rds file
  saveRDS(pkg_metadata, file = file.path(outputpath, paste0(outputprefix, ".metadata.rds")))
  
  ## save some memory
  rm(first_file, first_chrom, pkg_metadata)
  gc()

  for (pop in data_pops) {
    
    ## get chromosomes
    prefix2 <- sprintf("%s\\.%s\\.", prefix, pop)
    data_chroms <- list.files(path=inputpath,
                              pattern=prefix2)
    data_chroms <- gsub(prefix2, "", data_chroms)
    data_chroms <- gsub("\\.rds", "", data_chroms)
    
    ## create h5 file
    fname <- file.path(outputpath, paste0(outputprefix, ".", pop, ".h5"))
    h5createFile(fname)
    
    ## save each chromosome as a group inside the h5 file
    for (chr in data_chroms) {
      cat(sprintf("Processing chromosome %s from population %s.\n", chr, pop))

      h5createGroup(fname, group=chr)
      gsco <- readRDS(file.path(inputpath, sprintf("%s.%s.%s.rds", prefix, pop, chr)))
      maskREF <- RleArray(metadata(gsco)$maskREF, length(metadata(gsco)$maskREF))
      max_abs_error <- metadata(gsco)$max_abs_error
      scores <- RleArray(Rle(decode(gsco)), length(Rle(decode(gsco))))
      writeHDF5Array(x=maskREF, filepath=fname, name=paste0(chr,"/maskREF"),
                     level=9, verbose = TRUE)
      writeHDF5Array(x = max_abs_error, filepath=fname, name=paste0(chr,"/max_abs_error"),
                     level=9, verbose=TRUE)
      writeHDF5Array(x = scores, filepath=fname, name=paste0(chr, "/scores"),
                     level=9, verbose=TRUE)
      rm(gsco, maskREF, scores)
    }
  }
}
