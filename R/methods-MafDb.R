## constructor
MafDb <- function(provider, provider_version, download_url,
                  download_date, reference_genome,
                  rsIDgpSNVs, rsIDSNVs, rsIDidxSNVs,
                  data_pkgname, data_dirpath, data_serialized_objnames,
                  data_tag="") {
  data_cache <- new.env(hash=TRUE, parent=emptyenv())
  data_pops <- list.files(path=data_dirpath, pattern=data_pkgname)
  data_pops <- data_pops[grep("snv", data_pops, invert=TRUE)]
  data_pops <- gsub(paste0(data_pkgname, "."), "", data_pops)
  data_pops <- gsub("\\..+$", "", data_pops)
  data_pops <- sort(unique(data_pops))
  data_serialized_objnames <- sub(".rds", "",
                                  list.files(path=data_dirpath, pattern="*.rds"))

  assign(data_pkgname, list(), envir=data_cache)
  assign(sprintf("%s.nonsnvs", data_pkgname), GRangesList(), envir=data_cache)

  nov <- NA_integer_
  if (file.exists(file.path(data_dirpath, "nov.rds")))
    nov <- as.integer(readRDS(file.path(data_dirpath, "nov.rds")))

  new("MafDb", provider=provider,
               provider_version=provider_version,
               download_url=download_url,
               download_date=download_date,
               reference_genome=reference_genome,
               data_pkgname=data_pkgname,
               data_dirpath=data_dirpath,
               data_serialized_objnames=data_serialized_objnames,
               data_tag=data_tag,
               data_pops=data_pops,
               data_nov=nov,
               .data_cache=data_cache)
}

## accessors
setMethod("provider", "MafDb", function(x) x@provider)

setMethod("providerVersion", "MafDb", function(x) x@provider_version)

setMethod("referenceGenome", "MafDb", function(x) x@reference_genome)

setMethod("organism", "MafDb", function(object) organism(referenceGenome(object)))

setMethod("seqinfo", "MafDb", function(x) seqinfo(referenceGenome(x)))

setMethod("seqnames", "MafDb", function(x) seqnames(referenceGenome(x)))

setMethod("seqlengths", "MafDb", function(x) seqlengths(referenceGenome(x)))

setMethod("seqlevelsStyle", "MafDb", function(x) seqlevelsStyle(referenceGenome(x)))

setMethod("populations", "MafDb", function(x) x@data_pops)

citation.MafDb <- function(package, ...) {
  obj <- get(package@data_pkgname, envir=package@.data_cache)
  cit <- metadata(obj[[1]][[1]])$citation
  if (is.null(cit))
    cit <- bibentry()
  cit
}
setMethod("citation", signature="MafDb", citation.MafDb)

## adapted from VariantTools::extractCoverageForPositions()
.extractRawFromRleList <- function(rlelst, ranges) {
  if (any(!unique(seqnames(ranges)) %in% names(rlelst)))
    stop("Some sequence names from input positions are missing from rlelst")
  if (any(width(ranges) > 1L))
    stop("Some ranges are of width > 1")
  seqlevels(ranges) <- names(rlelst)
  ord <- order(seqnames(ranges))
  ans <- raw(length(ranges))
  ans[ord] <- unlist(mapply(function(v, p) {
    runValue(v)[findRun(p, v)]
  }, rlelst, split(start(ranges), seqnames(ranges)), SIMPLIFY=FALSE), use.names=FALSE)
  ans
}

## adapted from BSgenome:::.normarg_ranges()
.str2gr <- function(ranges) {
  if (class(ranges) == "character") {
    spl1 <- strsplit(ranges, ":", fixed=TRUE)
    if (!all(sapply(spl1, length) == 2L))
      stop("Genomic ranges given in a character string should have the format CHR:START[-END]\n")
    ranges <- sapply(spl1, function(ranges) {
                             ss <- strsplit(ranges[2], "-", fixed=TRUE)[[1]]
                             if (length(ss) == 1L)
                               ss <- c(ss, ss)
                             c(ranges[1], ss)
                           })
    ranges <- t(ranges)
    ranges <- GRanges(seqnames=ranges[, 1],
                      IRanges(start=as.integer(ranges[, 2]),
                              end=as.integer(ranges[, 3])))
  } else if (class(ranges) == "VRanges") {
    ranges <- as(ranges, "GRanges")
    mcols(ranges) <- NULL
  } else if (!is(ranges, "GenomicRanges"))
    stop("argument 'ranges' must be either a GenomicRanges object or a character string with the format CHR:START[-END]\n")

  ranges
}

.mafByOverlaps_snvs <- function(x, ranges, snames, pop, caching) {
  mafsnvs <- get(x@data_pkgname, envir=x@.data_cache)
  missingMask <- !pop %in% names(mafsnvs)
  for (popname in pop[missingMask])
    mafsnvs[[popname]] <- RleList(compress=FALSE)
  anyMissing <- any(missingMask)

  ans <- DataFrame(as.data.frame(matrix(NA_real_, nrow=length(ranges), ncol=length(pop),
                                        dimnames=list(NULL, pop))))
  for (popname in pop) {
    missingMask <- !snames %in% names(mafsnvs[[popname]])
    anyMissing <- anyMissing || any(missingMask)
    for (sname in snames[missingMask]) {
      fname <- file.path(x@data_dirpath,
                         sprintf("%s.%s.%s.rds", x@data_pkgname, popname, sname))
      if (file.exists(fname))
        mafsnvs[[popname]][[sname]] <- readRDS(fname)
      else {
        warning(sprintf("No MAF data for chromosome %s.", sname))
        mafsnvs[[popname]][[sname]] <- Rle(raw(0))
      }
    }
    ans[[popname]] <- rep(NA_real_, nrow(ans))
    if (length(mafsnvs[[popname]]) > 0) {
      q <- .extractRawFromRleList(mafsnvs[[popname]], ranges)
      ans[[popname]] <- metadata(mafsnvs[[popname]][[1]])$dqfun(q)
    }
  }

  if (anyMissing && caching)
    assign(x@data_pkgname, mafsnvs, envir=x@.data_cache)
  rm(mafsnvs)

  ans
}

.mafByOverlaps_nonsnvs <- function(x, ranges, snames, pop, caching) {
  mafnonsnvs <- get(sprintf("%s.nonsnvs", x@data_pkgname), envir=x@.data_cache)
  mcnames <- character(0)
  if (length(mafnonsnvs) > 0)
    mcnames <- colnames(mcols(mafnonsnvs[[1]]))
  missingMask <- !snames %in% names(mafnonsnvs)
  anyMissing <- any(missingMask)
  for (sname in snames[missingMask]) {
    fname <- file.path(x@data_dirpath,
                       sprintf("%s.GRnonsnv.%s.rds", x@data_pkgname, sname))
    obj <- GRanges()
    if (file.exists(fname)) {
      obj <- readRDS(fname)
      for (cname in mcnames) {
        fname <- file.path(x@data_dirpath,
                           sprintf("%s.RLEnonsnv.%s.%s.rds", x@data_pkgname, cname, sname))
        if (file.exists(fname)) {
          q <- readRDS(fname)
          mcols(obj)[[cname]] <- metadata(q)$dqfun(q)
        } else
          stop(sprintf("internal file %s not found", fname))
      }
    } else {
      warning(sprintf("No MAF data for chromosome %s.", sname))
      mcols(obj) <- DataFrame(as.data.frame(matrix(raw(0), nrow=0, ncol=length(mcnames),
                                                   dimnames=list(NULL, mcnames))))
    }
    mafnonsnvs[[sname]] <- obj
  }

  ans <- DataFrame(as.data.frame(matrix(NA_real_, nrow=length(ranges), ncol=length(pop),
                                        dimnames=list(NULL, pop))))
  missingMask <- !pop %in% mcnames
  anyMissing <- anyMissing || any(missingMask)
  tmp <- unlist(mafnonsnvs)
  names(tmp) <- NULL
  for (popname in pop[missingMask]) {
    mafvalues <- numeric(length(tmp))
    i <- 1
    ## b/c we're storing MAF values as metadata columns of GRanges
    ## populations need to be loaded for all already loaded chromosomes
    for (sname in names(mafnonsnvs)) {
      fname <- file.path(x@data_dirpath,
                         sprintf("%s.RLEnonsnv.%s.%s.rds", x@data_pkgname, popname, sname))
      if (file.exists(fname)) {
        q <- readRDS(fname)
        mafvalues[i:(i+length(mafnonsnvs[[sname]])-1)] <- metadata(q)$dqfun(q)
      } else {
        mafvalues[i:(i+length(mafnonsnvs[[sname]])-1)] <- rep(NA_real_, length(mafnonsnvs[[sname]]))
      }
      i <- i + length(mafnonsnvs)
    }
    mcols(tmp)[[popname]] <- mafvalues
  }
  mafnonsnvs <- split(tmp, seqnames(tmp), drop=TRUE)
  rm(tmp)

  ## b/c MAF data is imported from VCF files and these represent insertions and
  ## deletions using the nucleotide composition of the reference sequence we
  ## use here minoverlap=1L but this may change in the future
  ov <- findOverlaps(ranges, unlist(mafnonsnvs), minoverlap=1L)
  if (length(ov) > 0) {
    ans[queryHits(ov), pop] <- mcols(unlist(mafnonsnvs))[subjectHits(ov), pop]
    if (any(duplicated(queryHits(ov))))
      message("mafByOverlaps: more than one variant overlapping queried positions, reporting only the first hit.")
  }

  if (anyMissing && caching)
    assign(sprintf("%s.nonsnvs", x@data_pkgname), mafnonsnvs, envir=x@.data_cache)
  rm(mafnonsnvs)

  ans
}

setMethod("mafByOverlaps", signature="MafDb",
          function(x, ranges, pop="AF", type=c("snvs", "nonsnvs"),
                   maf.only=FALSE, caching=TRUE) {
            type <- match.arg(type)
            ranges <- .str2gr(ranges)

            if (class(pop) != "character")
              stop("argument 'pop' must be a character vector")

            snames <- unique(as.character(runValue(seqnames(ranges))))
            if (any(!snames %in% seqnames(x)))
              stop(sprintf("Sequence names %s in 'ranges' not present in MafDb object.",
                   paste(snames[!snames %in% seqnames(x)], collapse=", ")))

            if (any(!pop %in% populations(x)))
              stop(sprintf("population %s must be one of %s\n", pop, paste(populations(x), collapse=", ")))

            ans <- NULL
            if (type == "snvs")
              ans <- .mafByOverlaps_snvs(x, ranges, snames, pop, caching)
            else ## nonsnvs
              ans <- .mafByOverlaps_nonsnvs(x, ranges, snames, pop, caching)

            if (!maf.only) {
              mcols(ranges) <- ans
              ans <- ranges
            }
            ans
          })

setMethod("mafById", signature="MafDb",
          function(x, ids, pop="AF", maf.only=FALSE, caching=TRUE) {
            if (class(ids) != "character")
              stop("argument 'ids' must be a character string vector.")

            if (!exists("rsIDs", envir=x@.data_cache)) {
              if (file.exists(file.path(x@data_dirpath, "rsIDs.rds"))) {
                message("Loading first time annotations of rs identifiers to variants, produced by data provider.")
                rsIDs <- readRDS(file.path(x@data_dirpath, "rsIDs.rds"))
                assign("rsIDs", rsIDs, envir=x@.data_cache)
              } else {
                warning("The data provider did not produce annotations of rs identifiers to variants.")
                ans <- DataFrame(as.data.frame(matrix(NA_real_, nrow=length(ids), ncol=length(pop),
                                                      dimnames=list(NULL, pop))),
                                 row.names=ids)
                return(ans)
              }
            }

            rsIDs <- get("rsIDs", envir=x@.data_cache)
            if (is.integer(rsIDs)) {
              idsint <- rep(NA_integer_, length(ids))
              rsMask <- regexpr("^rs", ids) == 1
              idsint[rsMask] <- as.integer(sub(pattern="^rs", replacement="", x=ids[rsMask]))
              mt <- rep(NA_integer_, length(idsint))
              if (any(!is.na(idsint))) {
                idsintnoNAs <- idsint[!is.na(idsint)]
                ord <- order(idsintnoNAs)                        ## order ids to speed up
                mtfi <- findInterval(idsintnoNAs[ord], rsIDs)    ## call to findInterval()
                mtfi[ord] <- mtfi                                ## put matches into original order
                mt[!is.na(idsint)] <- mtfi                       ## integrate matches into result
              }
              mt[mt == 0] <- 1
              if (any(!is.na(mt))) {
                maskNAs <- idsint != rsIDs[mt]
                mt[maskNAs] <- NA
              }
              if (any(!is.na(mt))) {
                if (!exists("rsIDidx", envir=x@.data_cache)) {
                  rsIDidx <- readRDS(file.path(x@data_dirpath, "rsIDidx.rds"))
                  assign("rsIDidx", rsIDidx, envir=x@.data_cache)
                }
                rsIDidx <- get("rsIDidx", envir=x@.data_cache)
                mt <- rsIDidx[mt]
              }
            } else
              stop("internal object 'rsIDs' is not an integer vector.")

            if (any(!pop %in% populations(x)))
              stop(sprintf("population %s must be one of %s\n", pop, paste(populations(x), collapse=", ")))

            ans <- DataFrame(as.data.frame(matrix(NA_real_, nrow=length(ids), ncol=length(pop),
                                                  dimnames=list(NULL, pop))),
                             row.names=ids)

            if (any(!is.na(mt))) {
              if (!exists("rsIDgp", envir=x@.data_cache)) {
                rsIDgp <- readRDS(file.path(x@data_dirpath, "rsIDgp.rds"))
                assign("rsIDgp", rsIDgp, envir=x@.data_cache)
              }
              rsIDgp <- get("rsIDgp", envir=x@.data_cache)
              issnvmask <- rsIDgp$isSNV
              rsIDgp <- updateObject(rsIDgp, verbose=FALSE) ## temporary solution until GPos objects are updated
              rsIDgp$isSNV <- issnvmask
              rm(issnvmask)
              rng <- rsIDgp[mt[!is.na(mt)]]
              mask <- logical(length(mt))
              mask[!is.na(mt)] <- rng$isSNV
              if (any(mask))
                ans[mask, pop] <- mafByOverlaps(x, rng[rng$isSNV], pop, type="snvs",
                                                maf.only=TRUE, caching)
              mask <- logical(length(mt))
              mask[!is.na(mt)] <- !rng$isSNV
              if (any(mask))
                ans[mask, pop] <- mafByOverlaps(x, rng[!rng$isSNV], pop, type="nonsnvs",
                                                maf.only=TRUE, caching)
            }

            ans
          })

## moved to methods-GScores
## .pprintseqs <- function(x) {
##   y <- x
##   if (length(x) > 5)
##     y <- c(y[1:2], "...", y[length(y)])
##   y <- paste(y, collapse=", ")
##   y 
## }

## show method
setMethod("show", "MafDb",
          function(object) {
            snvsobj <- get(object@data_pkgname, envir=object@.data_cache)
            nonsnvsobj <- get(paste0(object@data_pkgname, ".nonsnvs"),
                              envir=object@.data_cache)
            loadedsnvspops <- loadedsnvsseqs <- "none"
            loadednonsnvspops <- loadednonsnvsseqs <- "none"
            if (length(snvsobj) > 0) {
              loadedsnvspops <- names(snvsobj)
              loadedsnvsseqs <- names(snvsobj[[1]])
            }
            if (ncol(mcols(nonsnvsobj)) > 0)
              loadednonsnspops <- colnames(mcols(nonsnvsobj))
            if (length(nonsnvsobj) > 0)
              loadednonsnvsseqs <- unique(seqnames(nonsnvsobj))

            cat("Minor allele frequency Db (MafDb) object\n",
                "# organism: ", organism(object), "\n",
                "# provider: ", provider(object), "\n",
                "# provider version: ", providerVersion(object), "\n",
                "# download date: ", object@download_date, "\n",
                "# loaded sequences (SNVs): ", .pprintseqs(loadedsnvsseqs), "\n",
                "# loaded sequences (nonSNVs): ", .pprintseqs(loadednonsnvsseqs), "\n",
                "# loaded populations (SNVs): ", .pprintseqs(loadedsnvspops), "\n",
                "# loaded populations (nonSNVs): ", .pprintseqs(loadednonsnvspops), "\n",
                "# nr. of variants: ", object@data_nov, "\n",
                sep="")
          })

## $ method
setMethod("$", signature(x="MafDb"),
          function(x, name) {
            switch(name,
                   tag=x@data_tag,
                   stop("uknown MafDb slot.")
                   )
          })
