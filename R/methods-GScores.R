## GScores constructor
GScores <- function(provider, provider_version, download_url,
                    download_date, reference_genome,
                    data_pkgname, data_dirpath,
                    data_serialized_objnames=character(0),
                    default_pop="default",
                    data_group=sub("\\..*$", "", data_pkgname),
                    data_tag=sub("\\..*$", "", data_pkgname),
                    data_nsites=NA_real_,
                    data_hdf5=FALSE) {
  data_cache <- new.env(hash=TRUE, parent=emptyenv())
  data_pops <- list.files(path=data_dirpath, pattern=data_pkgname)
  data_nonsnrs <- length(grep("nonsnv", data_pops)) > 0
  data_pops <- data_pops[grep("nonsnv", data_pops, invert=TRUE)]
  data_pops <- gsub(paste0(data_pkgname, "."), "", data_pops)
  data_pops <- gsub("\\..+$", "", data_pops)
  data_pops <- sort(unique(data_pops))
  if (all(data_pops %in% seqlevels(reference_genome))) ## no additional score populations
    data_pops <- default_pop
  else {
    if (any(data_pops %in% seqlevels(reference_genome)))
      data_pops <- c(default_pop, setdiff(data_pops, seqlevels(reference_genome)))
  }

  assign(data_pkgname, list(), envir=data_cache)
  assign(sprintf("%s.nonsnvs", data_pkgname), GRangesList(), envir=data_cache)

  if (file.exists(file.path(data_dirpath, "nsites.rds")))
    data_nsites <- as.numeric(readRDS(file.path(data_dirpath, "nsites.rds")))

  new("GScores", provider=provider,
                 provider_version=provider_version,
                 download_url=download_url,
                 download_date=download_date,
                 reference_genome=reference_genome,
                 data_pkgname=data_pkgname,
                 data_dirpath=data_dirpath,
                 data_serialized_objnames=data_serialized_objnames,
                 data_group=data_group,
                 data_tag=data_tag,
                 data_pops=data_pops,
                 default_pop=default_pop,
                 data_nonsnrs=data_nonsnrs,
                 data_nsites=data_nsites,
                 data_hdf5=data_hdf5,
                 .data_cache=data_cache)
}

## accessors
setMethod("name", "GScores", function(x) x@data_pkgname)

setMethod("type", "GScores", function(x) sub("\\..*$", "", name(x)))

setMethod("provider", "GScores", function(x) x@provider)

setMethod("providerVersion", "GScores", function(x) x@provider_version)

setMethod("genomeDescription", "GScores", function(x) x@reference_genome)

setMethod("organism", "GScores",
          function(object) organism(genomeDescription(object)))

setMethod("seqinfo", "GScores", function(x) seqinfo(genomeDescription(x)))

setMethod("seqnames", "GScores", function(x) seqnames(genomeDescription(x)))

setMethod("seqlengths", "GScores", function(x) seqlengths(genomeDescription(x)))

setMethod("seqlevelsStyle", "GScores",
          function(x) seqlevelsStyle(genomeDescription(x)))

setMethod("gscoresNonSNRs", "GScores",
          function(x) x@data_nonsnrs)

setMethod("populations", "GScores", function(x) x@data_pops)

setMethod("defaultPopulation", "GScores", function(x) x@default_pop)

setReplaceMethod("defaultPopulation", c("GScores", "character"),
                 function(x, value) {
                   if (length(value) > 1)
                     message("more than one default scores population name supplied, using the 1st one only.")
                   value <- value[1]
                   if (any(!value %in% populations(x)))
                     stop(sprintf("scores population %s is not part of the available scores populations. Use 'populations()' to figure out which ones are available.", value))
                   x@default_pop <- value
                   x
                 })

setMethod("gscoresTag", "GScores", function(x) x@data_tag)

setReplaceMethod("gscoresTag", c("GScores", "character"),
                 function(x, value) {
                   if (length(value) > 1)
                     message("more than one genomic scores tag name supplied, using the 1st one only.")
                   x@data_tag <- value[1]
                   x
                 })

setMethod("gscoresCategory", "GScores", function(x) x@data_group)

setReplaceMethod("gscoresCategory", c("GScores", "character"),
                 function(x, value) {
                   if (length(value) > 1)
                     message("more than one genomic scores category name supplied, using the 1st one only.")
                   x@data_group <- value[1]
                   x
                 })

setMethod("nsites", "GScores", function(x) x@data_nsites)

setMethod("hdf5Backend", "GScores", function(x) x@data_hdf5)

## convert digits in vector 'd', grouped by 'g', into numbers in base 'b'
.toBase <- function(d, g=rep(1, length(d)), b) {
  n <- tapply(d, g, function(d, b) sum(d * b^(0:(length(d)-1))), b)
  as.integer(n)
}

## convert each number in 'n' into 'd' digits in base 'b'
.fromBase <- function(n, d, b) {
  totaldigits <- length(n) * d
  digits <- rep(NA_integer_, times=totaldigits)
  i <- 0L
  while (i < d) {
    digits[seq(1L+i, totaldigits, by=d)] <- n %% b
    n <- floor(n / b)
    i <- i + 1L
  }
  digits
}

## this has been improved using RleViews as discussed in
## https://stat.ethz.ch/pipermail/bioconductor/2013-December/056409.html
.rleGetValues <- function(rlelst, gr, mdata, summaryFun,
                          quantized=FALSE, isHDF5=FALSE) {
  summaryFun <- match.fun(summaryFun)
  numericmean <- TRUE
  if (!identical(summaryFun, mean))
    numericmean <- FALSE

  whregions <- which(width(gr) > 1)
  if (quantized && length(whregions) > 0)
    stop("When 'quantized=TRUE' input 'ranges' must have width one because summarization can only be calculated with dequantized values.")

  .dequantizer <- mdata$dqfun
  dqargs <- mdata$dqfun_args
  seqlevels(gr) <- names(rlelst)
  ord <- order(seqnames(gr)) ## store ordering below in 'split()'
  startbyseq <- split(start(gr), seqnames(gr), drop=TRUE)
  x <- ans <- numeric(0)
  if (quantized)
    x <- decode(unlist(rlelst[startbyseq], use.names=FALSE))
  else {
    lappargs <- c(list(X=rlelst[startbyseq], FUN=.dequantizer), dqargs)
    x <- unlist(do.call("lapply", lappargs), use.names=FALSE)
  }
  valxpos <- mdata$valxpos
  if (is.null(valxpos))
    valxpos <- 1
  stopifnot(length(valxpos) == 1 && valxpos[1] > 0) ## QC

  if (valxpos == 1) {
    ans <- vector(mode=class(x), length=length(gr))
    ans[ord] <- x
    if (quantized)
      ans <- Rle(ans)
  } else {
    if (length(whregions) > 0)
      stop("This GScores object returns more than one value per position, and therefore, input ranges must have width one.")
    ans <- matrix(NA_real_, nrow=length(gr), ncol=valxpos)
    ans[ord, ] <- matrix(x, nrow=length(gr), ncol=valxpos, byrow=TRUE)
  }

  if (length(whregions) > 0) { ## regions comprising more than
    tmpans <- NA_real_         ## one position are summarized
    ord2 <- order(seqnames(gr)[whregions]) ## store ordering below in 'split()'
    if (numericmean && !isHDF5) {
      rngbyseq <- split(gr[whregions], seqnames(gr)[whregions], drop=TRUE)
      tmpans <- lapply(names(rngbyseq),
                       function(sname) {
                         coercedrle <- rlelst[[sname]]
                         dqargs2 <- c(list(q=runValue(coercedrle)), dqargs)
                         runValue(coercedrle) <- do.call(".dequantizer", dqargs2)
                         viewMeans(Views(coercedrle,
                                         start=start(rngbyseq[[sname]]),
                                         end=end(rngbyseq[[sname]])))
                       })
      tmpans <- unsplit(tmpans, seqnames(gr)[whregions])
    } else { ## this allows for other summary functions
             ## but it runs about 10-fold slower
      startbyseq <- split(start(gr)[whregions],
                          seqnames(gr)[whregions], drop=TRUE)
      widthbyseq <- split(width(gr)[whregions],
                          seqnames(gr)[whregions], drop=TRUE)
      f <- function(r, p, w)
             sapply(seq_along(p),
                    function(i, r, p, w) {
                      dqargs2 <- c(list(q=r[p[i]:(p[i]+w[i]-1)]), dqargs)
                      summaryFun(do.call(".dequantizer", dqargs2))
                    },
                    r, p, w)

      tmpans <- unlist(mapply(f,
                              rlelst[names(startbyseq)], startbyseq, widthbyseq,
                              SIMPLIFY=FALSE),
                       use.names=FALSE)
    }

    ans[ord[whregions]] <- tmpans[ord2]
  }
  ans
}

.rleGetMetaValues <- function(rlelst, gr, metadataID, isHDF5=FALSE) {
  stopifnot(all(width(gr) == 1)) ## QC
  if (isHDF5)
    rlelst <- endoapply(rlelst,
                        function(x)
                          HDF5Array(path(seed(x)),
                                    name=sprintf("%s/%s",
                                                 sub("/scores", "",
                                                     seed(x)@name, metadataID))))
  else
    rlelst <- endoapply(rlelst, function(x) metadata(x)[[metadataID]])

  if (all(sapply(rlelst, length) == 0L))
    return(Rle(rep(FALSE, length(gr))))

  seqlevels(gr) <- names(rlelst)
  ord <- order(seqnames(gr)) ## store ordering below in 'split()'
  startbyseq <- split(start(gr), seqnames(gr), drop=TRUE)
  x <- as.logical(decode(unlist(rlelst[startbyseq], use.names=FALSE)))
  ans <- vector(mode=class(x), length=length(gr))
  ans[ord] <- x
  ans <- Rle(ans)

  ans
}

setMethod("score", "GScores",
          function(x, ..., simplify=TRUE) {
            gsco <- gscores(x, ..., scores.only=TRUE)
            if (ncol(gsco) == 1 && simplify)
              gsco <- gsco[[1]]
            gsco
          })

setMethod("gscores", c("GScores", "GenomicRanges"),
          function(x, ranges, ...) {
            ## default non-generic arguments
            paramNames <- c("scores.only", "pop", "type",
                            "summaryFun", "quantized",
                            "ref", "alt", "minoverlap", "caching")
            pop <- defaultPopulation(x)
            type <- "snrs"
            scores.only <- FALSE
            summaryFun <- mean
            quantized <- FALSE
            ref <- character(0)
            alt <- character(0)
            minoverlap <- 1L ## only relevant for genomic scores associated with nonSNRs
            caching <- TRUE

            ## get arguments
            arglist <- list(...)
            mask <- nchar(names(arglist)) == 0
            if (any(mask))
              names(arglist)[mask] <- paste0("X", 1:sum(mask))

            mask <- names(arglist) %in% paramNames
            if (any(!mask))
                stop(sprintf("unused argument (%s)", names(arglist)[!mask]))
            list2env(arglist, envir=sys.frame(sys.nframe()))

            if (!type %in% c("snrs", "nonsnrs"))
              stop("argument 'type' must be either 'snrs' (default) or 'nonsnrs'.")

            mask <- pop %in% populations(x)
            if (any(!mask))
              stop(sprintf("scores population %s is not present in %s. Please use 'populations()' to find out the available ones.",
                           pop[!mask], name(x)))

            ## adapt to sequence style and genome version from the input
            ## GScores object, thus assuming positions are based on the same
            ## genome even though might be differently specified
            ## (i.e., hg19 vs hs37d5 or hg38 vs GRCh38)
            if (length(intersect(seqlevelsStyle(ranges), seqlevelsStyle(x))) == 0)
              seqlevelsStyle(ranges) <- seqlevelsStyle(x)[1]
            commonSeqs <- intersect(seqlevels(ranges), seqlevels(x))
            if (any(is.na(genome(ranges)))) {
              ## message(sprintf("assuming query ranges genome build is the one of the GScores object (%s).",
              ##                 unique(genome(x)[commonSeqs])))
              genome(ranges) <- genome(x)
            } else if (any(genome(ranges)[commonSeqs] != genome(x)[commonSeqs])) {
              message(sprintf("assuming %s represent the same genome build between query ranges and the GScores object, respectively.",
                              paste(c(unique(genome(ranges)[commonSeqs]),
                                      unique(genome(x)[commonSeqs])),
                                    collapse=" and ")))
              genome(ranges) <- genome(x)
            }

            ord <- 1:length(ranges)
            if (is.unsorted(ranges)) {
              ord <- order(ranges)
              ranges <- ranges[ord]
              if (length(ref) > 0 && length(alt) > 0) {
                ref <- ref[ord]
                alt <- alt[ord]
              }
            }

            ans <- NULL
            if (type == "snrs")
              ans <- .scores_snrs(x, ranges, pop, summaryFun, quantized,
                                  scores.only, ref, alt, caching)
            else
              ans <- .scores_nonsnrs(x, ranges, pop, quantized, scores.only,
                                     ref, alt, minoverlap, caching)
            if (is(ans, "DFrame"))
              ans[ord, ] <- ans
            else ## GRanges
              ans[ord] <- ans

            ans
          })

setMethod("gscores", c("GScores", "character"),
          function(x, ranges, ...) {
            ## default non-generic arguments
            paramNames <- c("scores.only", "pop",
                            "summaryFun", "quantized",
                            "ref", "alt", "minoverlap", "caching")
            pop <- defaultPopulation(x)
            scores.only <- FALSE
            summaryFun <- mean
            quantized <- FALSE
            ref <- character(0)
            alt <- character(0)
            minoverlap <- 1L ## only relevant for genomic scores associated with nonSNRs
            caching <- TRUE
            ids <- ranges

            ## get arguments
            arglist <- list(...)
            mask <- nchar(names(arglist)) == 0
            if (any(mask))
              names(arglist)[mask] <- paste0("X", 1:sum(mask))

            mask <- names(arglist) %in% paramNames
            if (any(!mask))
                stop(sprintf("unused argument (%s)", names(arglist)[!mask]))
            list2env(arglist, envir=sys.frame(sys.nframe()))

            mask <- pop %in% populations(x)
            if (any(!mask))
              stop(sprintf("scores population %s is not present in %s. Please use 'populations()' to find out the available ones.",
                           pop[!mask], name(x)))

            ans <- DataFrame(as.data.frame(matrix(NA_real_, nrow=length(ids),
                                                  ncol=length(pop),
                                                  dimnames=list(NULL, pop))),
                             row.names=ids)
            if (!exists("rsIDs", envir=x@.data_cache)) {
              if (file.exists(file.path(x@data_dirpath, "rsIDs.rds"))) {
                message("Loading first time annotations of identifiers to genomic positions, produced by data provider.")
                rsIDs <- readRDS(file.path(x@data_dirpath, "rsIDs.rds"))
                assign("rsIDs", rsIDs, envir=x@.data_cache)
              } else {
                message("The data provider did not produce annotations of identifiers to genomic positions.")
                return(ans)
              }
            }

            rsIDs <- get("rsIDs", envir=x@.data_cache)
            stopifnot(is.integer(rsIDs)) ## QC
            idsint <- rep(NA_integer_, length(ids))
            rsMask <- regexpr("^rs", ids) == 1
            idsint[rsMask] <- as.integer(sub(pattern="^rs", replacement="", x=ids[rsMask]))
            mt <- rep(NA_integer_, length(idsint))
            if (any(!is.na(idsint))) {
              idsintnoNAs <- idsint[!is.na(idsint)]
              ord <- order(idsintnoNAs)                     ## order ids to speed up
              mtfi <- findInterval(idsintnoNAs[ord], rsIDs) ## call to findInterval()
              mtfi[ord] <- mtfi                             ## put matches into original order
              mt[!is.na(idsint)] <- mtfi                    ## integrate matches into result
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

            ranges <- GRanges()
            mcols(ranges) <- DataFrame(as.data.frame(matrix(NA_real_, nrow=0,
                                                            ncol=length(pop),
                                                            dimnames=list(NULL, pop))))
            if (any(!is.na(mt))) {
              if (!exists("rsIDgp", envir=x@.data_cache)) {
                rsIDgp <- readRDS(file.path(x@data_dirpath, "rsIDgp.rds"))
                assign("rsIDgp", rsIDgp, envir=x@.data_cache)
              }
              rsIDgp <- get("rsIDgp", envir=x@.data_cache)
              rng <- rsIDgp[mt[!is.na(mt)]]
              mask <- logical(length(mt))
              mask[!is.na(mt)] <- rng$isSNV
              if (any(mask))
                ans[mask, pop] <- gscores(x, rng[rng$isSNV], pop=pop, type="snrs",
                                          summaryFun=summaryFun, quantized=quantized,
                                          scores.only=TRUE, ref=ref, alt=alt,
                                          minoverlap=minoverlap, caching=caching)
              mask <- logical(length(mt))
              mask[!is.na(mt)] <- !rng$isSNV
              if (any(mask))
                ans[mask, pop] <- gscores(x, rng[!rng$isSNV], pop=pop, type="nonsnrs",
                                          summaryFun=summaryFun, quantized=quantized,
                                          scores.only=TRUE, ref=ref, alt=alt,
                                          minoverlap=minoverlap, caching=caching)
              ranges <- as(rng, "GRanges")
              names(ranges) <- rownames(ans)[!is.na(mt)]
              mcols(ranges) <- ans[!is.na(mt), , drop=FALSE]
            }

            if (scores.only)
              return(ans[!is.na(mt), ])

            ranges
          })

## fetch genomic scores stored on disk for specific sequences given in 'snames'
## and add them to the 'gsco1pop'data structure, which is then returned back
.fetch_scores_snrs <- function(object, objectname, gsco1pop, pop, snames) {
  slengths <- seqlengths(object)
  avh5snames <- h5fpath <- character(0)

  if (hdf5Backend(object)) {
    fname <- sprintf("%s.%s.h5", object@data_pkgname, pop)
    h5fpath <- file.path(object@data_dirpath, fname)
    h5str <- h5ls(h5fpath)
    avh5snames <- h5str$name[which(h5str$group == "/")]
  }

  for (sname in snames) {

    snameExists <- FALSE

    if (hdf5Backend(object)) { ## HDF5 backend

      if (sname %in% avh5snames) {
        snameExists <- TRUE
        gsco1pop[[sname]] <- HDF5Array(h5fpath, paste0(sname, "/scores"))
      }

    } else {                   ## RDS-serialized Rle backend

      if (pop == "default")
        fname <- sprintf("%s.%s.rds", object@data_pkgname, sname)
      else
        fname <- sprintf("%s.%s.%s.rds", object@data_pkgname, pop, sname)
      if (length(object@data_serialized_objnames) > 0 &&
          fname %in% names(object@data_serialized_objnames))
        fname <- object@data_serialized_objnames[fname]
      fpath <- file.path(object@data_dirpath, fname)
      if (file.exists(fpath)) {
        snameExists <- TRUE
        gsco1pop[[sname]] <- readRDS(fpath)
      }

    }

    if (!snameExists) {
      message(sprintf("no %s scores for population %s in sequence %s from %s object %s (%s).",
                       type(object), pop, sname, class(object), objectname, name(object)))
      gsco1pop[[sname]] <- Rle(lengths=slengths[sname], values=as.raw(0L))
    }

  }

  gsco1pop
}

.check_ref_alt_args <- function(ref, alt) {
  if (length(ref) != length(alt))
    stop("'ref' and 'alt' arguments have different lengths.")

  if (length(ref) > 0) {
    if (!class(ref) %in% c("character", "DNAStringSet", "DNAStringSetList")) {
      stop("'ref' argument must be either a character vector, a DNAStringSet or a DNAStringSetList object.")
    } else if (class(ref) == "DNAStringSetList") {
      mask <- elementNROWS(ref) > 1
      if (any(mask)) {
        ## stop("'ref' argument must contain only a single nucleotide per position.")
        ref <- unstrsplit(CharacterList(ref), sep=",")
        ref[mask] <- NA_character_
      }
      ref <- unlist(ref)
    }

    if (!class(alt) %in% c("character", "DNAStringSet", "DNAStringSetList")) {
      stop("'alt' argument must be either a character vector, a DNAStringSet or a DNAStringSetList object.")
    } else if (class(alt) == "DNAStringSetList") {
      mask <- elementNROWS(alt) > 1
      if (any(mask)) {
        ## stop("'alt' argument must contain only a single nucleotide per position.")
        alt <- unstrsplit(CharacterList(alt), sep=",")
        alt[mask] <- NA_character_
      }
      alt <- unlist(alt)
    }
  }

  list(ref=ref, alt=alt)
}

## assumes 'object' and 'ranges' have the same sequence "styles"
.scores_snrs <- function(object, ranges, pop, summaryFun=mean, quantized=FALSE,
                         scores.only=FALSE, ref=character(0), alt=character(0),
                         caching=TRUE) {
  objectname <- deparse(substitute(object))
  if (length(ranges) == 0)
    return(numeric(0))

  ra <- .check_ref_alt_args(ref, alt)
  ref <- ra$ref
  alt <- ra$alt
  rm(ra)

  snames <- unique(as.character(runValue(seqnames(ranges))))
  if (any(!snames %in% seqnames(object)))
    stop(sprintf("Sequence names %s in 'ranges' not present in reference genome %s.",
                 paste(snames[!snames %in% seqnames(object)], collapse=", "),
                 providerVersion(genomeDescription(object))))

  gscopops <- get(object@data_pkgname, envir=object@.data_cache)
  if (all(names(gscopops) %in% seqlevels(object))) { ## temporary fix for annotations w/o populations
    tmp <- gscopops
    gscopops <- list()
    gscopops[[defaultPopulation(object)]] <- List() ## RleList(tmp, compress=FALSE)
    assign(object@data_pkgname, gscopops, envir=object@.data_cache)
  }

  missingMask <- !pop %in% names(gscopops)
  for (popname in pop[missingMask])
    gscopops[[popname]] <- List() ## RleList(compress=FALSE)
  anyMissing <- any(missingMask)

  ans <- DataFrame(as.data.frame(matrix(NA_real_, nrow=length(ranges),
                                        ncol=length(pop), dimnames=list(NULL, pop))))
  ans2 <- NULL

  for (popname in pop) {
    missingMask <- !snames %in% names(gscopops[[popname]])
    anyMissing <- anyMissing || any(missingMask)
    if (any(missingMask))
      gscopops[[popname]] <- .fetch_scores_snrs(object, objectname, gscopops[[popname]],
                                                popname, snames[missingMask])

    sco <- rep(NA_real_, length(ranges))
    mdata <- metadata(gscopops[[popname]]) ## metadata in List object of HDF5 pkgs
    if ("Rle" %in% class(gscopops[[popname]][[1]]))
      mdata <- metadata(gscopops[[popname]][[1]])
    if (length(mdata) > 0 && class(mdata$dqfun) == "function")
      sco <- .rleGetValues(gscopops[[popname]], ranges, mdata,
                           summaryFun=summaryFun, quantized=quantized,
                           isHDF5=hdf5Backend(object))
    else
      warning(sprintf("No dequantization function available for scores population %s. Scores will be all NA for this population.",
                      popname))

  if (length(ref) > 0 && length(alt) > 0) {
    if (is.matrix(sco)) { ## if it's a matrix, then we have multiple scores per position
      if (any(is.na(ref)) || any(is.na(alt)))
        warning("'ref' or 'alt' arguments contain more than one allele per position, NAs introduced.")
      mt.r <- match(ref, DNA_BASES)
      mt.a <- match(alt, DNA_BASES)
      mask <- mt.r < mt.a
      idxcol <- mt.a
      idxcol[mask] <- mt.a[mask] - 1L
      sco <- sco[cbind(1:nrow(sco), idxcol)]
    } else { # if it's not a matrix, then assume we have metadata associated with the scores
      maskREF <- as.logical(.rleGetMetaValues(gscopops[[popname]],
                                              ranges, "maskREF",
                                              isHDF5=hdf5Backend(object)))
      if (is.null(ans2) && length(maskREF) > 0) {
        ans2 <- ans
        colnames(ans) <- paste0(colnames(ans), "_REF")
        colnames(ans2) <- paste0(colnames(ans2), "_ALT")
      }
      popnameALT <- paste0(popname, "_ALT")
      popname <- paste0(popname, "_REF")
      ans2[[popnameALT]] <- sco
      ans2[[popnameALT]][maskREF] <- 1 - sco[maskREF]
      sco[!maskREF] <- 1 - sco[!maskREF]
    }
  }

  ans[[popname]] <- sco
  }
  if (anyMissing && caching)
    assign(object@data_pkgname, gscopops, envir=object@.data_cache)
  rm(gscopops)

  if (!is.null(ans2)) {
    ans <- cbind(ans, ans2)
    rm(ans2)
  }

  if (scores.only)
    return(ans)

  mcols(ranges) <- cbind(mcols(ranges), ans)
  seqlevels(ranges) <- seqlevels(object)
  seqinfo(ranges) <- seqinfo(object)

  ranges
}

## fetch genomic scores stored on disk for specific sequences given in 'snames'
## and add them to the 'gscononsnrs'data structure, which is then returned back
.fetch_scores_nonsnrs <- function(object, objectname, gscononsnrs, pop, snames) {
  popnames <- character(0)
  if (length(gscononsnrs) > 0)
    popnames <- colnames(mcols(gscononsnrs[[1]]))

  if (length(snames) > 0) {
    md <- list()
    slengths <- seqlengths(object)
    for (sname in snames) {
      fname <- file.path(object@data_dirpath,
                         sprintf("%s.GRnonsnv.%s.rds", name(object), sname))
      obj <- GRanges()
      if (file.exists(fname)) {
        obj <- readRDS(fname)
        for (popname in popnames) {
          fname <- file.path(object@data_dirpath,
                             sprintf("%s.RLEnonsnv.%s.%s.rds", name(object), popname, sname))
          if (file.exists(fname)) {
            rleobj <- readRDS(fname)
            md <- metadata(rleobj)
            mcols(obj)[[popname]] <- rleobj
          } else {
            message(sprintf("no %s scores for population %s nonSNRs in sequence %s from %s object %s.",
                            name(object), popname, sname, class(object), objectname))
            ## mcols(obj)[[popname]] <- rep(NA_real_, length(obj))
            mcols(obj)[[popname]] <- Rle(lengths=length(obj), values=as.raw(0L)) ## should this be NA?
          }
        }
      } else {
        message(sprintf("no %s scores for nonSNRs in sequence %s from %s object %s.",
                        name(object), sname, class(object), objectname))
        mcols(obj) <- DataFrame(as.data.frame(matrix(raw(0), nrow=0, ncol=length(popnames),
                                                     dimnames=list(NULL, popnames))))
      }
      gscononsnrs[[sname]] <- obj
    }
    metadata(gscononsnrs) <- md ## is this really necessary??
  }

  if (length(pop) > 0) {
    md <- list()
    tmp <- unlist(gscononsnrs)
    names(tmp) <- NULL
    for (popname in pop) {
      scovalues <- Rle(raw(length(tmp)))
      i <- 1
      ## b/c we're storing MAF values as metadata columns of GRanges
      ## populations need to be loaded for all already loaded chromosomes
      for (sname in names(gscononsnrs)) {
        fname <- file.path(object@data_dirpath,
                           sprintf("%s.RLEnonsnv.%s.%s.rds", name(object), popname, sname))
        if (file.exists(fname)) {
          obj <- readRDS(fname)
          md <- metadata(obj)
          scovalues[i:(i+length(gscononsnrs[[sname]])-1)] <- obj
        }
        i <- i + length(gscononsnrs)
      }
      mcols(tmp)[[popname]] <- scovalues
      ## add maskREF metadata to each population, if available
      metadata(mcols(tmp)[[popname]])$maskREF <- md$maskREF
    }
    gscononsnrs <- split(tmp, seqnames(tmp), drop=TRUE)
    metadata(gscononsnrs) <- md
    rm(tmp)
  }

  gscononsnrs
}

.scores_nonsnrs <- function(object, ranges, pop, quantized=FALSE,
                            scores.only=FALSE, ref=character(0), alt=character(0),
                            minoverlap=1L, caching=TRUE) {
  objectname <- deparse(substitute(object))
  gscononsnrs<- get(sprintf("%s.nonsnvs", name(object)), envir=object@.data_cache)

  ra <- .check_ref_alt_args(ref, alt)
  ref <- ra$ref
  alt <- ra$alt
  rm(ra)

  ## ## by now, only multiple SNR alleles considered
  ## if (length(ref) > 0 && length(alt) > 0)
  ##   message("arguments 'ref' and 'alt' are given but there is only one score per genomic position.")

  if (length(intersect(seqlevelsStyle(ranges), seqlevelsStyle(object))) == 0)
    seqlevelsStyle(ranges) <- seqlevelsStyle(object)[1]

  snames <- unique(as.character(runValue(seqnames(ranges))))
  if (any(!snames %in% seqnames(object)))
    stop(sprintf("Sequence names %s in 'ranges' not present in reference genome %s.",
                 paste(snames[!snames %in% seqnames(object)], collapse=", "),
                 providerVersion(genomeDescription(object))))

  missingSeqMask <- !snames %in% names(gscononsnrs)
  anyMissing <- any(missingSeqMask)
  allpopnames <- character(0)
  if (length(gscononsnrs) > 0)
    allpopnames <- colnames(mcols(gscononsnrs[[1]]))
  missingPopMask <- !pop %in% allpopnames
  anyMissing <- anyMissing || any(missingPopMask)

  if (any(missingSeqMask) || any(missingPopMask)) ## if genomic scores belong to queried sequences that are not cached, load them
    gscononsnrs <- .fetch_scores_nonsnrs(object, objectname, gscononsnrs,
                                         pop[missingPopMask], snames[missingSeqMask])

  ans <- ans2 <- NULL
  if (quantized)
    ans <- DataFrame(as.data.frame(matrix(raw(0), nrow=length(ranges), ncol=length(pop),
                                        dimnames=list(NULL, pop))))
  else
    ans <- DataFrame(as.data.frame(matrix(NA_real_, nrow=length(ranges), ncol=length(pop),
                                        dimnames=list(NULL, pop))))

  ## the default value of 'minoverlap=1L' assumes that the sought nonsnrs are
  ## stored as in VCF files, using the nucleotide composition of the reference sequence
  ov <- findOverlaps(ranges, unlist(gscononsnrs),
                     minoverlap=minoverlap, type="equal")
  qHits <- queryHits(ov)
  sHits <- subjectHits(ov)
  if (length(ov) > 0) {
    q <- mcols(unlist(gscononsnrs))[sHits, pop, drop=FALSE]
    if (quantized)
      ans[qHits, pop] <- DataFrame(lapply(q, decode))
    else {
      .dequantizer <- metadata(gscononsnrs)$dqfun
      dqargs <- metadata(gscononsnrs)$dqfun_args
      lappargs <- c(list(X=q, FUN=.dequantizer), dqargs)
      ans[qHits, pop] <- DataFrame(do.call("lapply", lappargs))
    }
    if (any(duplicated(qHits)))
      message("gscores: more than one genomic score overlapping queried positions, reporting only the first hit.")

    if (length(ref) > 0 && length(alt) > 0) {
      for (popname in pop) {
        maskREF <- as.logical(metadata(mcols(unlist(gscononsnrs))[[popname]])$maskREF[sHits])
        if (is.null(ans2) && length(maskREF) > 0) {
          ans2 <- ans
          colnames(ans) <- paste0(colnames(ans), "_REF")
          colnames(ans2) <- paste0(colnames(ans2), "_ALT")
        }
        popnameALT <- paste0(popname, "_ALT")
        popnameREF <- paste0(popname, "_REF")
        ans2[qHits, popnameALT] <- ans[qHits, popnameREF]
        ans2[qHits, popnameALT][maskREF] <- 1 - ans[qHits, popnameREF][maskREF]
        ans[qHits, popnameREF][!maskREF] <- 1 - ans[qHits, popnameREF][!maskREF]
      }
    }
  }

  if (anyMissing && caching)
    assign(sprintf("%s.nonsnvs", name(object)), gscononsnrs, envir=object@.data_cache)
  rm(gscononsnrs)

  if (!is.null(ans2)) {
    ans <- cbind(ans, ans2)
    rm(ans2)
  }

  if (scores.only)
    return(ans)

  mcols(ranges) <- cbind(mcols(ranges), ans)
  seqlevels(ranges) <- seqlevels(object)
  seqinfo(ranges) <- seqinfo(object)

  ranges
}

## getters qfun and dqfun
setMethod("qfun", "GScores",
          function(object, pop=defaultPopulation(object)) {
            if (!pop %in% populations(object))
              stop(sprintf("There is no score population %s in this GScores object. Use populations() to find out what score populations are available.", pop))
            obj <- get(object@data_pkgname, envir=object@.data_cache)
            fun <- NA
            if (hdf5Backend(object))
              fun <- metadata(obj[[pop]])$qfun
            else
              fun <- metadata(obj[[pop]][[1]])$qfun
            fun
          })

setMethod("dqfun", "GScores",
          function(object, pop=defaultPopulation(object)) {
            if (!pop %in% populations(object))
              stop(sprintf("There is no score population %s in this GScores object. Use populations() to find out what score populations are available.", pop))
            obj <- get(object@data_pkgname, envir=object@.data_cache)
            fun <- NA
            if (hdf5Backend(object))
              fun <- metadata(obj[[pop]])$dqfun
            else
              fun <- metadata(obj[[pop]][[1]])$dqfun
          })

citation <- function(package, ...) UseMethod("citation")
citation.character <- function(package, ...) {
  if (missing(package)) package <- "base"
  utils::citation(package, ...)
}
setMethod("citation", signature="missing", citation.character)
setMethod("citation", signature="character", citation.character)
citation.GScores <- function(package, ...) {
  obj <- get(package@data_pkgname, envir=package@.data_cache)
  cit <- bibentry()
  if ("Rle" %in% class(obj[[1]]) || length(obj[[1]]) == 0 || hdf5Backend(package))
    cit <- metadata(obj[[1]])$citation
  else
    cit <- metadata(obj[[1]][[1]])$citation
  if (is.null(cit))
    cit <- bibentry()
  cit
}
setMethod("citation", signature="GScores", citation.GScores)

.pprintseqs <- function(x) {
  y <- x
  if (length(x) > 5)
    y <- c(y[1:2], "...", y[length(y)])
  y <- paste(y, collapse=", ")
  y
}

.pprintnsites <- function(x) {
  y <- x
  if (x > 10e6)
    y <- sprintf("%.0f millions", x/1e6)
  else if (x > 1e6)
    y <- sprintf("%.1f millions", x/1e6)
  y
}

## show method
setMethod("show", "GScores",
          function(object) {
            snrobj <- get(name(object), envir=object@.data_cache) ## single-nucleotide ranges
            if (all(names(snrobj) %in% seqlevels(object))) { ## temporary fix will working w/ outdated annotations
              tmp <- snrobj
              snrobj <- list()
              snrobj[[defaultPopulation(object)]] <- RleList(tmp, compress=FALSE)
              assign(object@data_pkgname, snrobj, envir=object@.data_cache)
            }

            nonsnrobj <- get(paste0(object@data_pkgname, ".nonsnvs"),
                             envir=object@.data_cache)
            loadedsnrpops <- loadedsnrseqs <- "none"
            loadednonsnrpops <- loadednonsnrseqs <- "none"
            if (length(snrobj) > 0) {
              loadedsnrseqs <- names(snrobj)
              if (length(populations(object)) > 1) {
                loadedsnrseqs <- names(snrobj[[1]])
                loadedsnrpops <- names(snrobj)
              }
            }
            if (ncol(mcols(nonsnrobj)) > 0)
              loadednonsnrpops <- colnames(mcols(nonsnrobj))
            if (length(nonsnrobj) > 0)
              loadednonsnrseqs <- unique(seqnames(nonsnrobj))
              
            max.abs.error <- NA
            if (length(length(loadedsnrseqs)) > 0) {
              if (!hdf5Backend(object))
                max.abs.error <- max(unlist(sapply(lapply(snrobj[[defaultPopulation(object)]],
                                                          metadata), "[[", "max_abs_error"),
                                            use.names = FALSE))
              else
                max.abs.error <- max(unlist(lapply(snrobj[[defaultPopulation(object)]],
                                                   function(x)
                                                     as.numeric(HDF5Array(path(seed(x)),
                                                               name=sprintf("%s/max_abs_error",
                                                                            sub("/scores", "", seed(x)@name))))),
                                            use.names=FALSE))
            }

            cat(class(object), " object \n",
                "# organism: ", organism(object), " (", provider(genomeDescription(object)), ", ",
                                providerVersion(genomeDescription(object)), ")\n",
                "# provider: ", provider(object), "\n",
                "# provider version: ", providerVersion(object), "\n",
                "# download date: ", object@download_date, "\n", sep="")
            if (!hdf5Backend(object)) {
              if (gscoresNonSNRs(object)) {
                cat("# loaded sequences (SNRs): ", .pprintseqs(loadedsnrseqs), "\n",
                    "# loaded sequences (nonSNRs): ", .pprintseqs(loadednonsnrseqs), "\n", sep="")
                if (loadedsnrpops[1] != "none" && length(loadedsnrpops) > 1)
                  cat("# loaded populations (SNRs): ", .pprintseqs(loadedsnrpops), "\n",
                      "# loaded populations (nonSNRs): ", .pprintseqs(loadednonsnrpops), "\n", sep="")
              } else {
                cat("# loaded sequences: ", .pprintseqs(loadedsnrseqs), "\n", sep="")
                if (loadedsnrpops[1] != "none" && length(loadedsnrpops) > 1)
                  cat("# loaded populations: ", .pprintseqs(loadedsnrpops), "\n", sep="")
              }
            }
            if (defaultPopulation(object) != "default")
              cat("# default scores population: ", defaultPopulation(object), "\n", sep="")
            if (!is.na(nsites(object)) && nsites(object) > 0)
              cat("# number of sites: ", .pprintnsites(nsites(object)), "\n", sep="")
            if (!is.na(max.abs.error)) {
              if (defaultPopulation(object) != "default")
                cat("# maximum abs. error (def. pop.): ", signif(max.abs.error, 3), "\n", sep="")
              else
                cat("# maximum abs. error: ", signif(max.abs.error, 3), "\n", sep="")
            }
            if (length(citation(object)) > 0)
              cat("# use 'citation()' to cite these data in publications\n")
          })

.rgscoresRle <- function(object, gscopops, snames, n, pop, ranges) {
  nscores <- sapply(gscopops[[pop]][snames], function(x) {
                      mask <- runValue(x) != 0
                      sum(runLength(x)[mask])
                    })
  f <- factor(sample(snames, size=n, replace=TRUE, prob=nscores/sum(nscores)),
              levels=snames)
  nxs <- table(f)
  dqf <- dqfun(object)
  pxs <- lapply(setNames(snames, snames), function(s, l, n) {
                  gr <- GRanges(seqinfo=seqinfo(object))
                  if (n[s] > 0) {
                    mask <- runValue(l[[s]]) != 0
                    whmask <- which(mask)
                    rlcsum <- cumsum(runLength(l[[s]])[whmask])
                    rndpos <- sample(rlcsum[length(rlcsum)], size=n[s], replace=FALSE)
                    rndpos <- sapply(rndpos, function(x) sum(rlcsum <= x))
                    rndpos <- sapply(whmask[rndpos], function(j) sum(runLength(l[[s]])[1:j]))
                    rndsco <- l[[s]][rndpos]
                    gr <- GRanges(seqnames=rep(s, n[s]),
                                  ranges=IRanges(start=rndpos, width=1),
                                  seqinfo=seqinfo(object))
                    mcols(gr)[[pop]] <- dqf(rndsco)
                  }
                  gr
                }, gscopops[[pop]], nxs)

  unlist(do.call("GRangesList", pxs))
}

.rgscoresHDF5 <- function(object, gscopops, snames, n, pop, ranges) {
  nscores <- sapply(gscopops[[pop]][snames], function(x) {
                      mask <- x != 0
                      sum(mask)
                    })
  f <- factor(sample(snames, size=n, replace=TRUE, prob=nscores/sum(nscores)),
              levels=snames)
  nxs <- table(f)
  dqf <- dqfun(object)
  pxs <- lapply(setNames(snames, snames), function(s, l, n) {
                  gr <- GRanges(seqinfo=seqinfo(object))
                  if (n[s] > 0) {
                    mask <- l[[s]] != 0
                    whmask <- which(mask)
                    rndpos <- sample(whmask, size=n[s], replace=FALSE)
                    rndsco <- l[[s]][rndpos]
                    gr <- GRanges(seqnames=rep(s, n[s]),
                                  ranges=IRanges(start=rndpos, width=1),
                                  seqinfo=seqinfo(object))
                    mcols(gr)[[pop]] <- dqf(rndsco)
                  }
                  gr
                }, gscopops[[pop]], nxs)

  unlist(do.call("GRangesList", pxs))
}

setMethod("rgscores", signature(n="GScores", object="missing"),
          function(n, object, ...) {
            rgscores(n=1L, object=n, ...)
          })

setMethod("rgscores", signature(n="missing", object="GScores"),
          function(n, object, ...) {
            rgscores(n=1L, object, ...)
          })

setMethod("rgscores", signature(n="numeric", object="GScores"),
          function(n=1, object, ...) {
            rgscores(n=as.integer(n), object, ...)
          })

setMethod("rgscores", signature(n="integer", object="GScores"),
          function(n=1L, object, ...) {
            ## default non-generic arguments
            paramNames <- c("pop", "ranges", "scores.only")
            pop <- defaultPopulation(object)
            ranges <- keepStandardChromosomes(GRanges(seqinfo=seqinfo(object)))
            scores.only <- FALSE

            ## get arguments
            arglist <- list(...)
            mask <- nchar(names(arglist)) == 0
            if (any(mask))
              names(arglist)[mask] <- paste0("X", 1:sum(mask))

            mask <- names(arglist) %in% paramNames
            if (any(!mask))
                stop(sprintf("unused argument (%s)", names(arglist)[!mask]))
            list2env(arglist, envir=sys.frame(sys.nframe()))

            objectname <- deparse(substitute(object))
            stopifnot(length(pop) == 1) ## QC
            gscopops <- get(object@data_pkgname, envir=object@.data_cache)

            if (!is.logical(scores.only) || length(scores.only) > 1)
              stop("'scores.only' must be one logical value (TRUE or FALSE).")

            if (!is(ranges, "GRanges") && !is.character(ranges))
              stop("'ranges' must be either a 'GRanges' object or a character vector of sequence names.")
            if (is.character(ranges)) {
              mask <- !ranges %in% seqlevels(object)
              if (any(mask))
                stop(sprintf("Sequence names %s not in %s.",
                             paste(ranges[mask], collapse=", "), objectname))
              slen <- seqlengths(object)[ranges]
              ranges <- GRanges(seqnames=ranges, IRanges(rep(1, length(ranges)), slen))
            }

            snames <- seqlevels(keepStandardChromosomes(seqinfo(object)))
            if (length(ranges) > 0) {
              snames <- unique(as.character(runValue(seqnames(ranges))))
              if (any(!snames %in% seqnames(object)))
                stop(sprintf("Sequence names %s in 'ranges' not present in reference genome %s.",
                             paste(snames[!snames %in% seqnames(object)], collapse=", "),
                             providerVersion(genomeDescription(object))))
            }

            if (!pop %in% populations(object))
              stop(sprintf("Population %s is not in 'GScores' object %s. Use 'populations()' to find out the available ones.",
                           pop, objectname))

            if (!pop %in% names(gscopops))
              gscopops[[pop]] <- List() ## RleList(compress=FALSE)

            missingMask <- !snames %in% names(gscopops[[pop]])
            if (any(missingMask)) {
              gscopops[[pop]] <- .fetch_scores_snrs(object, objectname, gscopops[[pop]],
                                                    pop, snames[missingMask])
            }

            gsco <- GRanges(seqinfo=seqinfo(object))
            if (hdf5Backend(object))
              gsco <- .rgscoresHDF5(object, gscopops, snames, n, pop, ranges)
            else
              gsco <- .rgscoresRle(object, gscopops, snames, n, pop, ranges)

            if (scores.only)
              gsco <- mcols(gsco)[[pop]]

            gsco
          })
