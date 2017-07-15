## GScores constructor
GScores <- function(provider, provider_version, download_url,
                    download_date, reference_genome,
                    data_pkgname, data_dirpath,
                    data_serialized_objnames=character(0)) {
  data_cache <- new.env(hash=TRUE, parent=emptyenv())

  assign(data_pkgname, RleList(compress=FALSE), envir=data_cache)

  new("GScores", provider=provider,
                 provider_version=provider_version,
                 download_url=download_url,
                 download_date=download_date,
                 reference_genome=reference_genome,
                 data_pkgname=data_pkgname,
                 data_dirpath=data_dirpath,
                 data_serialized_objnames=data_serialized_objnames,
                 .data_cache=data_cache)
}

## accessors
setMethod("name", "GScores", function(x) x@data_pkgname)

setMethod("type", "GScores", function(x) sub("\\..*$", "", name(x)))

setMethod("provider", "GScores", function(x) x@provider)

setMethod("providerVersion", "GScores", function(x) x@provider_version)

setMethod("referenceGenome", "GScores", function(x) x@reference_genome)

setMethod("organism", "GScores",
          function(object) organism(referenceGenome(object)))

setMethod("seqinfo", "GScores", function(x) seqinfo(referenceGenome(x)))

setMethod("seqnames", "GScores", function(x) seqnames(referenceGenome(x)))

setMethod("seqlengths", "GScores", function(x) seqlengths(referenceGenome(x)))

setMethod("seqlevelsStyle", "GScores",
          function(x) seqlevelsStyle(referenceGenome(x)))


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
.rleGetValues <- function(rlelst, gr, summaryFun) {
  summaryFun <- match.fun(summaryFun)
  numericmean <- TRUE
  if (!identical(summaryFun, mean))
    numericmean <- FALSE

  .dequantizer <- metadata(rlelst[[1]])$dqfun
  dqargs <- metadata(rlelst[[1]])$dqfun_args
  seqlevels(gr) <- names(rlelst)
  ord <- order(seqnames(gr))
  startbyseq <- split(start(gr), seqnames(gr), drop=TRUE)
  lappargs <- c(list(X=rlelst[startbyseq], FUN=.dequantizer), dqargs)
  x <- unlist(do.call("lapply", lappargs), use.names=FALSE)
  ans <- numeric(0)
  valxpos <- metadata(rlelst[[1]])$valxpos
  if (is.null(valxpos))
    valxpos <- 1
  stopifnot(length(valxpos) == 1 && valxpos[1] > 0) ## QC
  whregions <- which(width(gr) > 1)
  if (valxpos == 1) {
    ans <- numeric(length(gr))
    ans[ord] <- x
  } else {
    if (length(whregions) > 0)
      stop("This GScores object returns more than one value per position, and therefore, input ranges can only have width one.")
    ans <- matrix(NA_real_, nrow=length(gr), ncol=valxpos)
    ans[ord, ] <- matrix(x, nrow=length(gr), ncol=valxpos, byrow=TRUE)
  }

  if (length(whregions) > 0) { ## regions comprising more than one position
    tmpans <- NA_real_         ## need to be summarized
    if (numericmean) {
      rngbyseq <- split(gr[whregions], seqnames(gr)[whregions])
      tmpans <- lapply(names(rngbyseq),
                       function(sname) {
                         coercedrle <- rlelst[[sname]]
                         dqargs2 <- c(list(q=runValue(coercedrle)), dqargs)
                         runValue(coercedrle) <- do.call(".dequantizer", dqargs2)
                         viewMeans(Views(coercedrle,
                                         start=start(rngbyseq[[sname]]),
                                         end=end(rngbyseq[[sname]])))
                       })
      tmpans <- unsplit(tmpans, as.factor(seqnames(gr)[whregions]))
    } else { ## this allows for other summary and coercion functions
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

    ans[ord[whregions]] <- tmpans
  }
  ans
}

setMethod("scores", c("GScores", "GenomicRanges"),
          function(object, ranges, scores.only=FALSE, ...) {
            objectname <- deparse(substitute(object))
            ## default non-generic arguments
            summaryFun <- mean
            caching <- TRUE

            ## get arguments
            arglist <- list(...)
            mask <- nchar(names(arglist)) == 0
            if (any(mask))
              names(arglist)[mask] <- paste0("X", 1:sum(mask))
            list2env(arglist, envir=sys.frame(sys.nframe()))

            if (length(ranges) == 0)
              return(numeric(0))

            if (length(intersect(seqlevelsStyle(ranges), seqlevelsStyle(object))) == 0)
              seqlevelsStyle(ranges) <- seqlevelsStyle(object)[1]

            snames <- unique(as.character(runValue(seqnames(ranges))))
            if (any(!snames %in% seqnames(object)))
              stop(sprintf("Sequence names %s in GRanges object not present in reference genome %s.",
                           paste(snames[!snames %in% seqnames(object)], collapse=", "),
                           providerVersion(referenceGenome(object))))

            scorlelist <- get(object@data_pkgname, envir=object@.data_cache)
            missingMask <- !snames %in% names(scorlelist)
            slengths <- seqlengths(object)
            for (sname in snames[missingMask]) {
              fname <- sprintf("%s.%s.rds", object@data_pkgname, sname)
              if (length(object@data_serialized_objnames) > 0 &&
                  fname %in% names(object@data_serialized_objnames))
                fname <- object@data_serialized_objnames[fname]
              fpath <- file.path(object@data_dirpath, fname)
              if (file.exists(fpath))
                scorlelist[[sname]] <- readRDS(fpath)
              else {
                warning(sprintf("No %s scores for sequence %s in %s object '%s'.",
                                object@data_pkgname, sname, class(object),
                                objectname))
                scorlelist[[sname]] <- Rle(lengths=slengths[sname],
                                           values=as.raw(0L))
              }
            }

            if (any(missingMask) && caching)
              assign(object@data_pkgname, scorlelist, envir=object@.data_cache)

            sco <- .rleGetValues(scorlelist, ranges, summaryFun=summaryFun)
            rm(scorlelist)

            if (scores.only) {
              if (is.matrix(sco))
                colnames(sco) <- paste0("scores", 1:ncol(sco))
              return(sco)
            }

            if (is.matrix(sco)) {
              colnames(sco) <- paste0("scores", 1:ncol(sco))
              mcols(ranges) <- cbind(mcols(ranges), DataFrame(as.data.frame(sco)))
            } else
              ranges$scores <- sco

            ranges
          })

## getter qfun and dqfun methods
## setMethod("qfun", "GScores",
##           function(object) {
##             obj <- get(object@data_pkgname, envir=object@.data_cache)
##             metadata(obj[[1]])$qfun
##           })
## 
## setMethod("dqfun", "GScores",
##           function(object) {
##             obj <- get(object@data_pkgname, envir=object@.data_cache)
##             metadata(obj[[1]])$dqfun
##           })
citation <- function(package, ...) UseMethod("citation")
citation.character <- function(package, ...) {
  if (missing(package)) package <- "base"
  utils::citation(package, ...)
}
setMethod("citation", signature="missing", citation.character)
setMethod("citation", signature="character", citation.character)
citation.GScores <- function(package, ...) {
  obj <- get(package@data_pkgname, envir=package@.data_cache)
  cit <- metadata(obj[[1]])$citation
  if (is.null(cit))
    cit <- bibentry()
  cit
}
setMethod("citation", signature="GScores", citation.GScores)

## show method
setMethod("show", "GScores",
          function(object) {
            obj <- get(name(object), envir=object@.data_cache)
            seqs <- names(obj)
            max.abs.error <- NA
            if (length(seqs) > 0) {
              max.abs.error <- max(sapply(lapply(obj, metadata), "[[", "max_abs_error"))
            } else
              seqs <- "none"
            if (length(seqs) > 5)
              seqs <- c(seqs[1:2],  "...", seqs[length(seqs)])
            seqs <- paste(seqs, collapse=", ")
            cat(class(object), " object \n",
                "# organism: ", organism(referenceGenome(object)), " (", provider(referenceGenome(object)), ", ", providerVersion(referenceGenome(object)), ")\n",
                "# provider: ", provider(object), "\n",
                "# provider version: ", providerVersion(object), "\n",
                "# download date: ", object@download_date, "\n",
                "# loaded sequences: ", seqs, "\n",
                "# maximum abs. error: ", signif(max.abs.error, 3), "\n", sep="")
            if (length(citation(object)) > 0)
              cat("# use 'citation()' to know how to cite these data in publications\n")
          })
