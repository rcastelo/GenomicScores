## GScores constructor
GScores <- function(provider, provider_version, download_url,
                    download_date, reference_genome,
                    data_pkgname, data_dirpath) {
  data_serialized_objnames <- data_pkgname
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

## this has been improved using RleViews as discussed in
## https://stat.ethz.ch/pipermail/bioconductor/2013-December/056409.html
rleGetValues <- function(rlelst, gr, summaryFun="mean",
                         coercionFun="as.numeric") {
  numericmean <- TRUE
  if (summaryFun != "mean" || coercionFun != "as.numeric")
    numericmean <- FALSE

  summaryFun <- match.fun(summaryFun)
  coercionFun <- match.fun(coercionFun)
  seqlevels(gr) <- names(rlelst)
  ord <- order(seqnames(gr))
  ans <- numeric(length(gr))
  startbyseq <- split(start(gr), seqnames(gr), drop=TRUE)
  ans[ord] <- unlist(lapply(rlelst[startbyseq], coercionFun), use.names=FALSE)
  whregions <- which(width(gr) > 1)
  if (length(whregions) > 0) { ## regions comprising more than one position
    tmpans <- NA_real_         ## need to be summarized
    if (numericmean) {
      rngbyseq <- split(gr[whregions], seqnames(gr)[whregions])
      tmpans <- lapply(names(rngbyseq),
                       function(sname) {
                         coercedrle <- rlelst[[sname]]
                         ## this coercion can take up to one
                         ## second for chromosome 1, we should
                         ## consider storing coerced version if
                         ## memory consumption is not an issue
                         runValue(coercedrle) <- as.numeric(runValue(coercedrle))
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
                    function(i, r, p, w)
                      summaryFun(coercionFun(r[p[i]:(p[i]+w[i]-1)])), r, p, w)

      tmpans <- unlist(mapply(f,
                              rlelst[names(startbyseq)], startbyseq, widthbyseq,
                              SIMPLIFY=FALSE),
                       use.names=FALSE)
    }

    ans[ord[whregions]] <- tmpans
  }
  ans
}

setMethod("scores", c("GScores", "GRanges"),
          function(object, gpos, summaryFun="mean", coercionFun="as.numeric",
                   caching=TRUE) {
            if (seqlevelsStyle(gpos) != seqlevelsStyle(object))
              seqlevelsStyle(gpos) <- seqlevelsStyle(object)

            snames <- unique(as.character(runValue(seqnames(gpos))))
            if (any(!snames %in% seqnames(object)))
              stop("Sequence names %s in GRanges object not present in GScores object.",
                   paste(snames[!snames %in% seqnames(object)], collapse=", "))

            scorlelist <- get(object@data_pkgname, envir=object@.data_cache)
            missingMask <- !snames %in% names(scorlelist)
            for (sname in snames[missingMask]) {
              fp <- file.path(object@data_dirpath,
                              sprintf("%s.%s.rds", object@data_pkgname, sname))
              if (file.exists(fp))
                scorlelist[[sname]] <- readRDS(fp)
              else {
                warning(sprintf("No %s scores for chromosome %s.",
                                object@data_pkgname, sname))
                scorlelist[[sname]] <- Rle(raw())
              }
            }

            if (any(missingMask) && caching)
              assign(object@data_pkgname, scorlelist, envir=object@.data_cache)

            sco <- rleGetValues(scorlelist, gpos, summaryFun=summaryFun,
                                coercionFun=coercionFun)
            sco[sco == 255L] <- NA_integer_
            sco <- sco / 10L
            rm(scorlelist)

            sco
          })

## show method
setMethod("show", "GScores",
          function(object) {
            cat(class(object), " object for ", organism(object), " (",
                provider(object), ")\n", sep="")
          })
