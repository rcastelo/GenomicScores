## This file explains the steps to download and freeze the MCAP scores
## for the human genome version hg19. If you use these data on your own
## research please cite the following publication:

## Jagadeesh K, Wenger A, Berger M, Guturu H, Stenson P, Cooper D,
## Bernstein J and Bejerano G. M-CAP eliminates a majority of variants
## with uncertain significance in clinical exomes at high sensitivity.
## Nature Genetics, 2016. DOI: 10.1038/ng.3703

## The data were downloaded, uncompressed, compressed again using
## bgzip and tabix indexed
##
## $ wget http://bejerano.stanford.edu/MCAP/downloads/dat/mcap_v1_0.txt.gz
## $ gzip -d mcap_v1_0.txt.gz
## $ bgzip -c mcap_v1_0.txt > mcap_v1.0.txt.gz
## $ tabix -p vcf mcap_v1_0.txt.gz

## The following R script processes the downloaded data to
## store the MCAP scores in integer-Rle objects

library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(doParallel)

downloadURL <- "http://bejerano.stanford.edu/MCAP/downloads/dat/mcap_v1_0.txt.gz"

registerDoParallel(cores=2)

## freeze the GenomeDescription data for Hsapiens

refgenomeGD <- GenomeDescription(organism=organism(Hsapiens),
                                 common_name=commonName(Hsapiens),
                                 provider=provider(Hsapiens),
                                 provider_version=providerVersion(Hsapiens),
                                 release_date=releaseDate(Hsapiens),
                                 release_name=releaseName(Hsapiens),
                                 seqinfo=Seqinfo(seqnames=seqnames(Hsapiens),
                                                 seqlengths=seqlengths(Hsapiens),
                                                 isCircular=isCircular(Hsapiens),
                                                 genome=releaseName(Hsapiens)))

saveRDS(refgenomeGD, file="refgenomeGD.rds")

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

## quantizer function
## x: values to quantize, x >= 0 & x <= 1, length(x) is multiple of d
## n: maximum number of quantized values
## d: number of digits in 'x' forming a value to quantize
.quantizer <- function(x, n=Inf, d=1L, na.zero=FALSE) {
  stopifnot(d > 0L) ## QC
  stopifnot(length(x) %% d == 0) ## QC
  q <- rep(NA_integer_, length(x))
  q[!is.na(x)] <- as.integer(sprintf("%.0f", 100*x[!is.na(x)]))
  if (max(q, na.rm=TRUE) > (n - 1))
    q[q > (n - 1)] <- n - 1
  base <- n
  if (na.zero) { ## should we recode NAs into 0s?
    q <- q + 1L
    q[is.na(q)] <- 0L
    base <- base + 1L
  }
  if (d > 1L)
    q <- .toBase(d=q, g=rep(1:(length(x)/d), each=d), b=base)
  q <- q + 1L
  if (any(q < 0 || q > .Machine$integer.max))
    stop("current number of quantized values cannot be stored into one integer")
  q
}
attr(.quantizer, "description") <- "round to 2 digits, transform into integer, store them using a given base"

## dequantizer function
.dequantizer <- function(q, d, b, na.zero=FALSE) {
  q <- q - 1L
  if (d > 1L)
    q <- .fromBase(q, d=d, b=b)
  if (na.zero) { ## should we decode 0s as NAs
    q[q == 0L] <- NA
    q <- q - 1L
  }
  q <- q/100L
  q
}
attr(.dequantizer, "description") <- "transform into base 10, divide by 100"

tbx <- open(TabixFile("mcap_v1_0.txt.gz"))
tbxchr <- sortSeqlevels(seqnamesTabix(tbx))
tbxchrUCSC <- tbxchr
seqlevelsStyle(tbxchrUCSC) <- "UCSC"
close(tbx)

tbxchrUCSC <- tbxchrUCSC[tbxchrUCSC %in% seqlevels(Hsapiens)]
stopifnot(length(tbxchrUCSC) > 0) ## QC
allchrgr <- GRanges(seqnames=seqnames(Hsapiens)[match(tbxchrUCSC, seqnames(Hsapiens))],
                    IRanges(1, seqlengths(Hsapiens)[tbxchrUCSC]),
                    seqinfo=seqinfo(Hsapiens)[tbxchrUCSC])
names(allchrgr) <- seqnames(allchrgr)
seqlevelsStyle(allchrgr) <- seqlevelsStyle(tbxchr)[1]

foreach (chr=seqlevels(allchrgr)) %dopar% {
  chr <- names(allchrgr)[which(seqlevels(allchrgr) == chr)]
  cat(chr, "\n")
  tryCatch({
    tbx <- open(TabixFile("mcap_v1_0.txt.gz"))
    rawscores <- scanTabix(tbx, param=allchrgr[chr])
    rawscores <- rawscores[[1]]
    rawscores <- do.call("rbind", strsplit(rawscores, split="\t"))
    rawscores <- data.frame(POSITION=as.integer(rawscores[, 2]),
                            REF=rawscores[, 3],
                            ALT=rawscores[, 4],
                            SCORE=as.numeric(rawscores[ ,5]),
                            stringsAsFactors=FALSE)
    gc()
    ## b/c there are no MCAP scores for every possible nucleotide change
    ## we need to add those missing changes and set those scores to NA
    uniqpos <- unique(rawscores$POSITION)
    uniqref <- rawscores$REF[!duplicated(rawscores$POSITION)]
    allrawscores <- data.frame(POSITION=rep(uniqpos, each=4),
                               REF=rep(uniqref, each=4),
                               ALT=rep(c("A", "C", "G", "T"), times=length(uniqpos)),
                               stringsAsFactors=FALSE)
    allrawscores <- allrawscores[allrawscores$REF != allrawscores$ALT, ]
    allrawscores$SCORE <- NA_real_
    rownames(allrawscores) <- paste(as.character(allrawscores[, 1]),
                                    allrawscores[, 2], allrawscores[, 3], sep="_")
    mt <- match(paste(as.character(rawscores$POSITION), rawscores$REF, rawscores$ALT, sep="_"),
                rownames(allrawscores))
    stopifnot(all(!is.na(mt))) ## QC
    allrawscores$SCORE[mt] <- rawscores$SCORE
    rawscores <- allrawscores
    rm(allrawscores)
    gc()

    q <- .quantizer(rawscores$SCORE, n=101L, d=3L, na.zero=TRUE)
    gr <- sprintf("%s:%s-%s", chr, uniqpos, uniqpos)
    gr <- GRanges(gr)
    seqinfo(gr) <- seqinfo(Hsapiens)[chr]
    obj <- coverage(gr, weight=q)[[1]]
    x <- .dequantizer(q, d=3L, b=102L, na.zero=TRUE)
    max.abs.error <- max(abs(rawscores$SCORE - x), na.rm=TRUE)
    rawscores <- rawscores$SCORE[!is.na(rawscores$SCORE)]
    gc()
    if (length(unique(rawscores)) <= 10000)
      Fn <- ecdf(rawscores)
    else ## to save space with more than 10,000 different values use sampling
      Fn <- ecdf(sample(rawscores, size=10000, replace=TRUE))
    metadata(obj) <- list(seqname=chr,
                          provider="Stanford",
                          provider_version="v1.0",
                          download_url=downloadURL,
                          download_date=format(Sys.Date(), "%b %d, %Y"),
                          reference_genome=refgenomeGD,
                          data_pkgname="mcap.v1.0.hg19",
                          qfun=.quantizer,
                          dqfun=.dequantizer,
                          ecdf=Fn,
                          max_abs_error=max.abs.error)
    saveRDS(obj, file=sprintf("mcap.hg19.%s.rds", chr))
    rm(rawscores)
    rm(obj)
    gc()
    close(tbx)
  }, error=function(err) {
    message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
