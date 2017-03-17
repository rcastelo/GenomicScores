## This file explains the steps to download and freeze the fitCons
## scores version 1.01 for human genome version hg19. If you use
## these data on your own research please cite the following publication:

## Gulko B, Gronau I, Hubisz MJ, Siepel A. Probabilities of fitness consequences
## for point mutations across the human genome. Nat. Genet. 2015 Aug;47:276-83.
## (http://www.nature.com/ng/journal/v47/n3/full/ng.3196.html)

## The data were downloaded from the CSHL mirror of the UCSC genome browser as follows
##
## wget http://compgen.cshl.edu/fitCons/0downloads/tracks/V1.01/i6/scores/fc-i6-0.bw
##
## further information about these data can be found in the README file at
##
## http://compgen.cshl.edu/fitCons/0downloads/tracks/V1.01/readme.txt
##
## and at the CSHL mirror of the UCSC genome browser at
##
## http://genome-mirror.cshl.edu

## The following R script processes the downloaded data to
## store the fitCons scores in raw-Rle objects

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(doParallel)
library(S4Vectors)

downloadURL <- "http://compgen.cshl.edu/fitCons/0downloads/tracks/V1.01/i6/scores/fc-i6-0.bw"

registerDoParallel(cores=4)

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

## transform BIGWIG into Rle objects coercing fitCons scores into
## 2-decimal digit raw-encoded values to reduce memory requirements
## in principle centiles of fitCons probabilities should give the
## necessary resolution for the purpose of filtering genetic variants
## with fitCons scores

## quantizer function. it maps input real-valued [0, 1]
## phastCons scores to non-negative integers [0, 255] so that
## each of them can be later coerced into a single byte (raw type).
## quantization is done by rounding to one decimal significant digit,
## and therefore, mapping is restricted to 101 different positive
## integers only [1-101], where the 0 value is kept to code for missingness
.quantizer <- function(x) {
  q <- as.integer(sprintf("%.0f", 100*x))
  q <- q + 1L
  q
}
attr(.quantizer, "description") <- "multiply by 100, round to nearest integer, add up one"

## dequantizer function. it maps input non-negative integers [0, 255]
## to real-valued phastCons scores, where 0 codes for NA
.dequantizer <- function(q) {
  x <- as.numeric(q)
  x[x == 0] <- NA
  x <- (x - 1) / 100
  x
}
attr(.dequantizer, "description") <- "subtract one integer unit, divide by 100"

foreach (chr=seqnames(Hsapiens)) %dopar% {
  cat(chr, "\n")
  tryCatch({
    rawscores <- import.bw(BigWigFile("fc-i6-0.bw"),
                           which=GRanges(seqnames=chr, IRanges(1, seqlengths(Hsapiens)[chr])))
    qscores <- .quantizer(rawscores$score)
    max.abs.error <- max(abs(rawscores$score - .dequantizer(qscores)))
    assign(sprintf("fitCons_%s", chr), coverage(rawscores, weight=qscores)[[chr]])
    assign(sprintf("fitCons_%s", chr),
           do.call("runValue<-", list(get(sprintf("fitCons_%s", chr)),
                                      as.raw(runValue(get(sprintf("fitCons_%s", chr)))))))
    obj <- get(sprintf("fitCons_%s", chr))
    Fn <- function(x) { warning("no ecdf() function available") ; numeric(0) }
    n <- length(unique(rawscores$score[!is.na(rawscores$score)]))
    if (n > 10) {
      if (n <= 10000) {
        Fn <- ecdf(rawscores$score)
      } else { ## to save space with more than 10,000 different values use sampling
        Fn <- ecdf(sample(rawscores$score[!is.na(rawscores$score)], size=10000, replace=TRUE))
      }
    }
    metadata(obj) <- list(seqname=chr,
                          provider="UCSC",
                          provider_version="21Aug2014", ## it'd better to grab the date from downloaded file
                          download_url=downloadURL,
                          download_date=format(Sys.Date(), "%b %d, %Y"),
                          reference_genome=refgenomeGD,
                          data_pkgname="fitCons.UCSC.hg38",
                          qfun=.quantizer,
                          dqfun=.dequantizer,
                          ecdf=Fn,
                          max_abs_error=max.abs.error)
    saveRDS(obj, file=sprintf("fitCons.UCSC.hg19.%s.rds", chr))
    rm(rawscores, obj)
    rm(list=sprintf("fitCons_%s", chr))
    gc()
  }, error=function(err) {
    message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
