## This file explains the steps to download and freeze the phyloP
## conservation scores for human genome version hg19, calculated on
## 100 vertebrate species. If you use these data on your own research
## please cite the following publication:

## Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. Detection of nonneutral
## substitution rates on mammalian phylogenies. Genome Res. 2010 Oct;20(1):110-21.
## (http://genome.cshlp.org/content/20/1/110)

## The data were downloaded from the UCSC genome browser with the Unix 'rsync'
## command as follows
##
## $ rsync -avz --progress \
##     rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way/ \
##     ./hg19.100way.phyloP

## The following R script processes the downloaded data to
## store the phyloP scores in raw-Rle objects

library(plyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(doParallel)

downloadURL <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way"

registerDoParallel(cores=4) ## each process may need up to 20Gb of RAM

## transform WIG to BIGWIG format
si <- Seqinfo(seqnames=seqnames(Hsapiens), seqlengths=seqlengths(Hsapiens))
foreach (chr=seqnames(Hsapiens)) %dopar% {
  cat(chr, "\n")
  wigToBigWig(file.path("hg19.100way.phyloP", sprintf("%s.phyloP100way.wigFix.gz", chr)), seqinfo=si)
}

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

saveRDS(refgenomeGD, file=file.path("hg19.100way.phyloP", "refgenomeGD.rds"))

## transform BIGWIG into Rle objects coercing phyloP scores into
## 1-decimal digit raw-encoded values to reduce memory requirements
## in principle deciles of phyloP probabilities should give the
## necessary resolution for the purpose of filtering genetic variants
## on conservation

## quantizer function. it maps input real-valued [-Inf, Inf]
## phyloP scores to non-negative integers [0, 255] so that
## each of them can be later coerced into a single byte (raw type).
## according to http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=cons100way
## phyloP scores were calculated using a likelihood ratio test (LRT) and
## according to Pollard et al. (2010), pg. 118, phyloP scores are derived from
## LRT p-values by transforming them as -log10(P) where P is a two-sided P-value
## and setting this score as positive or negative depending on whether the site
## was predicted to be conserved or fast-evolving (accelerated), respectively.
## quantization is done by taking the absolute value of the phyloP score,
## rounding to closest 0.5 those between 0 and 2, rounding to the
## closest integer those between 2 and 10, and setting to 10 those larger than 10.
## this quantization maps phyloP scores to 29 different positive integers, which
## take values [1-29], to keep the 0 value to code for missingness, and will code
## for zero and positive phyloP scores. negative phyloP scores will be mapped to integer
## values [30-57].
.quantizer <- function(x) {
  mask <- abs(x) < 2
  x[mask] <- round_any(x[mask], 0.5, round)
  q <- as.integer(sprintf("%.0f", 10*x))
  mask <- abs(q) > 100L
  q[mask] <- as.integer(sign(q[mask]))*100L
  mask <- abs(q) > 20L
  q[mask] <- as.integer(sprintf("%.0f", q[mask]/10)) + sign(q[mask])*18
  q[q < 0] <- abs(q[q < 0]) + 28L
  q <- q + 1L
  q
}
attr(.quantizer, "description") <- "abs(x) < 2 round to 0.5, abs(x) >= 2 round to closest integer, set to 10 or -10  values exceeding 10 or -10, respectively, add up one"

## dequantizer function. it maps input non-negative integers [0, 255]
## to real-valued phyloP scores, where 0 codes for NA
.dequantizer <- function(q) {
  x <- as.numeric(q)
  x[x == 0] <- NA
  x <- x - 1
  maskNAs <- is.na(x)
  mask <- !maskNAs & x > 28
  x[mask] <- -1*(x[mask]-28)
  mask <- !maskNAs & abs(x) <= 20
  x[mask] <- x[mask] / 10
  mask <- !maskNAs & abs(x) > 20
  x[mask] <- x[mask] - sign(x[mask])*18
  x
}
attr(.dequantizer, "description") <- "subtract one integer unit, set sign, abs(x) < 20 divide by 10, abs(x) > 20 subtract 18"

foreach (chr=seqnames(Hsapiens)) %dopar% {
  cat(chr, "\n")
  tryCatch({
    rawscores <- import.bw(BigWigFile(file.path("hg19.100way.phyloP", sprintf("%s.phyloP100way.bw", chr))))
    qscores <- .quantizer(rawscores$score)
    ## calculate maximum absolute error for absolute phyloP scores <= 10
    ## b/c those larger than 10 have been set to 10
    f <- cut(abs(rawscores$score), breaks=c(0, 2, max(abs(rawscores$score))), include.lowest=TRUE)
    err <- abs(rawscores$score - .dequantizer(qscores))
    mask10 <- abs(rawscores$score) <= 10
    max.abs.error <- tapply(err[mask10], f[mask10], max, na.rm=TRUE)
    assign(sprintf("phyloP100way_%s", chr), coverage(rawscores, weight=qscores)[[chr]])
    assign(sprintf("phyloP100way_%s", chr),
           do.call("runValue<-", list(get(sprintf("phyloP100way_%s", chr)),
                                      as.raw(runValue(get(sprintf("phyloP100way_%s", chr)))))))
    obj <- get(sprintf("phyloP100way_%s", chr))
    n <- length(unique(rawscores$score[!is.na(rawscores$score)]))
    Fn <- function(x) 0
    if (n > 10) {
      if (n <= 10000) {
        Fn <- ecdf(rawscores$score)
      } else { ## to save space with more than 10,000 different values use sampling
        Fn <- ecdf(sample(rawscores$score[!is.na(rawscores$score)], size=10000, replace=TRUE))
      }
    }
    metadata(obj) <- list(seqname=chr,
                          provider="UCSC",
                          provider_version="10Feb2014", ## it'd better to grab the date from downloaded file
                          download_url=downloadURL,
                          download_date=format(Sys.Date(), "%b %d, %Y"),
                          reference_genome=refgenomeGD,
                          data_pkgname="phyloP100way.UCSC.hg19",
                          qfun=.quantizer,
                          dqfun=.dequantizer,
                          ecdf=Fn,
                          max_abs_error=max.abs.error)
    saveRDS(obj, file=file.path("hg19.100way.phyloP", sprintf("phyloP100way.UCSC.hg19.%s.rds", chr)))
    rm(rawscores, obj)
    rm(list=sprintf("phyloP100way_%s", chr))
    gc()
  }, error=function(err) {
      message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
