## This file explains the steps to download and freeze the phastCons
## conservation scores for human genome version hg38, calculated on
## 7 vertebrate species. If you use these data on your own research
## please cite the following publication:

## Siepel A, Bejerano G, Pedersen JS, Hinrichs AS, Hou M, Rosenbloom K, Clawson H,
## Spieth J, Hillier LW, Richards S, et al. Evolutionarily conserved elements
## in vertebrate, insect, worm, and yeast genomes. Genome Res. 2005 Aug;15(8):1034-50.
## (http://www.genome.org/cgi/doi/10.1101/gr.3715005)

## The data were downloaded from the UCSC genome browser with the Unix 'rsync'
## command as follows
##
## $ rsync -avz --progress \
##     rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons7way/hg38.phastCons7way.wigFix.gz \
##     ./hg38.7way.phastCons

## The following R script processes the downloaded data to
## store the phastCons scores in raw-Rle objects

library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(doParallel)
library(S4Vectors)

downloadURL <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons7way"

registerDoParallel(cores=4) ## each process may need up to 20Gb of RAM

## transform WIG to BIGWIG format
si <- Seqinfo(seqnames=seqnames(Hsapiens), seqlengths=seqlengths(Hsapiens))
wigToBigWig(file.path("hg38.7way.phastCons", "hg38.phastCons7way.wigFix.gz"), seqinfo=si)

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

saveRDS(refgenomeGD, file=file.path("hg38.7way.phastCons", "refgenomeGD.rds"))

## transform BIGWIG into Rle objects coercing phastCons scores into
## 1-decimal digit raw-encoded values to reduce memory requirements
## in principle deciles of phastCons probabilities should give the
## necessary resolution for the purpose of filtering genetic variants
## on conservation

## quantizer function. it maps input real-valued [0, 1]
## phastCons scores to non-negative integers [0, 255] so that
## each of them can be later coerced into a single byte (raw type).
## quantization is done by rounding to one decimal significant digit,
## and therefore, mapping is restricted to 11 different positive
## integers only [1-11], where the 0 value is kept to code for missingness
.quantizer <- function(x) {
  q <- as.integer(sprintf("%.0f", 10*x))
  q <- q + 1L
  q
}
attr(.quantizer, "description") <- "multiply by 10, round to nearest integer, add up one"

## dequantizer function. it maps input non-negative integers [0, 255]
## to real-valued phastCons scores, where 0 codes for NA
.dequantizer <- function(q) {
  x <- as.numeric(q)
  x[x == 0] <- NA
  x <- (x - 1) / 10
  x
}
attr(.dequantizer, "description") <- "subtract one integer unit, divide by 10"

foreach (chr=seqnames(Hsapiens)) %dopar% {
  cat(chr, "\n")
  tryCatch({
    rawscores <- import.bw(BigWigFile(file.path("hg38.7way.phastCons", "hg38.phastCons7way.bw")),
                           which=GRanges(seqnames=chr, IRanges(1, seqlengths(Hsapiens)[chr])))
    qscores <- .quantizer(rawscores$score)
    max.abs.error <- max(abs(rawscores$score - .dequantizer(qscores)))
    assign(sprintf("phastCons7way_%s", chr), coverage(rawscores, weight=qscores)[[chr]])
    assign(sprintf("phastCons7way_%s", chr),
           do.call("runValue<-", list(get(sprintf("phastCons7way_%s", chr)),
                                      as.raw(runValue(get(sprintf("phastCons7way_%s", chr)))))))
    obj <- get(sprintf("phastCons7way_%s", chr))
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
                          provider_version="04Jun2014", ## it'd better to grab the date from downloaded file
                          download_url=downloadURL,
                          download_date=format(Sys.Date(), "%b %d, %Y"),
                          reference_genome=refgenomeGD,
                          data_pkgname="phastCons7way.UCSC.hg38",
                          qfun=.quantizer,
                          dqfun=.dequantizer,
                          ecdf=Fn,
                          max_abs_error=max.abs.error)
    saveRDS(obj, file=file.path("hg38.7way.phastCons", sprintf("phastCons7way.UCSC.hg38.%s.rds", chr)))
    rm(rawscores, obj)
    rm(list=sprintf("phastCons7way_%s", chr))
    gc()
  }, error=function(err) {
    message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
