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

downloadURL <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons7way"
citationdata <- bibentry(bibtype="Article",
                         author=c(person("Adam Siepel"), person("Gill Berejano"), person("Jakob S. Pedersen"),
                                  person("Angie S. Hinrichs"), person("Minmei Hou"), person("Kate Rosenbloom"),
                                  person("Hiram Clawson"), person("John Spieth"), person("LaDeana W. Hillier"),
                                  person("Stephen Richards"), person("George M. Weinstock"),
                                  person("Richard K. Wilson"), person("Richard A. Gibbs"),
                                  person("W. James Kent"), person("Webb Miller"), person("David Haussler")),
                         title="Evolutionarily conserved elements in vertebrate, insect, worm, and yeast genomes",
                         journal="Genome Research",
                         volume="15",
                         pages="1034-1050",
                         year="2005",
                         doi="10.1101/gr.3715005")

registerDoParallel(cores=4) ## each process may need up to 20Gb of RAM

## transform WIG to BIGWIG format
si <- Seqinfo(seqnames=seqnames(Hsapiens), seqlengths=seqlengths(Hsapiens))
wigToBigWig(file.path("hg38.7way.phastCons", "hg38.phastCons7way.wigFix.gz"), seqinfo=si)

pkgname <- "phastCons7way.UCSC.hg38"
dir.create(pkgname)

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

saveRDS(refgenomeGD, file=file.path(pkgname, "refgenomeGD.rds"))

## transform BIGWIG into Rle objects coercing phastCons scores into
## 1-decimal digit raw-encoded values to reduce memory requirements
## in principle deciles of phastCons probabilities should give the
## necessary resolution for the purpose of filtering genetic variants
## on conservation. We also include a version of the scores with
## 2-decimal digits.

## quantizer function for 1-decimal place. it maps input real-valued [0, 1]
## phastCons scores to non-negative integers [0, 255] so that
## each of them can be later coerced into a single byte (raw type).
## quantization is done by rounding to one decimal place,
## and therefore, mapping is restricted to 11 different positive
## integers only [1-11], where the 0 value is kept to code for missingness
## because there is no NA value in the raw class (Sec. 3.3.4 NA handling, R Language Definition)
.quantizer <- function(x) {
  q <- as.integer(sprintf("%.0f", 10*x))
  q <- q + 1L
  q
}
attr(.quantizer, "description") <- "multiply by 10, round to nearest integer, add up one"

## quantizer function for 2-decimal places. it maps input real-valued [0, 1]
## phastCons scores to non-negative integers [0, 255] so that
## each of them can be later coerced into a single byte (raw type).
## quantization is done by rounding to two decimal places,
## and therefore, mapping is restricted to 101 different positive
## integers only [1-101], where the 0 value is kept to code for missingness
## because there is no NA value in the raw class (Sec. 3.3.4 NA handling, R Language Definition)
.quantizerDP2 <- function(x) {
  q <- as.integer(sprintf("%.0f", 100*x))
  q <- q + 1L
  q
}
attr(.quantizerDP2, "description") <- "multiply by 100, round to nearest integer, add up one"

## dequantizer function for 1-decimal place. it maps input non-negative integers [0, 255]
## to real-valued phastCons scores, where 0 codes for NA
.dequantizer <- function(q) {
  x <- as.numeric(q)
  x[x == 0] <- NA
  x <- (x - 1) / 10
  x
}
attr(.dequantizer, "description") <- "subtract one integer unit, divide by 10"

## dequantizer function for 2-decimal places. it maps input non-negative integers [0, 255]
## to real-valued phastCons scores, where 0 codes for NA
.dequantizerDP2 <- function(q) {
  x <- as.numeric(q)
  x[x == 0] <- NA
  x <- (x - 1) / 100
  x
}
attr(.dequantizerDP2, "description") <- "subtract one integer unit, divide by 100"

nsites <- foreach (chr=seqnames(Hsapiens), .combine='c') %dopar% {
  cat(chr, "\n")
  tryCatch({
    rawscores <- import.bw(BigWigFile(file.path("hg38.7way.phastCons", "hg38.phastCons7way.bw")),
                           which=GRanges(seqnames=chr, IRanges(1, seqlengths(Hsapiens)[chr])))
    qscoresDP1 <- .quantizer(rawscores$score)
    qscoresDP2 <- .quantizerDP2(rawscores$score)
    max.abs.error.DP1 <- max(abs(rawscores$score - .dequantizer(qscoresDP1)))
    max.abs.error.DP2 <- max(abs(rawscores$score - .dequantizerDP2(qscoresDP2)))
    objDP1 <- coverage(rawscores, weight=qscoresDP1)[[chr]]
    objDP2 <- coverage(rawscores, weight=qscoresDP2)[[chr]]
    if (any(runValue(objDP1) > 0) && any(runValue(objDP2) > 0)) {
      runValue(objDP1) <- as.raw(runValue(objDP1))
      runValue(objDP2) <- as.raw(runValue(objDP2))
      Fn <- function(x) { warning("no ecdf() function available") ; numeric(0) }
      n <- length(unique(rawscores$score[!is.na(rawscores$score)]))
      if (n > 10) {
        if (n <= 10000) {
          Fn <- ecdf(rawscores$score)
        } else { ## to save space with more than 10,000 different values use sampling
          Fn <- ecdf(sample(rawscores$score[!is.na(rawscores$score)], size=10000, replace=TRUE))
        }
      }
      metadata(objDP1) <- list(seqname=chr,
                               provider="UCSC",
                               provider_version="04Jun2014", ## it'd better to grab the date from downloaded file
                               citation=citationdata,
                               download_url=downloadURL,
                               download_date=format(Sys.Date(), "%b %d, %Y"),
                               reference_genome=refgenomeGD,
                               data_pkgname=pkgname,
                               qfun=.quantizer,
                               dqfun=.dequantizer,
                               ecdf=Fn,
                               max_abs_error=max.abs.error.DP1)
      saveRDS(objDP1, file=file.path(pkgname, sprintf("phastCons7way.UCSC.hg38.%s.rds", chr)))
      metadata(objDP2) <- list(seqname=chr,
                               provider="UCSC",
                               provider_version="09Feb2014",
                               citation=citationdata,
                               download_url=downloadURL,
                               download_date=format(Sys.Date(), "%b %d, %Y"),
                               reference_genome=refgenomeGD,
                               data_pkgname=pkgname,
                               qfun=.quantizerDP2,
                               dqfun=.dequantizerDP2,
                               ecdf=Fn,
                               max_abs_error=max.abs.error.DP2)
      saveRDS(objDP2, file=file.path(pkgname, sprintf("phastCons7way.UCSC.hg38.DP2.%s.rds", chr)))
    }
    nsites <- as.numeric(sum(objDP1 > 0)) ## use 'numeric' to avoid integer overflow
    rm(rawscores, objDP1, objDP2)
    gc()
    nsites
  }, error=function(err) {
      message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}

saveRDS(sum(nsites), file=file.path(pkgname, "nsites.rds"))
