## This file explains the steps to download and freeze the LINSIGHT
## scores for human genome version hg19. If you use these data on
## your own research please cite the following publication:

## Huang Y-F, Gulko B, Siepel A. Fast, scalable prediction of deleterious
## noncoding variants from functional and population genomic data.
## Nat. Genet. 2015 Aug;47:276-83.
## (http://www.nature.com/articles/ng.3810)

## The data were downloaded from the CSHL mirror of the UCSC genome browser as follows
##
## wget http://compgen.cshl.edu/%7Eyihuang/tracks/LINSIGHT.bw
##
## further information about these data can be found at the following URL:
##
## http://compgen.cshl.edu/~yihuang/LINSIGHT
##
## and at the Adam Siepel Lab website
##
## http://siepellab.labsites.cshl.edu

## The following R script processes the downloaded data to
## store the fitCons scores in raw-Rle objects

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(doParallel)
library(S4Vectors)

downloadURL <- "http://compgen.cshl.edu/%7Eyihuang/tracks/LINSIGHT.bw"
citationdata <- bibentry(bibtype="Article",
                         author=c(person("Yi-Fei Huang"), person("Brad Gulko"),
                                  person("Adam Siepel")),
                         title="Fast, scalable prediction of deleterious noncoding variants from functional and population genomic data",
                         journal="Nature Genetics",
                         volume="49",
                         pages="618-624",
                         year="2017",
                         doi="10.1038/ng.3810")

registerDoParallel(cores=4)

pkgname <- "linsight.UCSC.hg19"
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

## transform BIGWIG into Rle objects coercing LINSIGHT scores into
## 2-decimal digit raw-encoded values to reduce memory requirements
## in principle centiles of LINSIGHT probabilities should give the
## necessary resolution for the purpose of filtering genetic variants
## with LINSIGHT scores

## quantizer function. it maps input real-valued [0, 1]
## LINSIGHT scores to non-negative integers [0, 255] so that
## each of them can be later coerced into a single byte (raw type).
## quantization is done by rounding to two decimal significant digits,
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

nsites <- foreach (chr=seqnames(Hsapiens), .combine='c') %dopar% {
  cat(chr, "\n")
  tryCatch({
    rawscores <- import.bw(BigWigFile("LINSIGHT.bw"),
                           which=GRanges(seqnames=chr, IRanges(1, seqlengths(Hsapiens)[chr])))
    qscores <- .quantizer(rawscores$score)
    max.abs.error <- max(abs(rawscores$score - .dequantizer(qscores)))
    obj <- coverage(rawscores, weight=qscores)[[chr]]
    if (any(runValue(obj) > 0)) {
      runValue(obj) <- as.raw(runValue(obj))
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
                            provider_version="19Aug2016",
                            citation=citationdata,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            data_group="Fitness",
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn,
                            max_abs_error=max.abs.error)
      saveRDS(obj, file=file.path(pkgname, sprintf("linsight.UCSC.hg19.%s.rds", chr)))
    }
    nsites <- sum(runValue(obj) > 0)
    rm(rawscores, obj)
    gc()
    nsites
  }, error=function(err) {
    message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}

saveRDS(sum(nsites), file=file.path(pkgname, "nsites.rds"))
