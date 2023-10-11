#!/usr/bin/env /projects_fg/soft/R/release/bin/Rscript

## This file explains the steps to download and freeze the CADD scores v1.6
## for the human genome version hg38. If you use these data on your own
## research please cite the following publication:

## Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J.
## A general framework for estimating the relative pathogenicity of human
## genetic variants. Nature Genetics, 46:310-315, 2014. DOI: 10.1038/ng.2892

## The data were downloaded from the following URL:
##
## $ wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
## $ wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
##
## The data were first splitted into tabix files per chromosome as follows:
##
## mkdir -p CADD_by_chr
## allchr=`tabix -l whole_genome_SNVs.tsv.gz`
## for chr in $allchr ; do {
##   echo chr$chr
##   tabix -h whole_genome_SNVs.tsv.gz $chr | bgzip -c > CADD_by_chr/chr$chr.tsv.gz
##   if [ -s CADD_by_chr/chr$chr.tsv.gz ] ; then
##     tabix -p vcf CADD_by_chr/chr$chr.tsv.gz
##   else
##     rm CADD_by_chr/chr$chr.tsv.gz
##   fi
## } done

## The following R script processes the downloaded and splitted data to
## store the CADD scores in raw-Rle objects

stopifnot(BiocManager::version() == "3.17")

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    stop("make-data_cadd.v1.6.GRCh38.R <PATH2TBXFILES>")
}
path2tbxfiles <- args[1]
if (!file.exists(path2tbxfiles))
          stop(sprintf("%s does not exist.", path2tbxfiles))

library(Rsamtools)
library(GenomeInfoDb)
library(GenomicRanges)
library(GenomicScores)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)

downloadURL <- "https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38"
caddVersion <- "v1.6"
datacitation <- bibentry(bibtype="Article",
                         author=c(person("Martin Kircher"), person("Daniela M. Witten"), person("Preti Jain"),
                                  person("Brian J. O'Roak"), person("Gregory M. Cooper"), person("Jay Shendure")),
                         title="A general framework for estimating the relative pathogenicity of human genetic variants",
                         journal="Nature Genetics",
                         volume="46",
                         pages="310-315",
                         year="2014",
                         doi="10.1038/ng.2892")

registerDoParallel(cores=25)

## freeze the GenomeDescription data for Hsapiens

refgenomeGD <- GenomeDescription(organism=organism(Hsapiens),
                                 common_name=commonName(Hsapiens),
                                 provider=provider(Hsapiens),
                                 provider_version=metadata(Hsapiens)$genome,
                                 release_date=releaseDate(Hsapiens),
                                 release_name=metadata(Hsapiens)$genome,
                                 seqinfo=keepStandardChromosomes(seqinfo(Hsapiens)))

saveRDS(refgenomeGD, file="refgenomeGD.rds")

## quantizer function
## x: values to quantize, x >= 0 & x <= 99, length(x) is multiple of d
## n: maximum number of quantized values
## d: number of digits in 'x' forming a value to quantize
.quantizer <- function(x, ...) {
  n <- Inf ; d <- 1L ; na.zero <- FALSE
  otherArgs <- list(...)
  for (i in seq_along(otherArgs))
    assign(names(otherArgs)[i], otherArgs[[i]])
  stopifnot(d > 0L) ## QC
  stopifnot(length(x) %% d == 0) ## QC
  q <- rep(NA_integer_, length(x))
  q[!is.na(x)] <- as.integer(sprintf("%.0f", x[!is.na(x)]/10L))
  if (max(q, na.rm=TRUE) > (n - 1)) ## cap at n-1
    q[q > (n - 1)] <- n - 1
  base <- n
  if (na.zero) { ## should we recode NAs into 0s?
    q <- q + 1L
    q[is.na(q)] <- 0L
    base <- base + 1L
  }
  if (d > 1)
    q <- GenomicScores:::.toBase(d=q, g=rep(1:(length(x)/d), each=d), b=base)
  q <- q + 1L ## 0 values correspond to genomic
              ## positions without scores
  q[is.na(q)] <- 0L ## any NA left is recoded into a genomic position without any score
  if (any(q > 255))
    stop("current number of quantized values > 255 and cannot be stored into one byte")
  q
}
attr(.quantizer, "description") <- "quantize PHRED scores into deciles and store them using a given base"

## dequantizer function
.dequantizer <- function(q, ...) {
  d <- 1L ; b <- 10L ; na.zero=FALSE
  otherArgs <- list(...)
  for (i in seq_along(otherArgs))
    assign(names(otherArgs)[i], otherArgs[[i]])
  x <- as.numeric(q)
  x[x == 0] <- NA ## 0 values correspond to genomic
  x <- x - 1      ## positions without scores
  if (d > 1L)
    x <- GenomicScores:::.fromBase(x, d=d, b=b)
  if (na.zero) { ## should we decode 0s as NAs
    x[x == 0L] <- NA
    x <- x - 1L
  }
  x <- 10L * x
  x
}
attr(.dequantizer, "description") <- "transform into base-10 digits, multiply by 10"

## allele normalization of genomic score data
## input: data.frame with POS, REF, ALT, SCO for
## a given chromosome
## output: data.frame with POS, REF, ALT, SCO
## for every possible ALT, inserting NAs in the
## SCO column for missing ALT and ALT are
## alphabetically ordered within each position
anorm_score_data <- function(dat) {
  stopifnot(all(dat$REF %in% c("A", "C", "G", "T"))) ## QC
  abyp <- split(dat$ALT, dat$POS) ## ALT-by-POS
  abyp <- abyp[order(as.integer(names(abyp)))] ## just in case split orders arbitrarily
  sbyp <- relist(dat$SCO, abyp) ## SCO-by-POS
  rbyp <- relist(dat$REF, abyp) ## REF-by-POS

  ## fill-in NAs on missing ALT
  abyp <- lapply(abyp, "length<-", 3)
  sbyp <- lapply(sbyp, "length<-", 3)
  rbyp <- lapply(rbyp, "length<-", 3)

  ## fill-in REF nucleotide on missing ALT
  urbyp <- unlist(rbyp, use.names=FALSE)
  wh <- which(is.na(urbyp))
  urbyp[wh] <- urbyp[wh-1]
  wh <- which(is.na(urbyp))
  urbyp[wh] <- urbyp[wh-1]
  rbyp <- relist(urbyp, rbyp)
  rm(urbyp)

  ## order ALT nucleotides and their scores,
  ## and put the result into a data.frame
  i <- 0
  ndat <- mapply(function(r, a, s, nt) {
                   i <<- i + 1
                   repnt <- setdiff(nt, c(a, r[1]))
                   stopifnot(sum(is.na(a)) == length(repnt)) ## QC
                   a[is.na(a)] <- repnt
                   p <- order(a)
                   data.frame(REF=r, ALT=a[p], SCO=s[p])
                 }, rbyp, abyp, sbyp,
                 MoreArgs=list(nt=c("A", "C", "G", "T")),
                 SIMPLIFY=FALSE)
  pos <- as.integer(names(ndat))
  names(ndat) <- NULL
  ndat <- do.call("rbind", ndat)
  ndat$POS <- rep(pos, each=3)
  ndat <- ndat[, c("POS", "REF", "ALT", "SCO")]
  ndat
}

tbx <- open(TabixFile(file.path(path2tbxfiles, "whole_genome_SNVs.tsv.gz")))
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
  tryCatch({
    tbx <- open(TabixFile(file.path(path2tbxfiles, "CADD_by_chr", sprintf("chr%s.tsv.gz", chr)),
                          yieldSize=1000000))
    obj <- Rle(integer(seqlengths(allchrgr)[chr]))
    max.abs.error <- 0
    allrawscores <- c()
    buffer <- data.frame(POS=integer(0), REF=character(0),
                         ALT=character(0), SCO=numeric(0))
    while (length(rawscores <- scanTabix(tbx)[[1]])) {
      rawscores <- do.call("rbind", strsplit(rawscores, split="\t"))
      rawscores <- data.frame(POS=as.integer(rawscores[, 2]),
                              REF=rawscores[, 3],
                              ALT=rawscores[, 4],
                              SCO=as.numeric(rawscores[ ,6]),
                              stringsAsFactors=FALSE)

      ## put together alleles from a common position splitted
      ## in two consecutive chunks
      chromposbuf <- paste(buffer$POS, buffer$REF, buffer$ALT, sep="_")
      chrompos <- paste(rawscores$POS[1:2], rawscores$REF[1:2], rawscores$ALT[1:2], sep="_")
      mt <- match(chromposbuf, chrompos)
      if (any(is.na(mt)))
          rawscores <- rbind(buffer[is.na(mt), , drop=FALSE], rawscores)
      rawscores <- anorm_score_data(rawscores)
      cat(sprintf("processing chr%s (~%.1f%%)\n", chr,
                  100*rawscores$POS[nrow(rawscores)]/seqlengths(allchrgr)[chr]))

      ## quantize raw PHRED scores into 6 different values, setting
      ## PHRED scores > 50 to 50, and convert tuples of 3 digits into
      ## a base-6 integer number < 255
      q <- .quantizer(rawscores$SCO, n=6L, d=3L, na.zero=FALSE)
      ## dequantize for sanity check, one extra base for NAs
      x <- .dequantizer(q=q, d=3L, b=6L, na.zero=FALSE)

      ## calculate maximum absolute error for PHRED scores <= 50
      ## b/c those > 50 have been set to 50
      mask50 <- rawscores$SCO <= 50
      max.abs.error <- max(c(max.abs.error, abs(rawscores$SCO[mask50] - x[mask50])), na.rm=TRUE)
      if (max.abs.error > 5) cat(sprintf("chr=%s max.abs.error=%.1f, pos=%d to %d\n", chr,
                                         max.abs.error, rawscores$POS[1], rawscores$POS[nrow(rawscores)]))
      allrawscores <- c(allrawscores, rawscores$SCO)
      uniqpos <- unique(rawscores$POS)

      mask <- is.na(rawscores$SCO)
      if (any(mask)) {
        buffer <- rawscores[mask, , drop=FALSE]
        rownames(buffer) <- NULL
        uniqpos <- uniqpos[-length(uniqpos)]
        q <- q[-length(q)]
      } else
        buffer <- data.frame(POS=integer(0), REF=character(0),
                             ALT=character(0), SCO=numeric(0))

      gr <- sprintf("%s:%s-%s", chr, uniqpos, uniqpos)
      gr <- GRanges(gr)
      seqinfo(gr) <- seqinfo(allchrgr)[chr]
      obj <- obj + coverage(gr, weight=q)[[1]]
    }
    close(tbx)
    runValue(obj) <- as.raw(runValue(obj))
    rm(rawscores)
    gc()
    if (length(unique(allrawscores)) <= 10000) {
      Fn <- ecdf(allrawscores)
    } else { ## to save space with more than 10,000 different values use sampling
      Fn <- ecdf(sample(allrawscores, size=10000, replace=TRUE))
    }
    chr <- names(allchrgr)[which(seqlevels(allchrgr) == chr)]
    metadata(obj) <- list(seqname=chr,
                          provider="UWashington",
                          provider_version=caddVersion,
                          citation=datacitation,
                          download_url=downloadURL,
                          download_date=format(Sys.Date(), "%b %d, %Y"),
                          reference_genome=refgenomeGD,
                          data_pkgname=sprintf("cadd.%s.%s", caddVersion, unname(genome(refgenomeGD)[1])),
                          qfun=.quantizer,
                          qfun_args=list(n=6L, d=3L, na.zero=FALSE),
                          dqfun=.dequantizer,
                          dqfun_args=list(d=3L, b=6L, na.zero=FALSE),
                          valxpos=3L,
                          ecdf=Fn,
                          max_abs_error=max.abs.error)
    saveRDS(obj, file=sprintf("cadd.%s.%s.%s.rds", metadata(obj)$provider_version,
                              unname(genome(refgenomeGD)[1]), chr))
    rm(allrawscores, obj)
    rm(obj)
    gc()
  }, error=function(err) {
    message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
