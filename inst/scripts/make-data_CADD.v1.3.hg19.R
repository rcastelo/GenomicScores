## This file explains the steps to download and freeze the CADD scores
## for the human genome version hg19. If you use these data on your own
## research please cite the following publication:

## Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J.
## A general framework for estimating the relative pathogenicity of human
## genetic variants. Nature Genetics, 46:310-315, 2014. DOI: 10.1038/ng.2892

## The data were downloaded from the following URL:
##
## $ wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz
## $ wget -c http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs.tsv.gz.tbi
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

library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(doParallel)

downloadURL <- "http://krishna.gs.washington.edu/download/CADD/v1.3"

registerDoParallel(cores=8)

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
## x: values to quantize, x >= 0 & x <= 99, length(x) is multiple of d
## n: maximum number of quantized values
## d: number of digits in 'x' forming a value to quantize
.quantizer <- function(x, n=Inf, d=1L, na.zero=FALSE) {
  stopifnot(d > 0L) ## QC
  stopifnot(length(x) %% d == 0) ## QC
  q <- as.integer(sprintf("%.0f", x / 10L))
  if (max(q, na.rm=TRUE) > (n - 1))
    q[q > (n - 1)] <- n - 1
  base <- n
  if (na.zero) { ## should we recode NAs into 0s?
    q <- q + 1L
    q[is.na(q)] <- 0L
    base <- base + 1L
  }
  if (d > 1)
    q <- .toBase(d=q, g=rep(1:(length(x)/d), each=d), b=base)
  q <- q + 1L
  if (any(q > 255))
    stop("current number of quantized values > 255 and cannot be stored into one byte")
  q
}
attr(.quantizer, "description") <- "quantize PHRED scores into deciles and store them using a given base"

## dequantizer function
.dequantizer <- function(q, d, b, na.zero=FALSE) {
  q <- q - 1L
  q <- .fromBase(q, d=d, b=b)
  if (na.zero)
    q[q == 0L] <- NA
  q <- 10L * q
  q
}
attr(.dequantizer, "description") <- "transform into base 10 digits, multiply by 10"

tbx <- open(TabixFile("whole_genome_SNVs.tsv.gz"))
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
    tbx <- open(TabixFile(file.path("CADD_by_chr", sprintf("chr%s.tsv.gz", chr)),
                          yieldSize=999999)) ## yieldSize must be multiple of 3!
    obj <- Rle(integer(seqlengths(allchrgr)[chr]))
    max.abs.error <- 0
    allrawscores <- c()
    i <- 0
    while (length(rawscores <- scanTabix(tbx)[[1]])) {
      if (length(rawscores) %% 3 != 0)
        stop(sprintf("CADD scores in chromosome %s are not multiple of 3", chr))
      i <- i + length(rawscores) / 3
      cat(sprintf("processing chr%s %d positions of about %d (~%.1f%%)\n",
                  chr, i, seqlengths(allchrgr)[chr], 100*i/seqlengths(allchrgr)[chr]))
      rawscores <- do.call("rbind", strsplit(rawscores, split="\t"))
      rawscores <- data.frame(POSITION=as.integer(rawscores[, 2]),
                              REF=rawscores[, 3],
                              ALT=rawscores[, 4],
                              SCORE=as.numeric(rawscores[ ,6]),
                              stringsAsFactors=FALSE)
      ## remove positions with missing reference allele (e.g., 3:60830534, 3:60830763, 3:60830764)
      rawscores <- rawscores[rawscores$REF %in% c("A", "C", "G", "T"), ]
      nbyp <- split(rawscores$ALT, rawscores$POSITION)
      ## remove positions with unsorted or missing alleles
      unsortedpos <- as.integer(names(nbyp)[!sapply(nbyp, isSorted)])
      nonmult3pos <- as.integer(names(nbyp)[elementNROWS(nbyp) != 3])
      rawscores <- rawscores[!rawscores$POSITION %in% unique(c(unsortedpos, nonmult3pos)), ]
      rm(nbyp, unsortedpos, nonmult3pos)
      gc()

      ## quantize raw PHRED scores into 6 different values, setting
      ## PHRED scores > 50 to 50, and convert tuples of 3 digits into
      ## a base-6 integer number < 255
      q <- .quantizer(rawscores$SCORE, n=6L, d=3L, na.zero=FALSE)
      x <- .dequantizer(q=q, d=3L, b=6L, na.zero=FALSE)
      ## calculate maximum absolute error for PHRED scores <= 50
      ## b/c those > 50 have been set to 50
      mask50 <- rawscores$SCORE <= 50
      max.abs.error <- max(c(max.abs.error, abs(rawscores$SCORE[mask50] - x[mask50])), na.rm=TRUE)
      if (max.abs.error > 5) cat(sprintf("chr=%s max.abs.error=%f, i=%d\n", chr, max.abs.error, i))
      allrawscores <- c(allrawscores, rawscores$SCORE)
      uniqpos <- unique(rawscores$POSITION)
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
                          provider_version="v1.3",
                          download_url=downloadURL,
                          download_date=format(Sys.Date(), "%b %d, %Y"),
                          reference_genome=refgenomeGD,
                          data_pkgname="cadd.v1.0.hg19",
                          qfun=.quantizer,
                          dqfun=.dequantizer,
                          ecdf=Fn,
                          max_abs_error=max.abs.error)
    saveRDS(obj, file=sprintf("cadd.hg19.%s.rds", chr))
    rm(allrawscores, obj)
    rm(obj)
    gc()
  }, error=function(err) {
    message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
