## This file explains the steps to download and freeze the AlphaMissense scores
## released in 2023 for the human genome version hg19. If you use these data on
## your own research please cite the following publication:

## Cheng J, Novati G, Pan J, Bycroft C, Žemgulytė A, Applebaum T, Pritzel A,
## Wong LH, Zielinski M, Sargeant T, Schneider RG, Senior, AW, Jumper J, Hassabis D,
## Kohli P, Avesec Ž. Accurate proteome-wide missense variant effect prediction
## with AlphaMissense. Science, 381:eadg7492, 2023. DOI: 10.1126/science.adg7492

## The data were downloaded from the following URL:
##
## $ wget -c https://zenodo.org/record/8208688/files/AlphaMissense_hg19.tsv.gz
##
## Build tabix index
##
## $ tabix -p vcf AlphaMissense_hg19.tsv.gz
##
## The data were first splitted into tabix files per chromosome as follows:
##
## mkdir -p AM_by_chr
## allchr=`tabix -l AlphaMissense_hg19.tsv.gz`
## for chr in $allchr ; do {
##   echo $chr
##   tabix -h AlphaMissense_hg19.tsv.gz $chr | bgzip -c > AM_by_chr/$chr.tsv.gz
##   if [ -s AM_by_chr/$chr.tsv.gz ] ; then
##     tabix -p vcf AM_by_chr/$chr.tsv.gz
##   else
##     rm AM_by_chr/$chr.tsv.gz
##   fi
## } done

## The following R script processes the downloaded and splitted data to
## store the AlphaMissense scores in raw-Rle objects

library(Rsamtools)
library(GenomeInfoDb)
library(GenomicRanges)
library(GenomicScores)
library(BSgenome.Hsapiens.UCSC.hg19)
library(doParallel)

downloadURL <- "https://zenodo.org/record/8208688/files/AlphaMissense_hg19.tsv.gz"
datacitation <- bibentry(bibtype="Article",
                         author=c(person("Jun Cheng"), person("Guido Novati"),
                                  person("Joshua Pan"), person("Clare Bycroft"),
                                  person("Akvilė Žemgulytė"),
                                  person("Taylor Applebaum"),
                                  person("Alexander Pritzel"),
                                  person("Lai Hong Wong"),
                                  person("Michal Zielinski"),
                                  person("Tobias Sargeant"),
                                  person("Rosalia G. Schneider"),
                                  person("Andrew W. Senior"),
                                  person("John Jumper"),
                                  person("Demis Hassabis"),
                                  person("Pushmeet Kohli"),
                                  person("Žiga Avsec")),
                         title="Accurate proteome-wide missense variant effect prediction with AlphaMissense",
                         journal="Science",
                         volume="381",
                         pages="eadg7492",
                         year="2023",
                         doi="10.1126/science.adg7492")

registerDoParallel(cores=10)

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
## x: values to quantize, x >= 0 & x <= 1, length(x) is multiple of d
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
  q[!is.na(x)] <- as.integer(sprintf("%.0f", 100*x[!is.na(x)]))
  if (max(q, na.rm=TRUE) > (n - 1))
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
  if (any(q < 0) || any(q > .Machine$integer.max))
    stop("current number of quantized values cannot be stored into one integer")
  q
}
attr(.quantizer, "description") <- "round to 2 digits, transform into integer, store multiple values per position using a given base"

## dequantizer function
.dequantizer <- function(q, ...) {
  d <- 1L ; b <- 10L ; na.zero <- FALSE
  otherArgs <- list(...)
  for (i in seq_along(otherArgs))
    assign(names(otherArgs)[i], otherArgs[[i]])
  x <- as.numeric(q)
  x[x == 0] <- NA ## 0 values correspond to genomic
  x <- x - 1      ## positions without scores
  if (d > 1L)
    x <- GenomicScores:::.fromBase(x, d=d, b=b)
  if (na.zero) { ## should we decode 0s as NAs
    x[x == 0] <- NA
    x <- x - 1L
  }
  x <- x / 100
  x
}
attr(.dequantizer, "description") <- "transform into base 10, multiply by 100"

## allele normalization of genomic score data
## input: data.frame with POS, REF, ALT, SCO for
## a given chromosome
## output: data.frame with POS, REF, ALT, SCO
## for every possible ALT, inserting NAs in the
## SCORE column for missing ALT and ALT are
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

tbx <- open(TabixFile("AlphaMissense_hg19.tsv.gz"))
tbxchr <- sortSeqlevels(seqnamesTabix(tbx))
close(tbx)

tbxchr <- tbxchr[tbxchr %in% seqlevels(Hsapiens)]
stopifnot(length(tbxchr) > 0) ## QC
allchrgr <- GRanges(seqnames=seqnames(Hsapiens)[match(tbxchr, seqnames(Hsapiens))],
                    IRanges(1, seqlengths(Hsapiens)[tbxchr]),
                    seqinfo=seqinfo(Hsapiens)[tbxchr])
names(allchrgr) <- seqnames(allchrgr)
seqlevelsStyle(allchrgr) <- seqlevelsStyle(tbxchr)[1]

foreach (chr=seqlevels(allchrgr)) %dopar% {
  tryCatch({
    cat(sprintf("processing chromosome %s\n", chr))
    ## read the whole TSV file for this chromosome into main memory
    tbx <- open(TabixFile(file.path("AM_by_chr", sprintf("%s.tsv.gz", chr))))
    obj <- Rle(integer(seqlengths(allchrgr)[chr]))
    rawscores <- scanTabix(tbx)[[1]]
    close(tbx)
    rawscores <- do.call("rbind", strsplit(rawscores, split="\t"))
    posrefalt <- sprintf("%s_%s_%s", rawscores[, 2], rawscores[, 3], rawscores[, 4])
    mask <- !duplicated(posrefalt) ## some rows duplicated due to transcript annotation
    rawscores <- data.frame(POS=as.integer(rawscores[mask, 2]),
                            REF=rawscores[mask, 3],
                            ALT=rawscores[mask, 4],
                            SCO=as.numeric(rawscores[mask, 9]),
                            stringsAsFactors=FALSE)
    nsites <- nrow(rawscores)

    rawscores <- anorm_score_data(rawscores)

    q <- .quantizer(rawscores$SCO, n=101L, d=3L, na.zero=TRUE)
    x <- .dequantizer(q=q, d=3L, b=102L, na.zero=TRUE)
    max.abs.error <- max(abs(x - rawscores$SCO), na.rm=TRUE)

    uniqpos <- unique(rawscores$POS)
    gr <- sprintf("%s:%s-%s", chr, uniqpos, uniqpos)
    gr <- GRanges(gr)
    seqinfo(gr) <- seqinfo(allchrgr)[chr]
    obj <- coverage(gr, weight=q)[[1]]
    rawscores <- rawscores$SCO[!is.na(rawscores$SCO)]
    if (length(unique(rawscores)) <= 10000) {
      Fn <- ecdf(rawscores)
    } else { ## to save space with more than 10,000 different values use sampling
      Fn <- ecdf(sample(rawscores, size=10000, replace=TRUE))
    }
    rm(rawscores)
    gc()
    chr <- names(allchrgr)[which(seqlevels(allchrgr) == chr)]
    metadata(obj) <- list(seqname=chr,
                          provider="Google DeepMind",
                          provider_version="v2023",
                          citation=datacitation,
                          download_url=downloadURL,
                          download_date=format(Sys.Date(), "%b %d, %Y"),
                          reference_genome=refgenomeGD,
                          data_pkgname="AlphaMissense.v2023.hg19",
                          qfun=.quantizer,
                          qfun_args=list(n=101L, d=3L, na.zero=TRUE),
                          dqfun=.dequantizer,
                          dqfun_args=list(d=3L, b=102L, na.zero=TRUE),
                          valxpos=3L,
                          ecdf=Fn,
                          max_abs_error=max.abs.error)
    saveRDS(obj, file=sprintf("AlphaMissense.%s.hg19.%s.rds",
                              metadata(obj)$provider_version, chr))
    saveRDS(nsites, file=sprintf("AlphaMissense.%s.hg19.%s.nsites.rds",
                                 metadata(obj)$provider_version, chr))
    rm(obj)
    gc()
  }, error=function(err) {
    message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}

nsitesfiles <- list.files(path=".", pattern="nsites.rds$", full.names=TRUE)
nsites <- 0
for (fname in nsitesfiles) {
  nsites <- nsites + readRDS(fname)
  unlink(fname)
}
saveRDS(nsites, file="nsites.rds")
