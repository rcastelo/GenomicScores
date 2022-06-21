## This file explains the steps to download and freeze the phyloP
## conservation scores for mouse genome version mm39, calculated on
## 35 vertebrate species. If you use these data on your own research
## please cite the following publication:

## Siepel A, Bejerano G, Pedersen JS, Hinrichs AS, Hou M, Rosenbloom K, Clawson H,
## Spieth J, Hillier LW, Richards S, et al. Evolutionarily conserved elements
## in vertebrate, insect, worm, and yeast genomes. Genome Res. 2005 Aug;15(8):1034-50.
## (http://www.genome.org/cgi/doi/10.1101/gr.3715005)

## The data were downloaded from the UCSC genome browser with the Unix 'rsync'
## command as follows
##
## $ rsync -avz --progress \
##     rsync://hgdownload.soe.ucsc.edu/goldenPath/mm39/phyloP35way/mm39.35way.phyloP/ \
##     ./mm39.35way.phyloP

## The following R script processes the downloaded data to
## store the phastCons scores in raw-Rle objects

library(plyr) ## for function round_any
library(BSgenome.Mmusculus.UCSC.mm39)
library(rtracklayer)
library(doParallel)

downloadURL <- "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/phyloP35way/mm39.35way.phyloP"
citationdata <- bibentry(bibtype="article",
                         author=c(person("Katherine S. Pollard"), person("Melissa J. Hubisz"),
                                  person("Kate R. Rosenbloom"), person("Adam Siepel")),
                         title="Detection of nonneutral substitution rates on mammalian phylogenies",
                         journal="Genome Research",
                         volume="20",
                         pages="110-121",
                         year="2010",
                         doi="10.1101/gr.097857.109")

registerDoParallel(cores=4) ## each process may need up to 20Gb of RAM

bsgenomeobj <- get("Mmusculus")

## transform WIG to BIGWIG format
si <- seqinfo(bsgenomeobj)
foreach (chr=seqnames(bsgenomeobj)) %dopar% {
  cat(chr, "\n")
  wigToBigWig(file.path("mm39.35way.phyloP", sprintf("%s.phyloP35way.wigFix.gz", chr)),
              seqinfo=si)
}

pkgname <- "phyloP35way.UCSC.mm39"
dir.create(pkgname)

## freeze the GenomeDescription data from the BSGenome object

refgenomeGD <- GenomeDescription(organism=organism(bsgenomeobj),
                                 common_name=commonName(bsgenomeobj),
                                 provider=provider(bsgenomeobj),
                                 provider_version=metadata(bsgenomeobj)$genome,
                                 release_date=releaseDate(bsgenomeobj),
                                 release_name=metadata(bsgenomeobj)$genome,
                                 seqinfo=si)

saveRDS(refgenomeGD, file=file.path(pkgname, "refgenomeGD.rds"))

## transform BIGWIG into Rle objects coercing phyloP scores into
## 1-decimal digit raw-encoded values to reduce memory requirements
## in principle deciles of phyloP probabilities should give the
## necessary resolution for the purpose of filtering genetic variants
## on conservation

## quantizer function. it maps input real-valued [-Inf, Inf]
## phyloP scores to non-negative integers [0, 255] so that
## each of them can be later coerced into a single byte (raw type).
## according to http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=652129995_dqrcEjjmKYIKOWyUPl4WNzQRnxjy&g=cons60way
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


nsites <- foreach (chr=seqnames(bsgenomeobj), .combine='c') %dopar% {
  cat(chr, "\n")
  tryCatch({
    rawscores <- import.bw(BigWigFile(file.path("mm39.35way.phyloP",
                                                sprintf("%s.phyloP35way.bw", chr))))
    qscores <- .quantizer(rawscores$score)
    ## calculate maximum absolute error for absolute phyloP scores <= 10
    ## b/c those larger than 10 have been set to 10
    f <- cut(abs(rawscores$score), breaks=c(0, 2, max(abs(rawscores$score))), include.lowest=TRUE)
    err <- abs(rawscores$score - .dequantizer(qscores))
    mask10 <- abs(rawscores$score) <= 10
    max.abs.error <- tapply(err[mask10], f[mask10], max, na.rm=TRUE)
    assign(sprintf("phyloP60way_%s", chr), coverage(rawscores, weight=qscores)[[chr]])
    assign(sprintf("phyloP60way_%s", chr),
           do.call("runValue<-", list(get(sprintf("phyloP60way_%s", chr)),
                                      as.raw(runValue(get(sprintf("phyloP60way_%s", chr)))))))
    obj <- get(sprintf("phyloP60way_%s", chr))
    if (any(runValue(obj) > 0)) {
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
                            provider_version="22Dec2020",
                            citation=citationdata,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn,
                            max_abs_error=max.abs.error)
      saveRDS(obj, file=file.path(pkgname, sprintf("%s.%s.rds", pkgname, chr)))
    }
    nsites <- as.numeric(sum(obj > 0)) ## use 'numeric' to avoid integer overflow
    rm(rawscores, obj)
    rm(list=sprintf("phyloP60way_%s", chr))
    gc()
    nsites
  }, error=function(err) {
      message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}

saveRDS(sum(nsites), file=file.path(pkgname, "nsites.rds"))
