## This file explains the steps to download and freeze the phastCons
## conservation scores for human genome version hg19, calculated on
## 100 vertebrate species. If you use these data on your own research
## please cite the following publication:

## Siepel A, Bejerano G, Pedersen JS, Hinrichs AS, Hou M, Rosenbloom K, Clawson H,
## Spieth J, Hillier LW, Richards S, et al. Evolutionarily conserved elements
## in vertebrate, insect, worm, and yeast genomes. Genome Res. 2005 Aug;15(8):1034-50.
## (http://www.genome.org/cgi/doi/10.1101/gr.3715005)

## The data were downloaded from the UCSC genome browser with the Unix 'rsync'
## command as follows
##
## $ rsync -avz --progress \
##     rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/ \
##     ./hg19.100way.phastCons

## The following R script processes the downloaded data to
## store the phastCons scores in raw-Rle objects

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(doParallel)
library(S4Vectors)

downloadURL <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons"

registerDoParallel(cores=4) ## each process may need up to 20Gb of RAM

## transform WIG to BIGWIG format
si <- Seqinfo(seqnames=seqnames(Hsapiens), seqlengths=seqlengths(Hsapiens))
foreach (chr=seqnames(Hsapiens)) %dopar% {
  cat(chr, "\n")
  wigToBigWig(file.path("hg19.100way.phastCons", sprintf("%s.phastCons100way.wigFix.gz", chr)), seqinfo=si)
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

saveRDS(refgenomeGD, file=file.path("hg19.100way.phastCons", "refgenomeGD.rds"))

## transform BIGWIG into Rle objects coercing phastCons scores into
## 1-decimal digit raw-encoded values to reduce memory requirements
## in principle deciles of phastCons probabilities should give the
## necessary resolution for the purpose of filtering genetic variants
## on conservation

## quantization function. it maps input real-valued [0, 1]
## phastCons scores to non-negative integers [0, 255] so that
## each of them can be later coerced into a single byte (raw type).
## quantization is done by rounding to one decimal significant digit,
## and therefore, mapping is restricted to 11 different non-negative
## integers only [0-10]
quantization <- function(x) {
  q <- 10*round(x, digits=1)
  q
}
attr(quantization, "description") <- "10*round(x, digits=1)"

## inverse quantization function. it maps input non-negative integers [0, 255]
## to real-valued phastCons scores
inversequantization <- function(q) {
  q/10
}
attr(inversequantization, "description") <- "q/10"

foreach (chr=seqnames(Hsapiens)) %dopar% {
  cat(chr, "\n")
  tryCatch({
    rawscores <- import.bw(BigWigFile(file.path("hg19.100way.phastCons", sprintf("%s.phastCons100way.bw", chr))),
                           as="RleList")[[chr]]
    qscores <- quantization(rawscores)
    assign(sprintf("phastCons100way_%s", chr), qscores)
    assign(sprintf("phastCons100way_%s", chr),
           do.call("runValue<-", list(get(sprintf("phastCons100way_%s", chr)),
                                      as.raw(runValue(get(sprintf("phastCons100way_%s", chr)))))))
    obj <- get(sprintf("phastCons100way_%s", chr))
    if (length(unique(rawscores)) <= 10000)
      Fn <- ecdf(decode(rawscores))
    else ## to save space with more than 10,000 different values use sampling
      Fn <- ecdf(sample(decode(rawscores), size=10000, replace=TRUE))
    metadata(obj) <- list(seqname=chr,
                          provider="UCSC",
                          provider_version="09Feb2014", ## it'd better to grab the date from downloaded file
                          download_url=downloadURL,
                          download_date=format(Sys.Date(), "%b %d, %Y"),
                          reference_genome=refgenomeGD,
                          data_pkgname="phastCons100way.UCSC.hg19",
                          qfun=quantization,
                          iqfun=inversequantization,
                          ecdf=Fn,
                          max_abs_error=max(abs(rawscores-inversequantization(qscores))))
    saveRDS(obj, file=file.path("hg19.100way.phastCons", sprintf("phastCons100way.UCSC.hg19.%s.rds", chr)))
    rm(rawscores)
    rm(obj)
    rm(list=sprintf("phastCons100way_%s", chr))
    gc()
  }, error=function(err) {
      message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
