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

downloadURL <- "http://compgen.cshl.edu/fitCons/0downloads/tracks/V1.01/i6/scores/fc-i6-0.bw"

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

## import BIGWIG file of fitCons scores as a GRanges object
fitConsGR <- import.bw(BigWigFile("fc-i6-0.bw"))

## get the regions without fitCons scores by complementing the genomic ranges
fitConsGRcompl <- gaps(fitConsGR)

## since fitCons scores have no strand we should eliminate ranges
## on the positive and negative strand spanning whole chromosomes
fitConsGRcompl <- fitConsGRcompl[strand(fitConsGRcompl) == "*"]

## set to NA genome positions without fitCons scores
fitConsGRcompl$score <- NA <- real <- 

## put together both, positions with and without fitCons scores to have every
## nucleotide covered with some value
fitConsGR <- sort(c(fitConsGR, fitConsGRcompl))

## transform BIGWIG into Rle objects coercing phastCons scores into
## 1-decimal digit raw-encoded values to reduce memory requirements
## in principle deciles of fitCons probabilities should give the
## necessary resolution for the purpose of filtering genetic variants

## round fitCons scores to one decimal digit, multiply them by 10 and
## coerce them into integer values
fitConsGR$score <- as.integer(10*round(fitConsGR$score, digits=1))

## because integer-coerced fitCons scores are later coerced into raw,
## and that coercion converts NA values into 0s, we set NA values to 255,
## so we can distinguish zero values from absence of fitCons scores (NAs) later on
fitConsGR$score[is.na(fitConsGR$score)] <- 255L

## store the fitCons scores from the GRanges object in a RleList object
fitConsRle <- coverage(fitConsGR, weight=fitConsGR$score)

stopifnot(sapply(fitConsRle, length) == seqlengths(fitConsGR)) ## QC (check that the Rle objects have the corresponding chromosome lengths)

## coerce the integer values of the Rle objects into 1-byte raw values
fitConsRle <- RleList(lapply(fitConsRle,
                             function(x) {
                               runValue(x) <- as.raw(runValue(x))
                               x
                             }),
                      compress=FALSE)

## save each raw-Rle object separately
for (chr in names(fitConsRle)) {
  tryCatch({
    metadata(fitConsRle[[chr]]) <- list(seqname=chr,
                                        provider="UCSC",
                                        provider_version="21Aug2014", ## it'd better to grab the date from downloaded file
                                        download_url=downloadURL,
                                        download_date=format(Sys.Date(), "%b %d, %Y"),
                                        reference_genome=refgenomeGD,
                                        data_pkgname="fitCons.UCSC.hg19")
    fname <- sprintf("fitCons.UCSC.hg19.%s.rds", chr)
    saveRDS(fitConsRle[[chr]], file=fname)
  }, error=function(err) {
      message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
