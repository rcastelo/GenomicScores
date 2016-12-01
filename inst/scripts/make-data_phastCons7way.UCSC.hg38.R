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

foreach (chr=seqnames(Hsapiens)) %dopar% {
  cat(chr, "\n")
  tryCatch({
    assign(sprintf("phastCons7way_%s", chr),
           10*round(import.bw(BigWigFile(file.path("hg38.7way.phastCons", "hg38.phastCons7way.bw")),
                              as="RleList", which=GRanges(seqnames=chr, IRanges(1, seqlengths(Hsapiens)[chr])))[[chr]], digits=1))
    assign(sprintf("phastCons7way_%s", chr),
           do.call("runValue<-", list(get(sprintf("phastCons7way_%s", chr)),
                                      as.raw(runValue(get(sprintf("phastCons7way_%s", chr)))))))
    metadata(get(sprintf("phastCons7way_%s", chr))) <- list(seqname=chr,
                                                            provider="UCSC",
                                                            provider_version="04Jun2014", ## it'd better to grab the date from downloaded file
                                                            download_url=downloadURL,
                                                            download_date=format(Sys.Date(), "%b %d, %Y"),
                                                            reference_genome=refgenomeGD,
                                                            data_pkgname="phastCons7way.UCSC.hg38")
    saveRDS(get(sprintf("phastCons7way_%s", chr)),
            file=file.path("hg38.7way.phastCons", sprintf("phastCons7way.UCSC.hg38.%s.rds", chr)))
    rm(list=sprintf("phastCons7way_%s", chr))
    gc()
  }, error=function(err) {
      message(chr, " ", conditionMessage(err), call.=TRUE)
  })
}
