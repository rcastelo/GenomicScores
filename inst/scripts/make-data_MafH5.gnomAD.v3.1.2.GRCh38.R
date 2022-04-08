#!/usr/bin/env Rscript

## This README file explains the steps to download and freeze in this
## annotation package the gnomAD genome allele frequencies. If you use
## these data please cite the following publication:

## Karczewski et al. The mutational constraint spectrum quantified from
## variation in 141,456 humans. Nature, 581:434-443, 2020.
## doi: https://doi.org/10.1038/s41586-020-2308-7

## See https://gnomad.broadinstitute.org/terms for further information
## on using this data for your own research

## The data were downloaded from as follows:
##
## allchr=`seq 1 22`" X Y"
## for chr in $allchr ; do {
##   wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr$chr.vcf.bgz
##   wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr$chr.vcf.bgz.tbi
## } done

## This script takes as argument a gnomAD v3.1.2 VCF file
## and process its alternate allele frequency (AF) values
## to calculate minor allele frequencies (MAF) from the
## global population and the maximum MAF throughout all
## other populations, and store their quantized compressed
## version into RDS files

stopifnot(BiocManager::version() == "3.15")

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop("make-data_MafH5.gnomAD.v3.1.2.GRCh38.R <VCFFILE>")
}

vcfFilename <- args[1]
if (!file.exists(vcfFilename))
  stop(sprintf("%s does not exist.", vcfFilename))

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))

genomeversion <- "GRCh38"
pkgname <- sprintf("MafH5.gnomAD.v3.1.2.%s", genomeversion)
dir.create(pkgname)

vcfHeader <- scanVcfHeader(vcfFilename)

namesstdchr <- standardChromosomes(Hsapiens)
commonseqnames <- intersect(seqnames(seqinfo(vcfHeader)), namesstdchr)
stopifnot(all(seqlengths(vcfHeader)[commonseqnames] == seqlengths(Hsapiens)[commonseqnames])) ## QC

downloadURL <- "https://gnomad.broadinstitute.org/downloads"
citationdata <- bibentry(bibtype="Article",
                         author=c(person("Konrad J Karczewski"), person("et al.")),
                         title="The mutational constraint spectrum quantified from variation in 141,456 humans",
                         journal="Nature",
                         volume="581",
                         pages="434-443",
                         year="2020",
                         doi="https://doi.org/10.1038/s41586-020-2308-7")

## quantizer function. it maps input real-valued [0, 1] allele frequencies
## to positive integers [1, 255] so that each of them can be later
## coerced into a single byte (raw type). allele frequencies between 0.1 and 1.0
## are quantized using two significant digits, while allele frequencies between
## 0 and 0.1 are quantized using one significant digit. quantized values are
## add up one unit to make them positive and keep the value 0 to encode for NAs
## since there is no NA value in the raw class (Sec. 3.3.4 NA handling, R Language Definition)
.quantizer <- function(x) {
  .ndec <- function(x) {
    spl <- strsplit(as.character(x+1), "\\.")
    spl <- sapply(spl, "[", 2)
    spl[is.na(spl)] <- ""
    nchar(spl)
  }

  maskNAs <- is.na(x)
  x[!maskNAs & x > 0.1] <- signif(x[!maskNAs & x > 0.1], digits=2)
  x[!maskNAs & x <= 0.1] <- signif(x[!maskNAs & x <= 0.1], digits=1)
  x[maskNAs] <- NA
  nd <- .ndec(x)
  x[nd < 3] <- x[nd < 3] * 100
  x[nd >= 3] <- x[nd >= 3] * 10^nd[nd >= 3] + 100 + (nd[nd >= 3] - 3) * 10
  x <- x + 1
  x[maskNAs] <- 0
  q <- as.integer(sprintf("%.0f", x)) ## coercing through character is necessary
                                      ## to deal with limitations of computer
                                      ## arithmetic such as with as.integer(0.29*100)
  if (any(q > 255))
    stop("current number of quantized values > 255 and cannot be stored into one byte")
  q
}
attr(.quantizer, "description") <- "quantize [0.1-1] with 2 significant digits and [0-0.1) with 1 significant digit"

## dequantizer function
.dequantizer <- function(q) {
  x <- q <- as.integer(q)
  x[x == 0L] <- NA
  x <- (x - 1L) / 100L
  maskNAs <- is.na(x)
  q <- q - 1L
  sel <- !maskNAs & q > 100L
  x[sel] <- ((q[sel] - 100L) %% 10L) / 10^(floor((q[sel] - 100L) / 10L) + 3L)
  x
}
attr(.dequantizer, "description") <- "dequantize [0-100] dividing by 100, [101-255] subtract 100, take modulus 10 and divide by the corresponding power in base 10"

## fill up SeqInfo data
si <- keepStandardChromosomes(seqinfo(vcfHeader))
isCircular(si) <- isCircular(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
genome(si) <- genome(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]

## save the GenomeDescription object
refgenomeGD <- GenomeDescription(organism=organism(Hsapiens),
                                 common_name=commonName(Hsapiens),
                                 provider=provider(Hsapiens),
                                 provider_version=metadata(Hsapiens)$genome,
                                 release_date=releaseDate(Hsapiens),
                                 release_name=metadata(Hsapiens)$genome,
                                 seqinfo=si)

## read INFO column data
infoCols <- rownames(info(vcfHeader))
AFcols <- infoCols[grep("^AF", infoCols)]
exclude <- grep("raw|male|max", AFcols)
if (length(exclude) > 0)
  AFcols <- AFcols[-exclude] ## remove columns such as raw and those specific for sex and maximum allele frequency among populations

## actually restrict to global frequencies and global non-cancer frequencies
AFcols <- AFcols[nchar(AFcols) == 6 | nchar(AFcols) == 2] ## select global AF and for main populations
stopifnot(all(AFcols %in% infoCols)) ## QC


## restrict VCF INFO columns to AF values
vcfPar <- ScanVcfParam(geno=NA,
                       fixed=c("ALT", "FILTER"),
                       info=AFcols)

message(sprintf("Reading VCF file %s", basename(vcfFilename)))

## read the whole VCF into main memory (up to 32Gb of RAM for chromosome 1)
vcf <- readVcf(vcfFilename, genome=genomeversion, param=vcfPar)

vcf <- keepStandardChromosomes(vcf)
seqinfo(vcf) <- si

## discard variants not passing all FILTERS
mask <- fixed(vcf)$FILTER == "PASS"
vcf <- vcf[mask, ]
gc()

## save the total number of variants
saveRDS(nrow(vcf),
        file=file.path(pkgname, paste0(gsub(".vcf.bgz", "", basename(vcfFilename)), ".nsites.rds")))

## mask variants where all alternate alleles are SNVs
evcf <- expand(vcf)
maskSNVs <- sapply(relist(isSNV(evcf), alt(vcf)), all)
rm(evcf)
gc()

## treat snvs and nonSNVs separately
vcfsnvs <- vcf[maskSNVs, ]
vcfnonsnvs <- vcf[!maskSNVs, ]

##
## SNVs
##

## fetch SNVs coordinates
rr <- rowRanges(vcfsnvs)

## fetch current chromosome
chr <- as.character(seqnames(rr)[1])

## clean up the ranges
mcols(rr) <- NULL
names(rr) <- NULL
gc()

## fill up missing SeqInfo data
isCircular(rr) <- isCircular(si)
seqlengths(rr) <- seqlengths(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
isCircular(rr) <- isCircular(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
genome(rr) <- genome(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]

## according to https://samtools.github.io/hts-specs/VCFv4.3.pdf
## "It is permitted to have multiple records with the same POS"
rrbypos <- split(rr, start(rr))
rr <- rr[!duplicated(rr)]
## put back the genomic order
mt <- match(as.character(start(rr)), names(rrbypos))
stopifnot(all(!is.na(mt))) ## QC
rrbypos <- rrbypos[mt]
rm(mt)
gc()

## fetch allele frequency data
afValues <- info(vcfsnvs)
clsValues <- sapply(afValues, class)
allpopmax <- rep(NA_real_, length(rr)) ## store maximum MAF among all populations

for (j in seq_along(AFcols)) {
  afCol <- AFcols[j]

  message(sprintf("Processing %s SNV allele frequencies from %s", afCol, basename(vcfFilename)))
  mafValuesCol <- afValues[[afCol]]
  if (clsValues[afCol] == "numeric" || clsValues[afCol] == "Numeric") {
    mafValuesCol <- as.numeric(mafValuesCol)
  } else if (clsValues[afCol] == "character") {
    mafValuesCol <- as.numeric(sub(pattern="\\.", replacement="0", mafValuesCol))
  } else if (clsValues[afCol] == "CompressedNumericList") { ## in multiallelic variants take
    mafValuesCol <- sapply(mafValuesCol, max)               ## the maximum allele frequency
  } else if (clsValues[afCol] == "CompressedCharacterList") {
    mafValuesCol <- sapply(NumericList(lapply(mafValuesCol, sub, pattern="\\.", replacement="0")), max)
  } else {
    stop(sprintf("%s: Uknown class for holding AF values (%s)", vcfFilename, clsValues[afCol]))
  }

  ## allele frequencies from gnomAD are calculated from alternative alleles,
  ## so for some of them we need to turn them into minor allele frequencies (MAF)
  ## for biallelic variants, in those cases the MAF comes from the REF allele
  maskREF <- !is.na(mafValuesCol) & mafValuesCol > 0.5
  if (any(maskREF))
    mafValuesCol[maskREF] <- 1 - mafValuesCol[maskREF]

  mafValuesCol <- relist(mafValuesCol, rrbypos)
  maskREF <- relist(maskREF, rrbypos)
  mafValuesCol <- sapply(mafValuesCol, max) ## in multiallelic variants
                                            ## take the maximum allele frequency
  maskREF <- sapply(maskREF, any) ## in multiallelic variants, when any of the
                                  ## alternate alleles has AF > 0.5, then we
                                  ## set to TRUE maskREF as if the MAF is in REF
  
  if (afCol != "AF") { ## for any population other than global, keep the maximum MAF
    allpopmax <- pmax(allpopmax, mafValuesCol, na.rm=TRUE)
  } else {
    q <- .quantizer(mafValuesCol)
    x <- .dequantizer(q)
    f <- cut(x, breaks=c(0, 10^c(seq(floor(min(log10(x[x!=0]), na.rm=TRUE)),
                                     ceiling(max(log10(x[x!=0]), na.rm=TRUE)), by=1))),
             include.lowest=TRUE)
    err <- abs(mafValuesCol-x)
    max.abs.error <- tapply(err, f, mean, na.rm=TRUE)

    ## build an integer-Rle object using the 'coverage()' function
    obj <- coverage(rr, weight=q)[[chr]]
    ## build an integer-Rle object of maskREF using the 'coverage()' function
    maskREFobj <- coverage(rr, weight=maskREF+0L)[[chr]]

    ## build ECDF of MAF values
    if (length(unique(mafValuesCol)) <= 10000) {
      Fn <- ecdf(mafValuesCol)
    } else {
      Fn <- ecdf(sample(mafValuesCol, size=10000, replace=TRUE))
    }

    ## coerce to raw-Rle, add metadata and save
    if (any(runValue(obj) != 0)) {
      runValue(obj) <- as.raw(runValue(obj))
      runValue(maskREFobj) <- as.raw(runValue(maskREFobj))
      metadata(obj) <- list(seqname=chr,
                            provider="BroadInstitute",
                            provider_version="v3.1.2",
                            citation=citationdata,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn,
                            max_abs_error=max.abs.error,
                            maskREF=maskREFobj)
      saveRDS(obj, file=file.path(pkgname, sprintf("%s.%s.%s.rds", pkgname, afCol, chr)))
    } else {
      warning(sprintf("No MAF values for SNVs in %s", vcfFilename))
    }
  }
}

## store allpopmax "population" values
q <- .quantizer(allpopmax)
x <- .dequantizer(q)
f <- cut(x, breaks=c(0, 10^c(seq(floor(min(log10(x[x!=0]), na.rm=TRUE)),
                                 ceiling(max(log10(x[x!=0]), na.rm=TRUE)), by=1))),
         include.lowest=TRUE)
err <- abs(allpopmax-x)
max.abs.error <- tapply(err, f, mean, na.rm=TRUE)

## build an integer-Rle object using the 'coverage()' function
obj <- coverage(rr, weight=q)[[chr]]
## build an integer-Rle object of maskREF using the 'coverage()' function
maskREFobj <- coverage(rr, weight=maskREF+0L)[[chr]]

## build ECDF of MAF values
if (length(unique(allpopmax)) <= 10000) {
  Fn <- ecdf(allpopmax)
} else {
  Fn <- ecdf(sample(allpopmax, size=10000, replace=TRUE))
}

## coerce to raw-Rle, add metadata and save
if (any(runValue(obj) != 0)) {
  runValue(obj) <- as.raw(runValue(obj))
  runValue(maskREFobj) <- as.raw(runValue(maskREFobj))
  metadata(obj) <- list(seqname=chr,
                        provider="BroadInstitute",
                        provider_version="v3.1.2",
                        citation=citationdata,
                        download_url=downloadURL,
                        download_date=format(Sys.Date(), "%b %d, %Y"),
                        reference_genome=refgenomeGD,
                        data_pkgname=pkgname,
                        qfun=.quantizer,
                        dqfun=.dequantizer,
                        ecdf=Fn,
                        max_abs_error=max.abs.error,
                        maskREF=maskREFobj)
  saveRDS(obj, file=file.path(pkgname, sprintf("%s.AF_allpopmax.%s.rds", pkgname, chr)))
} else {
  warning(sprintf("No MAF values for SNVs in %s", vcfFilename))
}

rm(vcfsnvs)
gc()

##
## nonSNVs
##

## fetch nonSNVs coordinates
rr <- rowRanges(vcfnonsnvs)

## fill up missing SeqInfo data
si <- seqinfo(vcf)
seqlengths(rr) <- seqlengths(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
isCircular(rr) <- isCircular(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
genome(rr) <- genome(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]

## clean up the GRanges
mcols(rr) <- NULL
names(rr) <- NULL
gc()

## re-order by chromosomal coordinates to deal with wrongly-ordered nonSNVs over multiple VCF lines
ord <- order(rr)
rr <- rr[ord]

## fetch allele frequency data
afValues <- info(vcfnonsnvs)
clsValues <- sapply(afValues, class)
rm(vcf)
rm(vcfnonsnvs)
gc()

## re-order by chromosomal coordinates
afValues <- afValues[ord, ]
rm(ord)

## according to https://samtools.github.io/hts-specs/VCFv4.3.pdf
## "It is permitted to have multiple records with the same POS"
posids <- paste(start(rr), end(rr), sep="-")
rrbypos <- split(rr, posids)
rr <- rr[!duplicated(rr)]
## put back the genomic order
posids <- paste(start(rr), end(rr), sep="-")
mt <- match(posids, names(rrbypos))
stopifnot(all(!is.na(mt))) ## QC
rrbypos <- rrbypos[mt]
saveRDS(rr, file=file.path(pkgname, sprintf("%s.GRnonsnv.%s.rds", pkgname, chr)))

rm(posids)
rm(mt)
gc()

allpopmax <- rep(NA_real_, length(rr)) ## store maximum MAF among all populations

for (j in seq_along(AFcols)) {
  afCol <- AFcols[j]

  message(sprintf("Processing %s nonSNV allele frequencies from %s", afCol, basename(vcfFilename)))
  mafValuesCol <- afValues[[afCol]]
  if (clsValues[afCol] == "numeric" || clsValues[afCol] == "Numeric") {
    mafValuesCol <- as.numeric(mafValuesCol)
  } else if (clsValues[afCol] == "character") {
    mafValuesCol <- as.numeric(sub(pattern="\\.", replacement="0", mafValuesCol))
  } else if (clsValues[afCol] == "CompressedNumericList") { ## in multiallelic variants take
    mafValuesCol <- sapply(mafValuesCol, max)               ## the maximum allele frequency
  } else if (clsValues[afCol] == "CompressedCharacterList") {
    mafValuesCol <- sapply(NumericList(lapply(mafValuesCol, sub, pattern="\\.", replacement="0")), max)
  } else {
    stop(sprintf("Uknown class for holding AF values (%s)", clsValues[afCol]))
  }

  ## allele frequencies from gnomAD genomes are calculated from alternative alleles,
  ## so for some of them we need to turn them into minor allele frequencies (MAF)
  maskREF <- !is.na(mafValuesCol) & mafValuesCol > 0.5
  if (any(maskREF))
    mafValuesCol[maskREF] <- 1 - mafValuesCol[maskREF]

  mafValuesCol <- relist(mafValuesCol, rrbypos)
  maskREF <- relist(maskREF, rrbypos)
  mafValuesCol <- sapply(mafValuesCol, max) ## in multiallelic variants
                                            ## take the maximum allele frequency
  maskREF <- sapply(maskREF, any) ## in multiallelic variants, when any of the
                                  ## alternate alleles has AF > 0.5, then we
                                  ## set to TRUE maskREF as if the MAF is in REF

  if (afCol != "AF") { ## for any population other than global, keep the maximum MAF
    allpopmax <- pmax(allpopmax, mafValuesCol, na.rm=TRUE)
  } else {
    q <- .quantizer(mafValuesCol)
    x <- .dequantizer(q)
    f <- cut(x, breaks=c(0, 10^c(seq(floor(min(log10(x[x!=0]), na.rm=TRUE)),
                                     ceiling(max(log10(x[x!=0]), na.rm=TRUE)), by=1))),
             include.lowest=TRUE)
    err <- abs(mafValuesCol-x)
    max.abs.error <- tapply(err, f, mean, na.rm=TRUE)

    ## build ECDF of MAF values
    if (length(unique(mafValuesCol)) <= 10000) {
      Fn <- ecdf(mafValuesCol)
    } else {
      Fn <- ecdf(sample(mafValuesCol, size=10000, replace=TRUE))
    }

    ## coerce the quantized value vector to an integer-Rle object
    obj <- Rle(q)
    ## coerce the maskREF vector to an integer-Rle object
    maskREFobj <- Rle(maskREF+0L)

    ## coerce to raw-Rle, add metadata and save
    if (any(runValue(obj) != 0)) {
      runValue(obj) <- as.raw(runValue(obj))
      runValue(maskREFobj) <- as.raw(runValue(maskREFobj))
      metadata(obj) <- list(seqname=chr,
                            provider="BroadInstitute",
                            provider_version="v3.1.2",
                            citation=citationdata,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn,
                            max_abs_error=max.abs.error,
                            maskREF=maskREFobj)
      saveRDS(obj, file=file.path(pkgname, sprintf("%s.RLEnonsnv.%s.%s.rds", pkgname, afCol, chr)))
    } else {
      warning(sprintf("No MAF values for nonSNVs in %s", basename(vcfFilename)))
    }
  }
}

## store allpopmax "population" values
q <- .quantizer(allpopmax)
x <- .dequantizer(q)
f <- cut(x, breaks=c(0, 10^c(seq(floor(min(log10(x[x!=0]), na.rm=TRUE)),
                                 ceiling(max(log10(x[x!=0]), na.rm=TRUE)), by=1))),
         include.lowest=TRUE)
err <- abs(allpopmax-x)
max.abs.error <- tapply(err, f, mean, na.rm=TRUE)

## build ECDF of MAF values
if (length(unique(allpopmax)) <= 10000) {
  Fn <- ecdf(allpopmax)
} else {
  Fn <- ecdf(sample(allpopmax, size=10000, replace=TRUE))
}

## coerce the quantized value vector to an integer-Rle object
obj <- Rle(q)
## coerce the maskREF vector to an integer-Rle object
maskREFobj <- Rle(maskREF+0L)

## coerce to raw-Rle, add metadata and save
if (any(runValue(obj) != 0)) {
  runValue(obj) <- as.raw(runValue(obj))
  runValue(maskREFobj) <- as.raw(runValue(maskREFobj))
  metadata(obj) <- list(seqname=chr,
                        provider="BroadInstitute",
                        provider_version="v3.1.2",
                        citation=citationdata,
                        download_url=downloadURL,
                        download_date=format(Sys.Date(), "%b %d, %Y"),
                        reference_genome=refgenomeGD,
                        data_pkgname=pkgname,
                        qfun=.quantizer,
                        dqfun=.dequantizer,
                        ecdf=Fn,
                        max_abs_error=max.abs.error,
                        maskREF=maskREFobj)
  saveRDS(obj, file=file.path(pkgname, sprintf("%s.RLEnonsnv.AF_allpopmax.%s.rds", pkgname, chr)))
} else {
  warning(sprintf("No MAF values for nonSNVs in %s", basename(vcfFilename)))
}

message(sprintf("Finished processing VCF file %s", basename(vcfFilename)))


## once the script above has been run for all VCF files, one
## should execute the following line to build HDF5 files from
## the previously generated RDS files, where 'pkgname' is the
## name of the directory storing the generated RDS files (i.e.,
## the variable 'pkgname' defined above).
##
## GenomicScores:::.makeMafH5(pkgname, pkgname, pkgname)
