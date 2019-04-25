## This README file explains the steps to download and store in this
## annotation package the allele frequencies of the Freeze5 version of
## NHLBI TOPMED project. If you use these data please include the following
## citation:

## Taliun et al. Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program.
## bioRxiv, doi:10.1101/563866, 2019.

## The data were downloaded on March 2019 from
## https://bravo.sph.umich.edu/freeze5/hg38/download/all

## The data were first splitted into tabix VCF files per chromosome as follows:

## mkdir -p topmed_by_chr
## allchr=`tabix -l ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz`
## for chr in $allchr ; do {
##   echo $chr
##   tabix -h ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz $chr | bgzip -c > topmed_by_chr/$chr.vcf.gz
## 
##   if [ -s topmed_by_chr/$chr.vcf.gz ] ; then
##     tabix -p vcf topmed_by_chr/$chr.vcf.gz
##   else
##     rm topmed_by_chr/$chr.vcf.gz
##   fi
## } done

## Because we are going to lift GRCh37 coordinates over to GRCh38 coordinates
## we need to download first the corresponding UCSC chain file
## download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
##               "hg38ToHg19.over.chain.gz", method="curl")

## The following R script processes the downloaded and splitted data
## to transform the allele frequencies into minor allele frequencies
## and store them using only one significant digit for values < 0.1,
## and two significant digits for values > 0.1, to reduce the memory
## footprint through RleList objects

library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19) ## this is not the assembly used
                                     ## by TOPMed but is the assembly
                                     ## to which hg38 coordinates will
                                     ## be lifted
library(rtracklayer)
library(doParallel)

downloadURL <- "https://bravo.sph.umich.edu/freeze5/hg38/download/all"
datacitation <- bibentry(bibtype="Article",
                         author=c(person("Daniel Taliun"),
                                  person("et al.")),
                         title="Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program",
                         journal="bioRxiv",
                         year="2019",
                         doi="10.1101/563866",
                         note="Allele frequency data accessed on Mar. 2018",
                         url="http://bravo.sph.umich.edu")

registerDoParallel(cores=2)

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

vcfFilename <- "ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz"
genomeversion <- "hg19"
pkgname <- sprintf("MafDb.TOPMed.freeze5.%s", genomeversion)
dir.create(pkgname)

vcfHeader <- scanVcfHeader(vcfFilename)
## no genome reference information in the VCF header
## stopifnot(all(seqlengths(vcfHeader)[1:25] == seqlengths(Hsapiens)[1:25])) ## QC

## save the GenomeDescription object
refgenomeGD <- GenomeDescription(organism=organism(Hsapiens),
                                 common_name=commonName(Hsapiens),
                                 provider=provider(Hsapiens),
                                 provider_version=providerVersion(Hsapiens),
                                 release_date=releaseDate(Hsapiens),
                                 release_name=releaseName(Hsapiens),
                                 seqinfo=seqinfo(Hsapiens))
saveRDS(refgenomeGD, file=file.path(pkgname, "refgenomeGD.rds"))

## read INFO column data
infoCols <- rownames(info(vcfHeader))
ANcols <- "AN"
ACcols <- "AC"
stopifnot(all(c(ANcols, ACcols) %in% infoCols)) ## QC

message("Starting to process variants")

## restrict VCF INFO columns to AC and AN values
vcfPar <- ScanVcfParam(geno=NA,
                       fixed=c("ALT", "FILTER"),
                       info=c(ACcols, ANcols))

## fetch sequence names
tbx <- open(TabixFile(vcfFilename))
tbxchr <- sortSeqlevels(seqnamesTabix(tbx))
close(tbx)

hg38tohg19chain <- import.chain("hg38ToHg19.over.chain")

nsites <- foreach (chr=tbxchr, .combine='c') %dopar% {

  ## the whole VCF for the chromosome into main memory
  vcf <- readVcf(sprintf("topmed_by_chr/%s.vcf.gz", chr), genome="hg38", param=vcfPar)

  ## override SeqInfo data because chromosomal positions
  ## in the VCF are hg38 but we are going to lift them to hg19
  seqinfo(vcf, new2old=match(seqlevels(Hsapiens), seqlevels(vcf)),
          pruning.mode="coarse") <- seqinfo(Hsapiens)

  ## discard variants not passing all FILTERS
  mask <- fixed(vcf)$FILTER == "PASS"
  vcf <- vcf[mask, ]
  gc()

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

  ## lift over coordinates to hg19
  ## we exclude nonuniquely mapping positions
  ## and positions mapping to a different chromosome
  rr2 <- liftOver(rr, hg38tohg19chain)
  lomask <- elementNROWS(rr2) == 1
  chrmask <- rep(FALSE, length(rr))
  chrmask[lomask] <- as.character(seqnames(unlist(rr2[lomask]))) == as.character(seqnames(rr[lomask]))
  lomask <- lomask & chrmask
  stopifnot(any(lomask)) ## QC
  rr <- dropSeqlevels(unlist(rr2[lomask]), setdiff(seqlevels(rr2), seqlevels(rr)))

  ## clean up the ranges
  mcols(rr) <- NULL
  names(rr) <- NULL
  gc()

  nsites <- as.numeric(sum(lomask))

  ## according to https://samtools.github.io/hts-specs/VCFv4.3.pdf
  ## "It is permitted to have multiple records with the same POS"
  ## in such a case we take the maximum MAF by looking at repeated positions
  rrbypos <- split(rr, start(rr))
  rr <- rr[!duplicated(rr)]

  ## fetch allele frequency data
  acanValues <- info(vcfsnvs)[lomask, ]
  clsValues <- sapply(acanValues, class)

  for (j in seq_along(ACcols)) {
    acCol <- ACcols[j]
    anCol <- ANcols[j]
    message(sprintf("Processing SNVs allele frequencies from chromosome %s", chr))
    acValuesCol <- acanValues[[acCol]]
    anValuesCol <- acanValues[[anCol]]
    if (clsValues[acCol] == "CompressedIntegerList") {    ## in multiallelic variants take the
      acValuesCol <- as.numeric(sapply(acValuesCol, max)) ## maximum allele count
    }

    mafValuesCol <- acValuesCol / anValuesCol

    ## allele frequencies from TOPMED are calculated from alternative alleles,
    ## so for some of them we need to turn them into minor allele frequencies (MAF)
    mask <- !is.na(mafValuesCol) & mafValuesCol > 0.5
    if (any(mask))
      mafValuesCol[mask] <- 1 - mafValuesCol[mask]

    mafValuesCol <- relist(mafValuesCol, rrbypos)
    mafValuesCol <- sapply(mafValuesCol, max) ## in multiallelic variants
                                              ## take the maximum allele frequency

    q <- .quantizer(mafValuesCol)
    x <- .dequantizer(q)
    f <- cut(x, breaks=c(0, 10^c(seq(floor(min(log10(x[x!=0]), na.rm=TRUE)),
                                     ceiling(max(log10(x[x!=0]), na.rm=TRUE)), by=1))),
             include.lowest=TRUE)
    err <- abs(mafValuesCol-x)
    max.abs.error <- tapply(err, f, mean, na.rm=TRUE)

    ## build an integer-Rle object using the 'coverage()' function
    obj <- coverage(rr, weight=q)[[chr]]
    if (length(unique(mafValuesCol)) <= 10000) {
      Fn <- ecdf(mafValuesCol)
    } else {
      Fn <- ecdf(sample(mafValuesCol, size=10000, replace=TRUE))
    }

    ## coerce to raw-Rle, add metadata and save
     if (any(runValue(obj) != 0)) {
       runValue(obj) <- as.raw(runValue(obj))
       metadata(obj) <- list(seqname=chr,
                             provider="NHLBI",
                             provider_version="Freeze5",
                             citation=datacitation,
                             download_url=downloadURL,
                             download_date=format(Sys.Date(), "%b %d, %Y"),
                             reference_genome=refgenomeGD,
                             data_pkgname=pkgname,
                             qfun=.quantizer,
                             dqfun=.dequantizer,
                             ecdf=Fn,
                             max_abs_error=max.abs.error)
       saveRDS(obj, file=file.path(pkgname, sprintf("%s.AF.%s.rds", pkgname, chr)))
     } else {
       warning(sprintf("No MAF values for SNVs in chromosome %s", chr))
     }
  }

  ##
  ## nonSNVs
  ##

  ## fetch nonSNVs coordinates
  rr <- rowRanges(vcfnonsnvs)

  ## lift over coordinates to hg38
  ## we exclude nonuniquely mapping positions
  ## and positions mapping to a different chromosome
  rr2 <- liftOver(rr, hg38tohg19chain)
  lomask <- elementNROWS(rr2) == 1
  chrmask <- rep(FALSE, length(rr))
  chrmask[lomask] <- as.character(seqnames(unlist(rr2[lomask]))) == as.character(seqnames(rr[lomask]))
  lomask <- lomask & chrmask
  stopifnot(any(lomask)) ## QC
  rr <- dropSeqlevels(unlist(rr2[lomask]), setdiff(seqlevels(rr2), seqlevels(rr)))

  ## clean up the GRanges and save it
  mcols(rr) <- NULL
  names(rr) <- NULL

  nsites <- nsites + as.numeric(sum(lomask))

  ## according to https://samtools.github.io/hts-specs/VCFv4.3.pdf
  ## "It is permitted to have multiple records with the same POS"
  ## in such a case we take the maximum MAF by looking at repeated position
  rrbypos <- split(rr, paste(start(rr), end(rr), sep="-"))
  rr <- rr[!duplicated(rr)]
  saveRDS(rr, file=file.path(pkgname, sprintf("%s.GRnonsnv.%s.rds", pkgname, chr)))

  ## fetch allele frequency data
  acanValues <- info(vcfnonsnvs)[lomask, ]
  clsValues <- sapply(acanValues, class)

  rm(vcf)
  rm(vcfsnvs)
  rm(vcfnonsnvs)
  gc()

  for (j in seq_along(ACcols)) {
    acCol <- ACcols[j]
    anCol <- ANcols[j]
    message(sprintf("Processing nonSNVs allele frequencies from chromosome %s", chr))
    acValuesCol <- acanValues[[acCol]]
    anValuesCol <- acanValues[[anCol]]
    if (clsValues[acCol] == "CompressedIntegerList") {    ## in multiallelic variants take the
      acValuesCol <- as.numeric(sapply(acValuesCol, max)) ## maximum allele count
    }

    mafValuesCol <- acValuesCol / anValuesCol

    ## allele frequencies from TOPMED are calculated from alternative alleles,
    ## so for some of them we need to turn them into minor allele frequencies (MAF)
    mask <- !is.na(mafValuesCol) & mafValuesCol > 0.5
    if (any(mask))
      mafValuesCol[mask] <- 1 - mafValuesCol[mask]

    mafValuesCol <- relist(mafValuesCol, rrbypos)
    mafValuesCol <- sapply(mafValuesCol, max) ## in multiallelic variants
                                              ## take the maximum allele frequency
    q <- .quantizer(mafValuesCol)
    x <- .dequantizer(q)
    f <- cut(x, breaks=c(0, 10^c(seq(floor(min(log10(x[x!=0]), na.rm=TRUE)),
                                     ceiling(max(log10(x[x!=0]), na.rm=TRUE)), by=1))),
             include.lowest=TRUE)
    err <- abs(mafValuesCol-x)
    max.abs.error <- tapply(err, f, mean, na.rm=TRUE)

    if (length(unique(mafValuesCol)) <= 10000) {
      Fn <- ecdf(mafValuesCol)
    } else {
      Fn <- ecdf(sample(mafValuesCol, size=10000, replace=TRUE))
    }

    ## coerce the quantized value vector to an integer-Rle object
    obj <- Rle(q)

    ## coerce to raw-Rle, add metadata and save
    if (any(runValue(obj) != 0)) {
      runValue(obj) <- as.raw(runValue(obj))
      metadata(obj) <- list(seqname=chr,
                            provider="NHLBI",
                            provider_version="Freeze5",
                            citation=datacitation,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn,
                            max_abs_error=max.abs.error)
      saveRDS(obj, file=file.path(pkgname, sprintf("%s.RLEnonsnv.AF.%s.rds", pkgname, chr)))
    } else {
      warning(sprintf("No MAF values for nonSNVs in chromosome %s", chr))
    }
  }

  nsites
}

saveRDS(sum(nsites), file=file.path(pkgname, "nsites.rds"))
