## This README file explains the steps to download and freeze in this
## annotation package the allele frequencies of the Phase 1 of the
## 1000 Genomes Project. If you use these data please cite the following publication:

## The 1000 Genomes Project Consortium. An integrated map of genetic variation
## from 1,092 human genomes. Nature, 491:56-65, 2012.
## doi: http://dx.doi.org/10.1038/nature11632

## The data were downloaded from the FTP server of the 1000 Genomes Project as follows:
##
## wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz
## wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz.tbi

## The data were first splitted into tabix VCF files per chromosome as follows:
##
## mkdir -p phase1_by_chr
## allchr=`tabix -l ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz`
## for chr in $allchr ; do {
##   echo chr$chr
##   tabix -h ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz $chr | bgzip -c > phase1_by_chr/chr$chr.vcf.gz
##
##   if [ -s phase1_by_chr/chr$chr.vcf.gz ] ; then
##     tabix -p vcf phase1_by_chr/chr$chr.vcf.gz
##   else
##     rm phase1_by_chr/chr$chr.vcf.gz
##   fi
## } done

## Because we are going to lift GRCh37 coordinates over to GRCh38 coordinates
## we need to download first the corresponding UCSC chain file
## download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
##               "hg19ToHg38.over.chain.gz", method="curl")

## The following R script processes the downloaded and splitted data
## to transform the allele frequencies into minor allele frequencies
## and store them using only one significant digit for values < 0.1,
## and two significant digits for values > 0.1, to reduce the memory
## footprint through RleList objects

library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(VariantAnnotation)
library(BSgenome.Hsapiens.NCBI.GRCh38) ## this is not the assembly used by
                                       ## the 1000 Genomes Project but is
                                       ## the assembly to which GRCh37
                                       ## coordinates will be lifted
library(rtracklayer)
library(doParallel)

downloadURL <- "http://www.internationalgenome.org/data"
citationdata <- bibentry(bibtype="Article",
                         author=person("The 1000 Genomes Project Consortium"),
                         title="An integrated map of genetic variation from 1,092 human genomes",
                         journal="Nature",
                         volume="491",
                         pages="56-65",
                         year="2012",
                         doi="10.1038/nature11632")

registerDoParallel(cores=8)

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

vcfFilename <- "ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz"
genomeversion <- "GRCh38"
pkgname <- sprintf("MafDb.1Kgenomes.phase1.%s", genomeversion)
dir.create(pkgname)

path2vcfs <- "/projects_fg/GenomicScores/1000G/Phase1"

vcfHeader <- scanVcfHeader(file.path(path2vcfs, vcfFilename))

## no seqinfo and genome reference information in the VCF header
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
AFcols <- infoCols[c(which(infoCols == "AF"), grep("_AF", infoCols))]

message("Starting to process variants")

## restrict VCF INFO columns to AC and AN values
vcfPar <- ScanVcfParam(geno=NA,
                       fixed=c("ALT", "FILTER"),
                       info=AFcols)

tbx <- open(TabixFile(file.path(path2vcfs, vcfFilename)))
tbxchr <- sortSeqlevels(seqnamesTabix(tbx))
close(tbx)

hg19tohg38chain <- import.chain("hg19ToHg38.over.chain")

foreach (chr=tbxchr) %dopar% {

  ## read the whole VCF for the chromosome into main memory
  vcf <- readVcf(sprintf("%s/phase1_by_chr/chr%s.vcf.gz", path2vcfs, chr),
                 genome="hs37d5", param=vcfPar)

  ## discard variants not passing all FILTERS
  mask <- fixed(vcf)$FILTER == "PASS"
  if (any(mask)) {
    vcf <- vcf[mask, ]
  } else {
    stop("No variants with FILTER=PASS")
  }
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

  ## lift over coordinates to hg38
  ## we exclude nonuniquely mapping positions
  ## and positions mapping to a different chromosome
  seqlevelsStyle(rr) <- "UCSC"
  rr2 <- liftOver(rr, hg19tohg38chain)
  lomask <- elementNROWS(rr2) == 1
  chrmask <- rep(FALSE, length(rr))
  chrmask[lomask] <- as.character(seqnames(unlist(rr2[lomask]))) == as.character(seqnames(rr[lomask]))
  lomask <- lomask & chrmask
  stopifnot(any(lomask)) ## QC
  rr <- dropSeqlevels(unlist(rr2[lomask]), setdiff(seqlevels(rr2), seqlevels(rr)))
  seqlevelsStyle(rr) <- "NCBI"

  ## clean up the ranges
  mcols(rr) <- NULL
  names(rr) <- NULL
  gc()

  ## according to https://samtools.github.io/hts-specs/VCFv4.3.pdf
  ## "It is permitted to have multiple records with the same POS"
  posids <- paste(seqnames(rr), start(rr), sep="-")
  rrbypos <- split(rr, posids)
  rr <- rr[!duplicated(rr)]
  ## put back the genomic order
  posids <- paste(seqnames(rr), start(rr), sep="-")
  mt <- match(posids, names(rrbypos))
  stopifnot(all(!is.na(mt))) ## QC
  rrbypos <- rrbypos[mt]
  rm(posids)
  rm(mt)
  gc()

  ## fill up SeqInfo data
  si <- seqinfo(vcf)
  seqlengths(rr) <- seqlengths(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
  isCircular(rr) <- isCircular(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
  genome(rr) <- genome(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]

  ## fetch allele frequency data
  afValues <- info(vcfsnvs)[lomask, ]
  clsValues <- sapply(afValues, class)

  for (j in seq_along(AFcols)) {
    afCol <- AFcols[j]

    message(sprintf("Processing %s SNVs from chromosome %s", afCol, chr))
    mafValuesCol <- afValues[[afCol]]
    if (clsValues[afCol] == "numeric" || clsValues[afCol] == "Numeric") {
      mafValuesCol <- as.numeric(mafValuesCol)
    } else if (clsValues[afCol] == "CompressedNumericList") { ## in multiallelic variants take
      mafValuesCol <- sapply(mafValuesCol, max)               ## the maximum allele frequency
    } else {
      stop(sprintf("Uknown class for holding AF values (%s)", clsValues[afCol]))
    }

    ## allele frequencies from 1000 genomes are calculated from alternative alleles,
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
                            provider="IGSR",
                            provider_version="Phase1",
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
  seqlevelsStyle(rr) <- "UCSC"
  rr2 <- liftOver(rr, hg19tohg38chain)
  lomask <- elementNROWS(rr2) == 1
  chrmask <- rep(FALSE, length(rr))
  chrmask[lomask] <- as.character(seqnames(unlist(rr2[lomask]))) == as.character(seqnames(rr[lomask]))
  lomask <- lomask & chrmask
  stopifnot(any(lomask)) ## QC
  rr <- dropSeqlevels(unlist(rr2[lomask]), setdiff(seqlevels(rr2), seqlevels(rr)))
  seqlevelsStyle(rr) <- "NCBI"

  ## fill up SeqInfo data
  si <- seqinfo(vcf)
  seqlengths(rr) <- seqlengths(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
  isCircular(rr) <- isCircular(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
  genome(rr) <- genome(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]

  ## clean up the GRanges and save it
  mcols(rr) <- NULL
  names(rr) <- NULL
  gc()

  ## re-order by chromosomal coordinates to deal with wrongly-ordered nonSNVs over multiple VCF lines
  ord <- order(rr)
  rr <- rr[ord]

  ## fetch allele frequency data
  afValues <- info(vcfnonsnvs)[lomask, ]
  clsValues <- sapply(afValues, class)
  rm(vcf)
  rm(vcfnonsnvs)
  gc()

  ## re-order by chromosomal coordinates
  afValues <- afValues[ord, ]

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

  rm(ord)
  rm(posids)
  rm(mt)
  gc()

  for (j in seq_along(AFcols)) {
    afCol <- AFcols[j]

    message(sprintf("Processing %s nonSNVs from chromosome %s", afCol, chr))
    mafValuesCol <- afValues[[afCol]]
    if (clsValues[afCol] == "numeric" || clsValues[afCol] == "Numeric") {
      mafValuesCol <- as.numeric(mafValuesCol)
    } else if (clsValues[afCol] == "CompressedNumericList") { ## in multiallelic variants take
      mafValuesCol <- sapply(mafValuesCol, max)               ## the maximum allele frequency
    } else {
      stop(sprintf("Uknown class for holding AF values (%s)", clsValues[afCol]))
    }

    ## allele frequencies from 1000 genomes are calculated from alternative alleles,
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
                            provider="IGSR",
                            provider_version="Phase1",
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
      warning(sprintf("No MAF values for nonSNVs in chromosome %s", chr))
    }
  }
}

## save rsIDs assignments from the 1000 genomes project
## streaming through the whole file
vcfPar <- ScanVcfParam(geno=NA,
                       fixed=c("ALT", "FILTER"),
                       info=NA)

tbx <- TabixFile(file.path(path2vcfs, vcfFilename), yieldSize=1000000)
open(tbx)

message("Starting to process variant identifiers")

rsIDs <- character(0)  ## to store rsIDs annotated by the 1000 genomes project
rsIDgp <- GPos()       ## to store positions of rsIDs
maskSNVs <- logical(0) ## to store a mask whether the variant is an SNV or not

nVar <- nremVar <- 0
while (nrow(vcf <- readVcf(tbx, genome="hs37d5", param=vcfPar))) {
  
  ## discard variants not passing all FILTERS
  mask <- fixed(vcf)$FILTER == "PASS"
  if (any(mask)) {
    vcf <- vcf[mask, ]
  } else {
    stop("No variants with FILTER=PASS")
  }
  gc()

  nVar <- nVar + nrow(vcf)
  rr <- rowRanges(vcf)

  ## lift over coordinates to hg38
  ## we exclude nonuniquely mapping positions
  ## and positions mapping to a different chromosome
  seqlevelsStyle(rr) <- "UCSC"
  rr2 <- liftOver(rr, hg19tohg38chain)
  lomask <- elementNROWS(rr2) == 1
  chrmask <- rep(FALSE, length(rr))
  chrmask[lomask] <- as.character(seqnames(unlist(rr2[lomask]))) == as.character(seqnames(rr[lomask]))
  lomask <- lomask & chrmask
  stopifnot(any(lomask)) ## QC
  rr <- dropSeqlevels(unlist(rr2[lomask]), setdiff(seqlevels(rr2), seqlevels(rr)))
  seqlevelsStyle(rr) <- "NCBI"
  nremVar <- nremVar + sum(!lomask)
  vcf <- vcf[lomask, ]

  whrsIDs <- grep("^rs", names(rr))
  evcf <- expand(vcf)
  maskSNVs <- c(maskSNVs, sapply(relist(isSNV(evcf), alt(vcf)), all)[whrsIDs])
  rm(evcf)
  gc()
  rrTmp <- rr[whrsIDs]
  mcols(rrTmp) <- NULL
  rrTmp <- resize(rrTmp, width=1, fix="start", ignore.strand=TRUE)
  idTmp <- names(rrTmp)
  names(rrTmp) <- NULL
  gpTmp <- as(rrTmp, "GPos")

  rsIDs <- c(rsIDs, idTmp)
  rsIDgp <- c(rsIDgp, gpTmp)
  rm(rrTmp)
  rm(gpTmp)
  rm(idTmp)
  gc()

  message(sprintf("%d variant identifiers processed", nVar))
}
close(tbx)

## save total number of variants
saveRDS(nVar-nremVar, file=file.path(pkgname, "nsites.rds"))

## store mask flagging SNVs
rsIDgp$isSNV <- Rle(maskSNVs)

## double check that all identifiers start with 'rs'
stopifnot(identical(grep("^rs", rsIDs), 1:length(rsIDs))) ## QC

## chop the 'rs' prefix and convert the character ids into integer values
rsIDs <- as.integer(sub(pattern="^rs", replacement="", x=rsIDs))

## calculate the indices that lead to an increasing values of the (integer) ids
rsIDidx <- order(rsIDs)

## order increasingly the integer values of rs IDs
rsIDs <- rsIDs[rsIDidx]

## save the objects that enable the search for rs ID
saveRDS(rsIDs, file=file.path(pkgname, "rsIDs.rds"))
saveRDS(rsIDidx, file=file.path(pkgname, "rsIDidx.rds"))
saveRDS(rsIDgp, file=file.path(pkgname, "rsIDgp.rds"))
