## This README file explains the steps to download and freeze in this
## annotation package the allele frequencies of the NHLBI Exome Sequencing
## Project (ESP). If you use these data please cite the following publication:

## Tennessen JA, et al. Evolution and functional impact of rare coding
## variation from deep sequencing of human exomes. Science, 337:64-69, 2012.
## doi: http://dx.doi.org/10.1126/science.1219240

## The data were downloaded from the Exome Variant Server (EVS) as follows:
##
## wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz

## The following R script process the downloaded data to
## transform the allele frequencies into minor allele frequencies
## and store them using only one significant digit for values < 0.1,
## and two significant digits for values > 0.1, to reduce the memory
## footprint through RleList objects

library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(VariantAnnotation)
library(BSgenome.Hsapiens.NCBI.GRCh38) ## this is not the assembly used by
                                       ## the NHLBI ESP but is the assembly
                                       ## to which GRCh37 coordinates were
                                       ## lifted and stored in one of the
                                       ## info columns of the VCF file
library(doParallel)

downloadURL <- "http://evs.gs.washington.edu/EVS"
citationdata <- bibentry(bibtype="Article",
                         author=c(person("Jacob A. Tennessen"), person("Abigail W. Bigham"),
                                  person("Timothy D. O'Connor"), person("Wenqing Fu"),
                                  person("Eimear E. Kenny"), person("Simon Gravel"),
                                  person("Sean McGee"), person("Xiaoming Liu"),
                                  person("Goo Jun"), person("Hyun Min Kang"),
                                  person("Daniel Jordan"), person("Suzanne M. Leal"),
                                  person("Stacey Gabriel"), person("Mark J. Rieder"),
                                  person("Goncalo Abecasis"), person("David Altshuler"),
                                  person("Deborah A. Nickerson"), person("Eric Boerwinkle"),
                                  person("Shamil Sunyaev"), person("Carlos D. Bustamante"),
                                  person("Michael J. Bamshad"), person("Joshua M. Akey"),
                                  person("Broad GO"), person("Seattle GO")),
                         title="Evoluution and functional impact of rare coding variation from deep sequencing of human exomes",
                         journal="Science",
                         volume="337",
                         pages="64-69",
                         year="2012",
                         doi="10.1126/science.1219240")

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

tarballFilename <- "ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz"
vcfFilenames <- untar(tarballFilename, list=TRUE)

## double check that the extracted files are VCF files
for (f in vcfFilenames) {
  tryCatch({
    vcfHeader <- scanVcfHeader(f)
  }, error=function(err) {
    stop(sprintf("extracted VCF file %s has not a valid VCF header.", f), call.=TRUE)
  })
}

genomeversion <- "GRCh38"
pkgname <- sprintf("MafDb.ESP6500SI.V2.SSA137.%s", genomeversion)
dir.create(pkgname)

vcfHeader <- scanVcfHeader(vcfFilenames[1])

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

## build SeqInfo object from Hsapiens
mt <- regexpr(pattern="chr[0-9XY]+", vcfFilenames)
chrs <- sub("chr", "",
            substr(vcfFilenames, mt, mt+attr(mt, "match.length")-1))
si <- seqinfo(Hsapiens)[seqnames(Hsapiens)[seqnames(Hsapiens) %in% chrs]]

vcfPar <- ScanVcfParam(geno=NA,
                       fixed="ALT",
                       info=c("MAF", "GRCh38_POSITION"))

message("Starting to process variants")

foreach (fname=vcfFilenames) %dopar% {
  mt <- regexpr(pattern="chr[0-9XY]+", fname)
  chr <- sub("chr", "", substr(fname, mt, mt+attr(mt, "match.length")-1))
  vcf <- readVcf(fname, genome=genomeversion, param=vcfPar)
  stopifnot(all(seqnames(vcf) == chr)) ## QC

  ## set SeqInfo information from Hsapiens
  seqinfo(vcf, new2old=match(seqlevels(si), seqlevels(vcf))) <- si

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

  ## fetch SNVs GRCh38 coordinates from hs37d5 lifted positions by NHLBI ESP
  rr <- info(vcfsnvs)$GRCh38_POSITION
  stopifnot(all(elementNROWS(rr) == 1)) ## QC
  rr <- unlist(rr, use.names=FALSE)

  ## discard variants whose hs37d5 coordinates could not be lifted to GRCh38
  maskNoLift <- rr == "-1"
  vcfsnvs <- vcfsnvs[!maskNoLift, ]
  rr <- rr[!maskNoLift]
  rr <- VariantFiltering:::.str2gr(rr)

  ## discard variants whose hs37d5 coordinates have been lifted to
  ## alternate loci in GRCh38
  maskNoLift <- !seqnames(rr) %in% seqlevels(si)
  vcfsnvs <- vcfsnvs[!maskNoLift, ]
  rr <- rr[!maskNoLift]

  ## set SeqInfo information from Hsapiens
  seqinfo(rr, new2old=match(seqlevels(si), seqlevels(rr))) <- si

  ## fetch minor allele frequency data
  ## according to info(scanVcfHeader(vcf)) the MAF column contains the
  ## "Minor Allele Frequency in percent in the order of EA,AA,All"
  mafValues <- matrix(as.numeric(unlist(info(vcfsnvs)$MAF)), byrow=TRUE, ncol=3) / 100
  colnames(mafValues) <- c("EA_AF", "AA_AF", "AF")

  for (afCol in colnames(mafValues)) {
    message(sprintf("Processing %s SNVs allele frequencies from chromosome %s", afCol, chr))
    mafValuesCol <- mafValues[, afCol]

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
                            provider="NHLBIESP",
                            provider_version="ESP6500SI-V2",
                            citation=citationdata,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn,
                            max_abs_error=max.abs.error)
      saveRDS(obj, file=file.path(pkgname, sprintf("%s.%s.%s.rds", pkgname, afCol, chr)))
    } else {
      warning(sprintf("No MAF values for SNVs in chromosome %s", chr))
    }
  }

  ##
  ## nonSNVs
  ##

  ## fetch nonSNVs GRCh38 coordinates from hs37d5 lifted positions by NHLBI ESP
  rr <- info(vcfnonsnvs)$GRCh38_POSITION
  stopifnot(all(elementNROWS(rr) == 1)) ## QC
  rr <- unlist(rr, use.names=FALSE)

  ## discard variants whose hs37d5 coordinates could not be lifted to GRCh38
  maskNoLift <- rr == "-1"
  vcfnonsnvs <- vcfnonsnvs[!maskNoLift, ]
  rr <- rr[!maskNoLift]
  rr <- VariantFiltering:::.str2gr(rr)

  ## discard variants whose hs37d5 coordinates have been lifted to
  ## alternate loci in GRCh38
  maskNoLift <- !seqnames(rr) %in% seqlevels(si)
  vcfnonsnvs <- vcfnonsnvs[!maskNoLift, ]
  rr <- rr[!maskNoLift]

  ## set SeqInfo information from Hsapiens
  seqinfo(rr, new2old=match(seqlevels(si), seqlevels(rr))) <- si

  if (all(width(rr) == 1))             ## because all lifted positions are width one
    rr <- GRanges(seqinfo=seqinfo(rr)) ## we're going to discard all nonSNVs, aparently
                                       ## they did not lift the whole range per nonSNV

  saveRDS(rr, file=file.path(pkgname, sprintf("%s.GRnonsnv.%s.rds", pkgname, chr)))

  ## fetch minor allele frequency data
  ## according to info(scanVcfHeader(vcf)) the MAF column contains the
  ## "Minor Allele Frequency in percent in the order of EA,AA,All"
  mafValues <- matrix(as.numeric(unlist(info(vcfnonsnvs)$MAF)), byrow=TRUE, ncol=3) / 100
  colnames(mafValues) <- c("EA_AF", "AA_AF", "AF")

  rm(vcf)
  rm(vcfsnvs)
  rm(vcfnonsnvs)
  gc()

  for (afCol in colnames(mafValues)) {
    message(sprintf("Processing %s nonSNVs allele frequencies from chromosome %s", afCol, chr))
    mafValuesCol <- mafValues[, afCol]

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
                            provider="NHLBIESP",
                            provider_version="ESP6500SI-V2",
                            citation=citationdata,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn,
                            max_abs_error=max.abs.error)
      saveRDS(obj, file=file.path(pkgname, sprintf("%s.RLEnonsnv.%s.%s.rds", pkgname, afCol, chr)))
    } else {
      warning(sprintf("No MAF values for nonSNVs in chromosome %s", chr))
    }
  }
}

## save rsIDs assignments from NHLBI ESP
## iterating through every chromosome VCF file
vcfPar <- ScanVcfParam(geno=NA,
                       fixed="ALT",
                       info="GRCh38_POSITION")

message("Starting to process variant identifiers")

rsIDs <- character(0)  ## to store rsIDs annotated by the 1000 genomes project
rsIDgp <- GPos()       ## to store positions of rsIDs
maskSNVs <- logical(0) ## to store a mask whether the variant is an SNV or not

nTotalVar <- 0
for (chr in seqlevels(si)) {
  fname <- sprintf("ESP6500SI-V2-SSA137.GRCh38-liftover.chr%s.snps_indels.vcf", chr)
  vcf <- readVcf(fname, genome=genomeversion, param=vcfPar)
  stopifnot(all(seqnames(vcf) == chr)) ## QC
  nVar <- nrow(vcf)
  nTotalVar <- nTotalVar + nVar
  ## fetch original hs37d5 positions and their rs identifier assignments
  origrr <- rowRanges(vcf)
  ## fetch nonSNVs GRCh38 coordinates from hs37d5 lifted positions by NHLBI ESP
  rr <- info(vcf)$GRCh38_POSITION
  stopifnot(all(elementNROWS(rr) == 1)) ## QC
  rr <- unlist(rr, use.names=FALSE)

  ## discard variants whose hs37d5 coordinates could not be lifted to GRCh38
  maskNoLift <- rr == "-1"
  vcf <- vcf[!maskNoLift, ]
  origrr <- origrr[!maskNoLift]
  rr <- rr[!maskNoLift]
  rr <- VariantFiltering:::.str2gr(rr)

  ## discard variants whose hs37d5 coordinates have been lifted to
  ## alternate loci in GRCh38
  maskNoLift <- !seqnames(rr) %in% seqlevels(si)
  vcf <- vcf[!maskNoLift, ]
  rr <- rr[!maskNoLift]
  origrr <- origrr[!maskNoLift]

  ## set SeqInfo information from Hsapiens
  seqinfo(origrr, new2old=match(seqlevels(si), seqlevels(origrr))) <- si
  seqinfo(rr, new2old=match(seqlevels(si), seqlevels(rr))) <- si

  whrsIDs <- grep("^rs", names(origrr))
  evcf <- expand(vcf)
  maskSNVs <- c(maskSNVs, sapply(relist(isSNV(evcf), alt(vcf)), all)[whrsIDs])
  rm(evcf)
  gc()
  rrTmp <- rr[whrsIDs]
  origrrTmp <- origrr[whrsIDs]
  rrTmp <- resize(rrTmp, width=1, fix="start", ignore.strand=TRUE)
  origrrTmp <- resize(origrrTmp, width=1, fix="start", ignore.strand=TRUE)
  idTmp <- names(origrrTmp)
  names(origrrTmp) <- names(rrTmp) <- NULL
  gpTmp <- as(rrTmp, "GPos")

  stopifnot(all(identical(grep("^rs[0-9]+$", idTmp), 1:length(idTmp)))) ## QC

  rsIDs <- c(rsIDs, idTmp)
  rsIDgp <- c(rsIDgp, gpTmp)
  rm(rrTmp)
  rm(gpTmp)
  rm(idTmp)
  gc()

  message(sprintf("%d variant identifiers processed from chromosome %s", nVar, chr))
}
message(sprintf("%d variants processed in total", nTotalVar))

## save total number of variants
saveRDS(nTotalVar, file=file.path(pkgname, "nov.rds"))

## store mask flagging SNVs
rsIDgp$isSNV <- Rle(maskSNVs)

## double check that all identifiers start with 'rs'
stopifnot(identical(grep("^rs", rsIDs), 1:length(rsIDs))) ## QC

## chop the 'rs' prefix and convert the character ids into integer values
rsIDs <- as.integer(sub(pattern="^rs", replacement="", x=rsIDs))
gc()

## calculate the indices that lead to an increasing values of the (integer) ids
rsIDidx <- order(rsIDs)

## order increasingly the integer values of rs IDs
rsIDs <- rsIDs[rsIDidx]

## save the objects that enable the search for rs ID
saveRDS(rsIDs, file=file.path(pkgname, "rsIDs.rds"))
saveRDS(rsIDidx, file=file.path(pkgname, "rsIDidx.rds"))
saveRDS(rsIDgp, file=file.path(pkgname, "rsIDgp.rds"))
