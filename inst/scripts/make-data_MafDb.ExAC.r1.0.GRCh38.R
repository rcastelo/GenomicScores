## This README file explains the steps to download and freeze in this
## annotation package the ExAC allele frequencies. If you use these data
## please cite the following publication:

## Lek et al. Analysis of protein-coding genetic variation in 60,706 humans.
## Nature, 536:285-291, 2016.
## doi: http://dx.doi.org/10.1038/nature19057

## The data was downloaded from the FTP server of the Broad Institute as
## follows:
##
## wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
## wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz.tbi

## Because we are going to lift GRCh37 coordinates over to GRCh38 coordinates
## we need to download first the corresponding UCSC chain file
## download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
##               "hg19ToHg38.over.chain.gz", method="curl")

## The following R script processes the downloaded data to
## transform the allele frequencies into minor allele frequencies
## and store them using only one significant digit for values < 0.1
## and two significant digits for values > 0.1, to reduce the memory
## footprint through RleList objects

library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(VariantAnnotation)
library(BSgenome.Hsapiens.NCBI.GRCh38) ## this is not the assembly used by
                                       ## ExAC but is the assembly to which
                                       ## GRCh37 coordinates will be lifted
library(rtracklayer)

downloadURL <- "ftp://ftp.broadinstitute.org/pub/ExAC_release/release1"
citationdata <- bibentry(bibtype="Article",
                         author=c(person("Monkol Lek"), person("Konrad J. Karczewski"),
                                  person("Eric V. Minikel"), person("Kaitlin E. Samocha"),
                                  person("Eric Banks"), person("Timothy Fennell"),
                                  person("Anne H. O'Donnell-Luria"), person("James S. Ware"),
                                  person("Andrew J. Hill"), person("Beryl B. Cummings"),
                                  person("Taru Tukiainen"), person("Daniel P. Birnbaum"),
                                  person("Jack A. Kosmicki"), person("Laramie E. Duncan"),
                                  person("Karol Estrada"), person("Fengmei Zhao"),
                                  person("James Zou"), person("Emma Pierce-Hoffman"),
                                  person("Joanne Berghout"), person("David N. Cooper"),
                                  person("Nicole Deflaux"), person("Mark DePristo"),
                                  person("Ron Do"), person("Jason Flannick"),
                                  person("Menachem Fromer"), person("Laura Gauthier"),
                                  person("Jackie Goldstein"), person("Namrata Gupta"),
                                  person("Daniel Howrigan"), person("Adam Kiezun"),
                                  person("Mitja L. Kurki"), person("Ami Levy Moonshine"),
                                  person("Pradeep Natarajan"), person("Lorena Orozco"),
                                  person("Gina M. Peloso"), person("Ryan Poplin"),
                                  person("Manuel A. Rivas"), person("Valentin Ruano-Rubio"),
                                  person("Samuel A. Rose"), person("Douglas M. Ruderfer"),
                                  person("Khalid Shakir"), person("Peter D. Stenson"),
                                  person("Christine Stevens"), person("Brett P. Thomas"),
                                  person("Grace Tiao"), person("Maria T. Tusie-Luna"),
                                  person("Ben Weisburd"), person("Hong-Hee Won"),
                                  person("Dongmei Yu"), person("David M. Altshuler"),
                                  person("Diego Ardissino"), person("Michael Boehnke"),
                                  person("John Danesh"), person("Stacey Donnelly"),
                                  person("Roberto Elosua"), person("Jose C. Florez"),
                                  person("Stacey B. Gabriel"), person("Gad Getz"),
                                  person("Stephen J. Glatt"), person("Christina M. Hultman"),
                                  person("Sekar Kathiresan"), person("Markku Laakso"),
                                  person("Steven McCarroll"), person("Mark L. McCarthy"),
                                  person("Dermot McGovern"), person("Ruth McPherson"),
                                  person("Benjamein M. Neale"), person("Aarno Palotie"),
                                  person("Shaun M. Purcell"), person("Danish Saleheen"),
                                  person("Jeremiah M. Scharf"), person("Pamela Sklar"),
                                  person("Patrick F. Sullivan"), person("Jaakko Tuomilehto"),
                                  person("Ming T. Tsuang"), person("Hugh C. Watkins"),
                                  person("James G. Wilson"), person("Mark J. Daly"),
                                  person("Daniel G. MacArthur"), person("Exome Aggregation Consortium")),
                         title="Analysis of protein-coding genetic variation in 60,706 humans",
                         journal="Nature",
                         volume="536",
                         pages="285-291",
                         year="2016",
                         doi="10.1038/nature19057")

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

vcfFilename <- "ExAC.r1.sites.vep.vcf.gz"
genomeversion <- "GRCh38"
pkgname <- sprintf("MafDb.ExAC.r1.0.%s", genomeversion)
dir.create(pkgname)

vcfHeader <- scanVcfHeader(vcfFilename)

## seqinfo in vcf is from GRCh37
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
ANcols <- c("AN", infoCols[grep("AN_", infoCols)])
exclude <- grep("POPMAX", ANcols)
if (length(exclude) > 0)
    ANcols <- ANcols[-exclude] ## remove colums such as those for maximum allele frequency values among populations
ACcols <- sub("AN", "AC", ANcols)
pops <- substr(ACcols, 4, 20)
AFcols <- gsub("^_", "", paste(pops, "AF", sep="_"))

hg19tohg38chain <- import.chain("hg19ToHg38.over.chain")

message("Starting to process variants")

## restrict VCF INFO columns to AC and AN values
vcfPar <- ScanVcfParam(geno=NA,
                       fixed=c("ALT", "FILTER"),
                       info=c(ACcols, ANcols))

## read the whole VCF file into main memory. the resulting object 'vcf' takes 2.5Gb of RAM
vcf <- readVcf(vcfFilename, genome="hs37d5", param=vcfPar)

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
seqinfo(rr, new2old=match(seqnames(seqinfo(Hsapiens)), seqnames(seqinfo(rr)))) <- seqinfo(Hsapiens)

## store number of SNVs
nsites <- length(rr)

## fetch allele frequency data
acanValues <- info(vcfsnvs)[lomask, ]
clsValues <- sapply(acanValues, class)

for (j in seq_along(ACcols)) {
  acCol <- ACcols[j]
  anCol <- ANcols[j]
  afCol <- AFcols[j]

  message(sprintf("Processing %s SNVs", afCol))
  acValuesCol <- acanValues[[acCol]]
  anValuesCol <- acanValues[[anCol]]
  if (clsValues[acCol] == "numeric" || clsValues[acCol] == "Numeric" ||
      clsValues[acCol] == "integer" || clsValues[acCol] == "Integer") {
    acValuesCol <- as.numeric(acValuesCol)
  } else if (clsValues[acCol] == "character") {
    acValuesCol <- as.numeric(sub(pattern="\\.", replacement="0", acValuesCol))
  } else if (clsValues[acCol] == "CompressedNumericList" || ## in multiallelic variants take the
           clsValues[acCol] == "CompressedIntegerList") {   ## maximum allele count
    acValuesCol <- as.numeric(sapply(acValuesCol, max))
  } else if (clsValues[acCol] == "CompressedCharacterList") {
    acValuesCol <- sapply(NumericList(lapply(acValuesCol, sub, pattern="\\.", replacement="0")), max)
  }

  if (clsValues[anCol] == "character")
    anValuesCol <- as.numeric(sub(pattern="\\.", replacement="0", anValuesCol))

  mafValuesCol <- acValuesCol / anValuesCol

  ## allele frequencies from ExAC are calculated from allele counts from alternative alleles,
  ## so some of them we need to turn them into minor allele frequencies (MAF)
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

  ## build an integer-RleList object using the 'coverage()' function
  rlelst <- coverage(rr, weight=q)
  ## build an integer-Rle object of maskREF using the 'coverage()' function
  maskREFobj <- coverage(rr, weight=maskREF+0L)

  seqs <- names(rlelst)
  idxbychr <- split(seq_len(length(mafValuesCol)), decode(seqnames(rr)))
  max.abs.error <- lapply(idxbychr,
                          function(i, err, f) tapply(err[i], f[i], mean, na.rm=TRUE),
                          err, f)

  ## build ECDF of MAF values
  Fn <- lapply(idxbychr,
               function(i, x) {
                 fn <- NA
                 if (any(!is.na(x[i]))) {
                   if (length(unique(x[i])) <= 10000) {
                     fn <- ecdf(x[i])
                   } else {
                     fn <- ecdf(sample(x[i], size=10000, replace=TRUE))
                   }
                 }
                 fn
               }, mafValuesCol)

  ## coerce to raw-Rle, add metadata and save
  for (k in seq_along(seqs)) {
    if (any(runValue(rlelst[[seqs[k]]]) != 0)) {
      obj <- rlelst[[seqs[k]]]
      objref <- maskREFobj[[seqs[k]]]
      runValue(obj) <- as.raw(runValue(obj))
      runValue(objref) <- as.raw(runValue(objref))
      metadata(obj) <- list(seqname=seqs[k],
                            provider="BroadInstitute",
                            provider_version="r1.0",
                            citation=citationdata,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn[[seqs[k]]],
                            max_abs_error=max.abs.error[[seqs[k]]],
                            maskREF=objref)
      saveRDS(obj, file=file.path(pkgname, sprintf("%s.%s.%s.rds", pkgname, afCol, seqs[k])))
    }
  }
}

rm(vcfsnvs)
gc()

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

## clean up the ranges
mcols(rr) <- NULL
names(rr) <- NULL
gc()

## fill up SeqInfo data
seqinfo(rr, new2old=match(seqnames(seqinfo(Hsapiens)), seqnames(seqinfo(rr)))) <- seqinfo(Hsapiens)

## fetch allele frequency data
acanValues <- info(vcfnonsnvs)[lomask, ]
clsValues <- sapply(acanValues, class)

## re-order by chromosomal coordinates to deal with wrongly-ordered nonSNVs over multiple VCF lines
ord <- order(rr)
rr <- rr[ord]
acanValues <- acanValues[ord, ]

## according to https://samtools.github.io/hts-specs/VCFv4.3.pdf
## "It is permitted to have multiple records with the same POS"
posids <- paste(seqnames(rr), start(rr), end(rr), sep="-")
rrbypos <- split(rr, posids)
rr <- rr[!duplicated(rr)]
## put back the genomic order
posids <- paste(seqnames(rr), start(rr), end(rr), sep="-")
mt <- match(posids, names(rrbypos))
stopifnot(all(!is.na(mt))) ## QC
rrbypos <- rrbypos[mt]

rm(ord)
rm(posids)
rm(mt)
gc()

## store number of SNVs
nsites <- nsites + length(rr)

## split genomic ranges by chromosome
stopifnot(identical(order(seqnames(rr)), 1:length(rr))) ## QC
rrbychr <- split(rr, seqnames(rr))
stopifnot(identical(unlist(rrbychr, use.names=FALSE), rr)) ## QC

seqs <- names(rrbychr)
for (j in seq_along(seqs))
  if (length(rrbychr[[seqs[j]]]) > 0)
    saveRDS(rrbychr[[seqs[j]]],
            file=file.path(pkgname, sprintf("%s.GRnonsnv.%s.rds", pkgname, seqs[j])))

for (j in seq_along(ACcols)) {
  acCol <- ACcols[j]
  anCol <- ANcols[j]
  afCol <- AFcols[j]

  message(sprintf("Processing %s nonSNVs", afCol))
  acValuesCol <- acanValues[[acCol]]
  anValuesCol <- acanValues[[anCol]]
  if (clsValues[acCol] == "numeric" || clsValues[acCol] == "Numeric" ||
      clsValues[acCol] == "integer" || clsValues[acCol] == "Integer") {
    acValuesCol <- as.numeric(acValuesCol)
  } else if (clsValues[acCol] == "character") {
    acValuesCol <- as.numeric(sub(pattern="\\.", replacement="0", acValuesCol))
  } else if (clsValues[acCol] == "CompressedNumericList" || ## in multiallelic variants take the
           clsValues[acCol] == "CompressedIntegerList") {   ## maximum allele count
    acValuesCol <- as.numeric(sapply(acValuesCol, max))
  } else if (clsValues[acCol] == "CompressedCharacterList") {
    acValuesCol <- sapply(NumericList(lapply(acValuesCol, sub, pattern="\\.", replacement="0")), max)
  } else {
    stop(sprintf("Uknown class for holding AF values (%s)", clsValues[acCol]))
  }

  if (clsValues[anCol] == "character")
    anValuesCol <- as.numeric(sub(pattern="\\.", replacement="0", anValuesCol))

  mafValuesCol <- acValuesCol / anValuesCol

  ## allele frequencies from ExAC are calculated from allele counts from alternative alleles,
  ## so some of them we need to turn them into minor allele frequencies (MAF)
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

  idxbychr <- relist(seq_len(length(mafValuesCol)), rrbychr)
  max.abs.error <- lapply(idxbychr,
                          function(i, err, f) tapply(err[i], f[i], mean, na.rm=TRUE),
                          err, f)

  ## build ECDF of MAF values
  Fn <- lapply(idxbychr,
               function(i, x) {
                 fn <- NA
                 if (any(!is.na(x[i]))) {
                   if (length(unique(x[i])) <= 10000) {
                     fn <- ecdf(x[i])
                   } else {
                     fn <- ecdf(sample(x[i], size=10000, replace=TRUE))
                   }
                 }
                 fn
               }, mafValuesCol)
  qbychr <- relist(q, rrbychr)
  stopifnot(identical(sapply(rrbychr, length), sapply(qbychr, length))) ## QC
  refbychr <- relist(maskREF, rrbychr)
  stopifnot(identical(sapply(rrbychr, length), sapply(refbychr, length))) ## QC

  ## coerce to raw-Rle, add metadata and save
  for (k in seq_along(seqs)) {
    if (length(rrbychr[[seqs[k]]]) > 0) {
      obj <- Rle(qbychr[[seqs[k]]])
      objref <- Rle(refbychr[[seqs[k]]])
      runValue(obj) <- as.raw(runValue(obj))
      runValue(objref) <- as.raw(runValue(objref))
      metadata(obj) <- list(seqname=seqs[k],
                            provider="BroadInstitute",
                            provider_version="r1.0",
                            citation=citationdata,
                            download_url=downloadURL,
                            download_date=format(Sys.Date(), "%b %d, %Y"),
                            reference_genome=refgenomeGD,
                            data_pkgname=pkgname,
                            qfun=.quantizer,
                            dqfun=.dequantizer,
                            ecdf=Fn[[seqs[k]]],
                            max_abs_error=max.abs.error[[seqs[k]]],
                            maskREF=objref)
      saveRDS(obj, file=file.path(pkgname, sprintf("%s.RLEnonsnv.%s.%s.rds", pkgname, afCol, seqs[k])))
    }
  }
}

## save the total number of variants
saveRDS(nsites, file=file.path(pkgname, "nsites.rds"))

## save rsID assignments from ExAC
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

## clean up the ranges
mcols(rr) <- NULL
gc()

## fill up SeqInfo data
seqinfo(rr, new2old=match(seqnames(seqinfo(Hsapiens)), seqnames(seqinfo(rr)))) <- seqinfo(Hsapiens)

whrsIDs <- grep("^rs", names(rr))
rsIDrr <- rr[whrsIDs]
mcols(rsIDrr) <- NULL
rsIDrr <- resize(rsIDrr, width=1, fix="start", ignore.strand=TRUE)
rsIDs <- names(rsIDrr)
names(rsIDrr) <- NULL
rsIDgp <- as(rsIDrr, "GPos")
rsIDgp$isSNV <- Rle(maskSNVs[whrsIDs]) ## store mask flagging SNVs

## some rsID assignments consist of several rsIDs separated by semicolon ';',
## assign just the first one
rsIDs <- strsplit(rsIDs, ";")
rsIDs <- sapply(rsIDs, "[", 1)

## double check that all identifiers start with 'rs' and end with a number
stopifnot(identical(grep("^rs[0-9]+$", rsIDs), 1:length(rsIDs))) ## QC

## chop the 'rs' prefix and convert the character ids into integer values
rsIDs <- as.integer(sub(pattern="^rs", replacement="", x=rsIDs))

## calculate the indices that lead to an increasing values of the (integer) ids
rsIDidx <- order(rsIDs)

## order increasingly the integer values of rsIDs
rsIDs <- rsIDs[rsIDidx]

## save the objects that enable the search for rsID
saveRDS(rsIDs, file=file.path(pkgname, "rsIDs.rds"))
saveRDS(rsIDidx, file=file.path(pkgname, "rsIDidx.rds"))
saveRDS(rsIDgp, file=file.path(pkgname, "rsIDgp.rds"))
