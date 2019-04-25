## This README file explains the steps to download and freeze in this
## annotation package the gnomAD exome allele frequencies. If you use
## these data please cite the following publication:

## Lek et al. Analysis of protein-coding genetic variation in 60,706 humans.
## Nature, 536:285-291, 2016.
## doi: http://dx.doi.org/10.1038/nature19057

## See http://gnomad.broadinstitute.org/terms for further information
## on using this data for your own research

## The data were downloaded from the FTP server at Ensembl as follows:
##
## wget --recursive --no-parent  --reject "index.html" ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/r2.1

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
                                       ## gnomAD but is the assembly to which
                                       ## GRCh37 coordinates were lifted
library(doParallel)

downloadURL <- "ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/r2.1"
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

vcfFilename <- "gnomad.exomes.r2.1.sites.grch38.chr21_noVEP.vcf.gz"
genomeversion <- "GRCh38"
pkgname <- sprintf("MafDb.gnomADex.r2.1.%s", genomeversion)
dir.create(pkgname)

path2vcfs <- "/projects_fg/GenomicScores/gnomad/GRCh38/r2.1/exomes"

vcfHeader <- scanVcfHeader(file.path(path2vcfs, vcfFilename))

namesstdchr <- standardChromosomes(Hsapiens)
stopifnot(all(seqlengths(vcfHeader)[namesstdchr] == seqlengths(Hsapiens)[namesstdchr])) ## QC

## fill up SeqInfo data
si <- keepStandardChromosomes(seqinfo(vcfHeader))
isCircular(si) <- isCircular(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]
genome(si) <- genome(seqinfo(Hsapiens))[match(seqnames(si), seqnames(seqinfo(Hsapiens)))]

## save the GenomeDescription object
refgenomeGD <- GenomeDescription(organism=organism(Hsapiens),
                                 common_name=commonName(Hsapiens),
                                 provider=provider(Hsapiens),
                                 provider_version=providerVersion(Hsapiens),
                                 release_date=releaseDate(Hsapiens),
                                 release_name=releaseName(Hsapiens),
                                 seqinfo=si)
saveRDS(refgenomeGD, file=file.path(pkgname, "refgenomeGD.rds"))

## read INFO column data
infoCols <- rownames(info(vcfHeader))
AFcols <- infoCols[grep("^AF", infoCols)]
exclude <- grep("raw|POPMAX|popmax", AFcols)
if (length(exclude) > 0)
    AFcols <- AFcols[-exclude] ## remove columns such as those for maximum allele frequency among populations

message("Starting to process variants")

## restrict VCF INFO columns to AF values
vcfPar <- ScanVcfParam(geno=NA,
                       fixed="ALT",
                       info=AFcols)

foreach (chr=setdiff(namesstdchr, "MT")) %dopar% { ## no VCF file for MT

  message(sprintf("Processing chromosome %s", chr))

  ## read the whole VCF file into main memory.
  vcf <- readVcf(file.path(path2vcfs, sprintf("gnomad.exomes.r2.1.sites.grch38.chr%s_noVEP.vcf.gz", chr)),
                 genome=genomeversion, param=vcfPar)

  vcf <- keepStandardChromosomes(vcf)
  seqinfo(vcf) <- si

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
  ## and according to the release notes of gnomAD v2.1
  ## "all multi-allelic sites have been split. This means that
  ## multiple lines now have the same chromosome and position."
  ## in such a case we take the maximum MAF by looking at repeated positions
  rrbypos <- split(rr, start(rr))
  rr <- rr[!duplicated(rr)]

  ## fetch allele frequency data
  afValues <- info(vcfsnvs)
  clsValues <- sapply(afValues, class)

  for (j in seq_along(AFcols)) {
    afCol <- AFcols[j]

    message(sprintf("Processing %s SNVs allele frequencies from chromosome %s", afCol, chr))
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

    ## allele frequencies from gnomAD are calculated from alternative alleles,
    ## so some of them we need to turn them into minor allele frequencies (MAF)
    mask <- !is.na(mafValuesCol) & mafValuesCol > 0.5
    if (any(mask))
      mafValuesCol[mask] <- 1 - mafValuesCol[mask]

    mafValuesCol <- relist(mafValuesCol, rrbypos)
    mafValuesCol <- sapply(mafValuesCol, max) ## in multiallelic variants
                                              ## take the maximum allele frequenc

    q <- .quantizer(mafValuesCol)
    x <- .dequantizer(q)
    f <- cut(x, breaks=c(0, 10^c(seq(floor(min(log10(x[x!=0]), na.rm=TRUE)),
                                     ceiling(max(log10(x[x!=0]), na.rm=TRUE)), by=1))),
             include.lowest=TRUE)
    err <- abs(mafValuesCol-x)
    max.abs.error <- tapply(err, f, mean, na.rm=TRUE)

    ## build an integer-RleList object using the 'coverage()' function
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
                            provider="BroadInstitute",
                            provider_version="r2.1",
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

  ## clean up the ranges
  mcols(rr) <- NULL
  names(rr) <- NULL
  gc()

  ## re-order by chromosomal coordinates to deal with wrongly-ordered nonSNVs over multiple VCF lines, e.g.,
  ##   rs766632658        Y 5037579-5037580      * |         <NA>             CA
  ##   rs766632658        Y         5037579      * |         <NA>              C
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
  ## and according to the release notes of gnomAD v2.1
  ## "all multi-allelic sites have been split. This means that
  ## multiple lines now have the same chromosome and position."
  ## in such a case we take the maximum MAF by looking at repeated positions
  rrbypos <- split(rr, paste(start(rr), end(rr), sep="-"))
  rr <- rr[!duplicated(rr)]
  saveRDS(rr, file=file.path(pkgname, sprintf("%s.GRnonsnv.%s.rds", pkgname, chr)))

  for (j in seq_along(AFcols)) {
    afCol <- AFcols[j]

    message(sprintf("Processing %s nonSNVs allele frequencies from chromosome %s", afCol, chr))
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

    ## allele frequencies from gnomAD are calculated from alternative alleles,
    ## so some of them we need to turn them into minor allele frequencies (MAF)
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
                            provider="BroadInstitute",
                            provider_version="r2.1",
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

## save rsIDs assignments from gnomAD
## iterating through every chromosome VCF file
vcfPar <- ScanVcfParam(geno=NA,
                       fixed="ALT",
                       info=NA)

message("Starting to process variant identifiers")

rsIDs <- character(0)  ## to store rsIDs annotated by the 1000 genomes project
rsIDgp <- GPos()       ## to store positions of rsIDs
maskSNVs <- logical(0) ## to store a mask whether the variant is an SNV or not

nTotalVar <- 0
for (chr in setdiff(namesstdchr, "MT")) { ## no VCF file for MT
  vcf <- readVcf(file.path(path2vcfs, sprintf("gnomad.exomes.r2.1.sites.grch38.chr%s_noVEP.vcf.gz", chr)),
                 genome=genomeversion, param=vcfPar)
  nVar <- nrow(vcf)
  nTotalVar <- nTotalVar + nVar
  rr <- rowRanges(vcf)
  mcols(rr) <- NULL
  gc()
  ord <- order(rr)
  rr <- rr[ord]
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

  ## just in case there are multiple rsID assignments separated by
  ## semicolons like it happens with ExAC, assign just the first one
  idTmp <- strsplit(idTmp, ";")
  idTmp <- sapply(idTmp, "[", 1)

  rsIDs <- c(rsIDs, idTmp)
  rsIDgp <- c(rsIDgp, gpTmp)
  rm(rrTmp)
  rm(gpTmp)
  rm(idTmp)
  gc()

  message(sprintf("%d variant identifiers processed from chromosome %s", nVar, chr))
}
message(sprintf("%d variants processed in total", nTotalVar))

## save the total number of variants
saveRDS(nTotalVar, file=file.path(pkgname, "nsites.rds"))

## store mask flagging SNVs
rsIDgp$isSNV <- Rle(maskSNVs)

## double check that all identifiers start with 'rs' and end with a number
stopifnot(identical(grep("^rs[0-9]+$", rsIDs), 1:length(rsIDs))) ## QC

## chop the 'rs' prefix and convert the character ids into integer values
rsIDs <- as.integer(sub(pattern="^rs", replacement="", x=rsIDs))
gc()

## calculate the indices that lead to an increasing values of the (integer) ids
rsIDidx <- order(rsIDs)

## order increasingly the integer values of rsIDs
rsIDs <- rsIDs[rsIDidx]

## save the objects that enable the search for rsID
saveRDS(rsIDs, file=file.path(pkgname, "rsIDs.rds"))
saveRDS(rsIDidx, file=file.path(pkgname, "rsIDidx.rds"))
saveRDS(rsIDgp, file=file.path(pkgname, "rsIDgp.rds"))
