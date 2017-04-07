.GenomicScoresMetadataFromUrl <- function(baseUrl, track) {
  require(GenomeInfoDb)

  ## genus+species to NCBI taxon identifier from GenomeInfoDb (NOT ANYMORE 17/3/17!)
  ## load(system.file("data", "speciesMap.rda", package="GenomeInfoDb"))
  ## speciesMap <- get("speciesMap") ## just to avoid a NOTE during R CMD check

  metadf <- data.frame(title=character(0), species=character(0),
                       taxonomyId=integer(0), genome=character(0),
                       sourceUrl=character(0), sourceVersion=character(0),
                       description=character(0), rDataPath=character(0))

  gd <- readRDS(gzcon(url(sprintf("%s%s/refgenomeGD.rds", baseUrl, track),
                          open="rb")))
  scRDSfiles <- GenomicScores:::.getSubDirs(sprintf("%s%s", baseUrl, track))
  scRDSfiles <- scRDSfiles[grep(track, scRDSfiles)]
  ## taxonId <- as.integer(subset(speciesMap, species == organism(gd))$taxon)
  taxonId <- 9606L
  rDataPath <- sprintf("%s/%s", track, scRDSfiles)
  chr <- sub(".rds", "", sub(paste0(track, "."), "", scRDSfiles))
  stopifnot(all(chr %in% seqnames(gd))) ## QC
  metadf <- data.frame(title=scRDSfiles,
                       species=rep(organism(gd), length(scRDSfiles)),
                       taxonomyId=rep(taxonId, length(scRDSfiles)),
                       genome=rep(providerVersion(gd), length(scRDSfiles)),
                       sourceUrl=sprintf("%s%s/%s", baseUrl, track, scRDSfiles),
                       sourceVersion=rep("3.5.0", length(scRDSfiles)),
                       description=sprintf("CADD scores v1.3 for %s on %s",
                                           organism(gd), chr),
                       rDataPath=rDataPath,
                       stringsAsFactors=FALSE)
  metadf <- rbind(metadf, metadf)
  rownames(metadf) <- NULL
  metadf
}

makeMetadata_CADD.v1.3.hg19 <- function()
{
  baseUrl <- "http://functionalgenomics.upf.edu/annotationhub/cadd/"
  meta <- .GenomicScoresMetadataFromUrl(baseUrl, "cadd.v1.3.hg19")
  n <- nrow(meta)
  data.frame(
    BiocVersion=rep("3.5", n),
    Description=meta$description,
    Genome=meta$genome,
    SourceUrl=meta$sourceUrl,
    SourceType=rep("VCF", n),
    SourceVersion=meta$sourceVersion,
    Species=meta$species,
    TaxonomyId=meta$taxonomyId,
    Title=meta$title,
    ResourceName=meta$title,
    RDataPath=meta$rDataPath,
    Coordinate_1_based=rep(TRUE, n),
    DataProvider=rep("UWashington", n),
    Maintainer=rep("Robert Castelo <robert.castelo@upf.edu>", n),
    RDataClass=rep("Rle", n),
    DispatchClass=rep("RDS", n),
    Location_Prefix=baseUrl,
    Tags=rep(paste("CADD", "GScores", sep=","), n))
}
metadata <- makeMetadata_CADD.v1.3.hg19()
write.csv(metadata, file="../extdata/metadata_CADD.v1.3.hg19.csv", row.names=FALSE)
