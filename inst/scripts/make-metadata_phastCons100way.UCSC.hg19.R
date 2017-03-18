.GenomicScoresMetadataFromUrl <- function(baseUrl, track) {
  require(GenomeInfoDb)

  ## genus+species to NCBI taxon identifier from GenomeInfoDb (NOT ANYMORE! 17/3/17)
  ## load(system.file("data", "speciesMap.rda", package="GenomeInfoDb"))
  ## speciesMap <- get("speciesMap") ## just to avoid a NOTE during R CMD check

  pcmetadf <- data.frame(title=character(0), species=character(0),
                         taxonomyId=integer(0), genome=character(0),
                         sourceUrl=character(0), sourceVersion=character(0),
                         description=character(0), rDataPath=character(0))

  gd <- readRDS(gzcon(url(sprintf("%s%s/refgenomeGD.rds", baseUrl, track),
                          open="rb")))
  pcRDSfiles <- AnnotationForge:::.getSubDirs(sprintf("%s%s", baseUrl, track))
  pcRDSfiles <- pcRDSfiles[grep(track, pcRDSfiles)]
  ## taxonId <- as.integer(subset(speciesMap, species == organism(gd))$taxon)
  taxonId <- 9606L
  rDataPath <- sprintf("%s/%s", track, pcRDSfiles)
  chr <- sub(".rds", "", sub(paste0(track, "."), "", pcRDSfiles))
  stopifnot(all(chr %in% seqnames(gd))) ## QC
  metadf <- data.frame(title=pcRDSfiles,
                       species=rep(organism(gd), length(pcRDSfiles)),
                       taxonomyId=rep(taxonId, length(pcRDSfiles)),
                       genome=rep(providerVersion(gd), length(pcRDSfiles)),
                       sourceUrl=sprintf("%s%s/%s", baseUrl, track, pcRDSfiles),
                       sourceVersion=rep("3.5.0", length(pcRDSfiles)),
                       description=sprintf("phastCons scores for %s on %s",
                                           organism(gd), chr),
                       rDataPath=rDataPath,
                       stringsAsFactors=FALSE)
  pcmetadf <- rbind(pcmetadf, metadf)
  rownames(pcmetadf) <- NULL
  pcmetadf
}

makeMetadata_phastCons100way.UCSC.hg19 <- function()
{
  baseUrl="http://functionalgenomics.upf.edu/annotationhub/phastCons/"
  meta <- .GenomicScoresMetadataFromUrl(baseUrl, "phastCons100way.UCSC.hg19")
  n <- nrow(meta)
  data.frame(
    BiocVersion=rep("3.5", n),
    Description=meta$description,
    Genome=meta$genome,
    SourceUrl=meta$sourceUrl,
    SourceType=rep("BigWig", n),
    SourceVersion=meta$sourceVersion,
    Species=meta$species,
    TaxonomyId=meta$taxonomyId,
    Title=meta$title,
    ResourceName=meta$title,
    RDataPath=meta$rDataPath,
    Coordinate_1_based=rep(TRUE, n),
    DataProvider=rep("UCSC", n),
    Maintainer=rep("Robert Castelo <robert.castelo@upf.edu>", n),
    RDataClass=rep("Rle", n),
    DispatchClass=rep("RDS", n),
    Location_Prefix=baseUrl,
    Tags=rep(paste("phastCons", "GScores", sep=","), n))
}
metadata <- makeMetadata_phastCons100way.UCSC.hg19()
write.csv(metadata, file="../extdata/metadata_phastCons100way.UCSC.hg19.csv", row.names=FALSE)
