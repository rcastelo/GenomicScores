.GenomicScoresMetadataFromUrl <- function(baseUrl, track) {
  require(GenomeInfoDb)

  ## genus+species to NCBI taxon identifier from GenomeInfoDb
  load(system.file("extdata", "assembly_accessions.rda", package="GenomeInfoDb"))
  assembly_accessions <- get("assembly_accessions") ## just to avoid a NOTE during R CMD check

  gd <- readRDS(gzcon(url(sprintf("%s%s/refgenomeGD.rds", baseUrl, track),
                          open="rb")))
  rdsFiles <- GenomicScores:::.getSubDirs(sprintf("%s%s", baseUrl, track))
  rdsFiles <- rdsFiles[grep(track, rdsFiles)]
  taxonId <- as.integer(subset(assembly_accessions, organism_name == organism(gd))$taxid)[1]
  obj <- readRDS(gzcon(url(sprintf("%s%s/%s", baseUrl, track, rdsFiles[1]),
                           open="rb")))
  rDataPath <- sprintf("%s/%s", track, rdsFiles)
  chr <- sub(".rds", "", sub(paste0(track, "."), "", rdsFiles))
  stopifnot(all(chr %in% seqnames(gd))) ## QC
  metadf <- data.frame(title=rdsFiles,
                       species=rep(organism(gd), length(rdsFiles)),
                       taxonomyId=rep(taxonId, length(rdsFiles)),
                       genome=rep(providerVersion(gd), length(rdsFiles)),
                       sourceUrl=sprintf("%s%s/%s", baseUrl, track, rdsFiles),
                       sourceVersion=rep(metadata(obj)$provider_version, length(rdsFiles)),
                       description=sprintf("phyloP scores for %s on %s", organism(gd), chr),
                       rDataPath=rDataPath,
                       stringsAsFactors=FALSE)
  rownames(metadf) <- NULL
  metadf
}

makeMetadata_phyloP35way.UCSC.mm39 <- function()
{
  biocver <- "3.16"
  baseUrl <- "https://functionalgenomics.upf.edu/annotationhub/phyloP/"
  meta <- .GenomicScoresMetadataFromUrl(baseUrl, "phyloP35way.UCSC.mm39")
  n <- nrow(meta)
  data.frame(
    BiocVersion=rep(biocver, n),
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
    Tags=rep(paste(c("phyloP", "UCSC", "GScores"), collapse=":"), n))
}
metadata <- makeMetadata_phyloP35way.UCSC.mm39()
write.csv(metadata, file="../extdata/metadata_phyloP35way.UCSC.mm39.csv", row.names=FALSE)
