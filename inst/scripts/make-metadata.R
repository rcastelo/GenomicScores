.phastConsMetadataFromUrl <- function(baseUrl) {
  baseUrl <- 'http://functionalgenomics.upf.edu/annotationhub/phastCons'
  subDirs <- AnnotationForge:::.getSubDirs(baseUrl)
  subDirs <- subDirs[!subDirs %in% "/annotationhub/"]
  phastConsTracks <- sub("/", "", subDirs)

  ## genus+species to NCBI taxon identifier from GenomeInfoDb
  load(system.file("data", "speciesMap.rda", package="GenomeInfoDb"))
  speciesMap <- get("speciesMap") ## just to avoid a NOTE during R CMD check

  pcmetadf <- data.frame(title=character(0), species=character(0),
                         taxonomyId=integer(0), genome=character(0),
                         sourceUrl=character(0), sourceVersion=character(0),
                         description=character(0), rDataPath=character(0))

  for (pct in phastConsTracks) {
    gd <- readRDS(gzcon(url(sprintf("%s/%s/refgenomeGD.rds", baseUrl, pct),
                            open="rb")))
    pcRDSfiles <- AnnotationForge:::.getSubDirs(sprintf("%s/%s", baseUrl, pct))
    pcRDSfiles <- pcRDSfiles[grep(pct, pcRDSfiles)]
    taxonId <- as.integer(subset(speciesMap, species == organism(gd))$taxon)
    rDataPath <- sprintf("%s/%s", pct, pcRDSfiles)
    chr <- sub(".rds", "", sub(paste0(pct, "."), "", pcRDSfiles))
    stopifnot(all(chr %in% seqnames(gd))) ## QC
    metadf <- data.frame(title=pcRDSfiles,
                         species=rep(organism(gd), length(pcRDSfiles)),
                         taxonomyId=rep(taxonId, length(pcRDSfiles)),
                         genome=rep(providerVersion(gd), length(pcRDSfiles)),
                         sourceUrl=sprintf("%s/%s/%s", baseUrl, pct, pcRDSfiles),
                         sourceVersion=rep("3.4.0", length(pcRDSfiles)),
                         description=sprintf("phastCons scores for %s on %s",
                                             organism(gd), chr),
                         rDataPath=rDataPath,
                         stringsAsFactors=FALSE)
    pcmetadf <- rbind(pcmetadf, metadf)
  }
  rownames(pcmetadf) <- NULL
  pcmetadf
}

makeMetadata <- function() 
{
  baseUrl="http://functionalgenomics.upf.edu/annotationhub/phastCons"
  meta <- .phastConsMetadataFromUrl(baseUrl)
  n <- nrow(meta)
  data.frame(
    BiocVersion=rep("3.4", n),
    Description=meta$description,
    Genome=meta$genome,
    SourceUrl=meta$sourceUrl,
    SourceType=rep("FIXME", n),
    SourceVersion=meta$sourceVersion,
    Species=meta$species,
    TaxonomyId=meta$taxonomyId,
    Title=meta$title,
    RDataPath=meta$rDataPath,
    Coordinate_1_based=rep(TRUE, n),
    DataProvider=rep("UCSC", n),
    Maintainer=rep("Robert Castelo <robert.castelo@upf.edu>", n),
    RDataClass=rep("Rle", n),
    Recipe=rep(NA, n),
    DispatchClass=rep("GenomicScores", n),
    Location_Prefix=meta$sourceUrl,
    Tags=rep(paste("phastCons", "Annotation", sep=","), n))
}
metadata <- makeMetadata()
write.csv(metadata, file="../extdata/metadata.csv", row.names=FALSE)
