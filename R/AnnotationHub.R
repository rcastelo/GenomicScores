
## adapted from AnnotationHubData:::.inparanoidMetadataFromUrl
.phastConsMetadataFromUrl <- function(baseUrl, justRunUnitTest) {
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

  if (justRunUnitTest)
    pcmetadf <- pcmetadf[1:2, ]

  pcmetadf
}

## here should go the transformation from the Rle objects to the GScores object
.RleToGScoresRecipe <- function(ahm) {

  AnnotationHubData::outputFile(ahm)
}


## adapted from https://www.bioconductor.org/packages/3.3/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHubRecipes.html
makePhastConsToAHM <-
  function(justRunUnitTest=FALSE,
           BiocVersion=BiocInstaller::biocVersion(),
           baseUrl="http://functionalgenomics.upf.edu/annotationhub/phastCons",
           ...) {

  ## fetch metadata
  meta <- .phastConsMetadataFromUrl(baseUrl, justRunUnitTest, ...)

  ## create list of AnnotationHubMetadata objects
  Map(AnnotationHubMetadata,
      Description=meta$description,
      Genome=meta$genome,
      SourceUrl=meta$sourceUrl,
      SourceType="RData",
      SourceVersion=meta$sourceVersion,
      Species=meta$species,
      TaxonomyId=meta$taxonomyId,
      Title=meta$title,
      RDataPath=meta$rDataPath,
      MoreArgs=list(
        Coordinate_1_based=TRUE,
        DataProvider="UCSC",
        Maintainer="Robert Castelo <robert.castelo@upf.edu>",
        RDataClass="Rle",
        RDataDateAdded=Sys.time(),
        Recipe="GenomicScores:::.RleToGScoresRecipe",
        DispatchClass="GenomicScores",
        Location_Prefix=baseUrl,
        Tags=c("phastCons", "Annotation")))
}
