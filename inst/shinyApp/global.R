# Global options, make DT::rendertable print NA values as 'NA'
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

options <- GenomicScores::availableGScores()

AnnotationHub::setAnnotationHubOption("MAX_DOWNLOADS", 600)

##### General functions #######

## imports BED files uploaded by the user through
## the shiny app. it only reads the first three
## columns: chromosome, start(0-based), end(1-based)
## it skips comments and track line
readBed <- function(filename) {

  if (!file.exists(filename))
    stop(sprintf("%s does not exist."))

  con <- file(filename, "r")
  l <- readLines(con, 1, warn=FALSE)
  i <- 0
  while (length(grep("^ *#", l)) > 0 || length(grep("^track", l)) > 0) {
    l <- readLines(con=con, n=1, warn=FALSE)
    i <- i + 1L
  }
  if (length(grep("^ *#", l)) == 0 && length(grep("^track", l)) == 0)
    pushBack(l, con)

  bed <- read.table(con, sep="\t", colClasses=c("character", "integer", "integer"),
                    stringsAsFactors=FALSE)
  close(con)

  GRanges(seqnames=bed[[1]],
          IRanges(bed[[2]]+1L, bed[[3]]))
}

## Returns a dataframe from a GScore object
## ready to be exported as a BED file
createBed <- function(gs) {
  df <- data.frame(
    seqnames=seqnames(gs),
    starts=start(gs)-1,
    ends=end(gs),
    names=c(rep(".", length(gs))),
    scores=elementMetadata(gs)[,1]
  )
  df
}



# Creates a Shiny download button using a GScore object and a 
# string containing the type of document the dwnbtn will generate ('bed' or 'csv')
downloadFile <- function(dt, type){
  downloadHandler(
    filename = function(){
      paste("file", type, sep=".")
    },
    content = function(file){
      df <- createBed(dt)
      switch(type,
             csv = write.table(df, file, quote=F, sep=",", row.names=F, col.names=F),
             bed = write.table(df, file, quote=F, sep="\t", row.names=F, col.names=F))
    }
  )
}



# Validates if a Shiny input is in fact an integer number
not_empty_or_char <- function(input){
  if(is.na(suppressWarnings(as.numeric(input)))){
    "ERROR: Genomic ranges must be numeric and cannot be null"
  } else { NULL }
}



# Validates if Shiny inputs for genomic ranges are integer numbers
is_smaller <- function(rStart, rEnd){
  if(as.numeric(rEnd) < as.numeric(rStart)){
    "ERROR: Ending position must be bigger or equal than Starting position"
  } else { NULL }
}



# Validates if a GRange object created with Shiny inputs (name, rStart and rEnd)
# is in range within the annotation package (phast)
is_within_range <- function(granges, phast){
  annot.pkg <- GRanges(seqnames=seqnames(phast), 
                    IRanges(rep(1,length(seqnames(phast))), seqlengths(seqinfo(phast))))
  seqlevelsStyle(granges) <- seqlevelsStyle(phast)[1]
  if(!identical(granges,subsetByOverlaps(granges, annot.pkg, type="within"))){
    "ERROR: The query genomic ranges are outside the boundaries of the genomic scores object"
  } else { NULL }
}


# Checks with a name string if a package is loaded or attached
.isPkgLoaded <- function(name) {
  paste("package:", name, sep="") %in% search()
}


# Loads a library into the package namespace, mainly used to load annotation
# packages previously installed in user's machine. Uses as parameters the
# package name and its type (GScore is hardcoded in server.R)
.loadAnnotationPackageObject <- function(pkgName) {
  
  annotObj <- NULL
  
  if(options[row.names(options)==pkgName,"Installed"]){
    if (!.isPkgLoaded(pkgName)) {
      suppressPackageStartupMessages(require(pkgName, character.only=TRUE))
    }
    annotObj <- get(pkgName)
  }
  
  if(options[row.names(options)==pkgName,"Cached"]){
    annotObj <- getGScores(pkgName)
  }
  annotObj
}

.installAnnotPkg <-  function(pkgName){
  
  if(options[row.names(options)==pkgName,"BiocManagerInstall"]){
    BiocManager::install(pkgName, update=FALSE)
  } else {
    getGScores(pkgName)
  }
  
  options <<- GenomicScores::availableGScores()
}
