# Global options, make DT::rendertable print NA values as 'NA'
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

##### General functions #######

# Returns GScore packages installed in user's machine
avAnnotations <- function(){
  pp <- system.file("scripts", package="GenomicScores")
  mkdatafnames <- list.files(pp, pattern="make-data_*")
  gspkgnames <- sub("make-data_", "", mkdatafnames, fixed=TRUE)
  gspkgnames <- sub(".R", "", gspkgnames, fixed=TRUE)
  ip <- installed.packages()
  gspkgnames[gspkgnames %in% ip]
}



# Returns a dataframe from a GScore object
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
is_within_range <- function(name, rStart, rEnd, phast){
  grgsco <- GRanges(seqnames=seqnames(seqinfo(phast)), 
                    IRanges(rep(1,length(seqnames(phast))), seqlengths(seqinfo(phast))))
  gr <- GRanges(seqnames = name, IRanges(start = as.numeric(rStart), end = as.numeric(rEnd)))
  seqlevelsStyle(gr) <- seqlevelsStyle(phast)
  if(!identical(gr,subsetByOverlaps(gr, grgsco, type="within"))){
    "ERROR: The query genomic ranges are outside the boundaries of the genomic scores object"
  } else { NULL }
}


# Checks with a name string if a package is loaded or attached
.isPkgLoaded <- function(name) {
  (paste("package:", name, sep="") %in% search()) ||
    (name %in% loadedNamespaces())
}


# Loads a library into the package namespace, mainly used to load annotation
# packages previously installed in user's machine. Uses as parameters the
# package name and its type (GScore is hardcoded in server.R)
.loadAnnotationPackageObject <- function(pkgName, pkgType, verbose=TRUE) {
  
  callobj <- match.call()
  annotObj <- NULL
  
  if (is.character(pkgName)) {
    if (!pkgName %in% installed.packages(noCache=TRUE)[, "Package"])
      stop(sprintf("please install the Bioconductor package %s.", pkgName))
    if (!.isPkgLoaded(pkgName)) {
      if (verbose)
        message("Loading ", pkgType, " annotation package ", pkgName)
      if (!suppressPackageStartupMessages(require(pkgName,
                                                  character.only=TRUE)))
        stop(sprintf("package %s could not be loaded.", pkgName))
    }
    tryCatch({
      annotObj <- get(pkgName)
    }, error=function(err) {
      stop(sprintf("The annotation package %s should automatically load
an %s object with the same name as the package.", pkgName, pkgType))
    })
  } else if (class(pkgName) != pkgType)
    stop(sprintf("'%s' is not the name of an '%s' annotation package or
an '%s' annotation object itself.",
                 pkgName, pkgType, pkgType))
  else
    annotObj <- pkgName
  
  if (!is(annotObj, pkgType))
    stop(sprintf("The object loaded with name %s is not an '%s' object.",
                 ifelse(is.character(pkgName), pkgName,
                        gettext(callobj)[2])), pkgType)
  
  annotObj
}
