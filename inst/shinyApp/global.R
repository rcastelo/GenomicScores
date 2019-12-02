source("library.R")

# shiny options, make DT::rendertable print NA values as 'NA'
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

##### General functions #######

#' Checks GScore packages installed in user's machine
#' @return installed gscore packages
avAnnotations <- function(){
  avgs <- readRDS(system.file("extdata", "avgs.rds", package="GenomicScores"))
  ip <- installed.packages()
  avgs[avgs %in% ip]
}



#' Creates a dataframe from a GScore object
#' @param gs A GScore object
#' @return a data.frame
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



#' Creates a Shiny download button
#' @param dt a GScore object
#' @param type a string containing the type of document the dwnbtn will generate ('bed' or 'csv')
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



#' Validates if a Shiny input is in fact an integer number
#' @param input A string from a Shiny input
not_empty_or_char <- function(input){
  if(is.na(suppressWarnings(as.numeric(input)))){
    "ERROR: Genomic ranges must be numeric and cannot be null"
  } else { NULL }
}



#' Validates if Shiny inputs for genomic ranges are integer numbers
#' @param rStart A string from a Shiny input, starting range
#' @param rEnd A string from a Shiny input, ending range
is_smaller <- function(rStart, rEnd){
  if(as.numeric(rEnd) < as.numeric(rStart)){
    "ERROR: Ending position must be bigger or equal than Starting position"
  }  else { NULL }
}



#' Validates if a GRange object created with Shiny inputs is in 
#' range within the annotation package
#' @param name A string from a Shiny input, chromosome
#' @param rStart A string from a Shiny input, starting range
#' @param rEnd A string from a Shiny input, ending range
is_within_range <- function(name, rStart, rEnd){
  grgsco <- GRanges(seqnames=seqnames(seqinfo(phast)), 
                    IRanges(rep(1,length(seqnames(phast))), seqlengths(seqinfo(phast))))
  gr <- GRanges(seqnames = name, IRanges(start = as.numeric(rStart), end = as.numeric(rEnd)))
  seqlevelsStyle(gr) <- seqlevelsStyle(phast)
  if(!identical(gr,subsetByOverlaps(gr, grgsco, type="within"))){
    "ERROR: The query genomic ranges are outside the boundaries of the genomic scores object"
  }else { NULL }
}