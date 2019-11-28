source("library.R")

##### General functions #######

"Confirms if GenomicScores package is installed. If not, returns () so input$selector
has no options. If it is, then checks which gscore packages does the user have installed
in its machine, from the 'avgs.rds' file in the same package"

avAnnotations <- function(){
  if(!require("GenomicScores")) return()
  avgs <- readRDS(system.file("extdata", "avgs.rds", package="GenomicScores"))
  ip <- installed.packages()
  avgs[avgs %in% ip]
}

"Receives a GRange object with scores and creates a data.frame from it"

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

"General download button function, receives the GRange object with scores
and the type of file the user wants (bed or csv)"
  
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
