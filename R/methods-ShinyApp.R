igscore <- function() {
  appDir <- system.file("shinyApp", package="GenomicScores")
  if (appDir == "")
    stop("The GenomicScores Shiny app cannot be found within the package.")
  
  runApp(appDir)
}