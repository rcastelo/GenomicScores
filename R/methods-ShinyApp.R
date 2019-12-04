# Calls the Shiny app function
igscores <- function() {
  appDir <- system.file("shinyApp", package="GenomicScores")
  if (appDir == "")
    stop("The GenomicScores shiny app cannot be found within the package.")
  
  shiny::runApp(appDir)
}
