# Calls the Shiny app function
igscores <- function() {

  shinydeps <- c("shiny", "shinyjs", "shinydashboard", "magrittr",
                 "shinycustomloader", "data.table", "DT")
  maskshinydeps <- shinydeps %in% installed.packages()
  if (any(!maskshinydeps))
    stop(sprintf("Please install the following packages to use the GenomiScores WebApp:\n\n  %s\n",
                 paste(shinydeps[!maskshinydeps], collapse=", ")))

  appDir <- system.file("shinyApp", package="GenomicScores")
  if (appDir == "")
    stop("The GenomicScores shiny app cannot be found within the package.")
  
  shiny::runApp(appDir)
}
