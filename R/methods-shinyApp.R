# Calls the Shiny app function
igscores <- function() {

  shinydeps <- c("shiny", "shinythemes", "shinyjs",
                 "shinycustomloader", "data.table", "DT")
  maskshinydeps <- shinydeps %in% installed.packages()
  if (any(!maskshinydeps))
    stop(sprintf("Please install the following packages to use the GenomiScores WebApp:\n\n  %s\n",
                 paste(shinydeps[!maskshinydeps], collapse=", ")))

  ## import everything from shiny except 'renderDataTable' and
  ## 'dataTableOutput', which clash with DT
  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("shiny"),
                      except=c("renderDataTable", "dataTableOutput"))

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("shinythemes"),
                      vars="shinytheme")

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("shinyjs"),
                      vars=c("useShinyjs", "showElement", "hideElement"))

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("shinycustomloader"),
                      vars="withLoader")

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("data.table"),
                      vars="as.data.table")

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("DT"),
                      vars=c("renderDataTable", "dataTableOutput"))

  appDir <- system.file("shinyApp", package="GenomicScores")
  if (appDir == "")
    stop("The GenomicScores shiny app cannot be found within the package.")
  
  runWebApp <- get("runApp", mode="function")
  runWebApp(appDir)
}
