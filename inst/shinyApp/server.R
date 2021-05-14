server <- function(input, output, session) {
  
  # This will funciton as a flag for the download buttons:
  # if the gscores object has errors, the buttons will not be printed
  areThereErrors <- reactiveValues("Yes" = FALSE)
  
  ################## UI HIDE AND SHOW  ##################

  # change inputs for 'web' or 'bed' options
  observe({
    if(input$webOrBed == 'web'){
      shinyjs::showElement("webOptions")
      shinyjs::hideElement("upload")
    } else {
      shinyjs::hideElement("webOptions")
      shinyjs::showElement("upload")
    }
  })
  
  # deactivate 'run' button until there is any input in 
  # input$granges and input$annotPackage or 
  # input$upload and input%annotPackage
  observe({
    if((isTruthy(input$granges) && isTruthy(input$annotPackage)) ||
       (isTruthy(input$upload) && isTruthy(input$annotPackage))){
      shinyjs::enable("run")
    } else {
      shinyjs::disable("run")
    }
  })


  ################## GENERATE INPUTS  ######################## 
  
  # Input for annotation pkgs, it's updated with choices in 
  # 'Organism' and 'Category' selectInput
  
  observe({
    organism <- input$organism
    category <- input$category
    
    options <- if(organism=="All") options else options[which(options$Organism==organism),]
    options <- if(category=="All") options else options[which(options$Category==category),]
    
    updateSelectInput(session, "annotPackage",
                      choices = c("Choose a package" = "", row.names(options))
    )
  })
  
  # this section generates the necessary css style in order to
  # programmatically change the annot.pkgs colors
  
  output$css.apkgs <- renderUI({
    req(options)
    names <- row.names(options)
    tags$style(
      HTML(unlist(
        lapply(names, function(x){
          if(options[row.names(options)==x,]$Installed || options[row.names(options)==x,]$Cached) {
            sprintf(
              "#cssref .selectize-dropdown-content > .option[data-value='%s']
              { color: green; font-weight: bold; }", x)
          } else {
            sprintf(
            "#cssref .selectize-dropdown-content > .option[data-value='%s']
            { color: red; font-weight: bold; }", x)
          }
        })
      )
      )
    )
  })

  # Observer for input$annotPackage: if user clicks a red option, a new
  # modal window appears that let's the user install the annot.pkg
  
  observeEvent(input$annotPackage, {
    name <- input$annotPackage
    if(!name=="" && !(options[row.names(options)==name,]$Installed ||
                      options[row.names(options)==name,]$Cached))
      {      
      showModal(
        modalDialog(
          title = "Annotation Packages",
          div(id="install.text", "The annotation package you chose is not installed and cannot be used yet.
          Would you like to install it?"),
          div(p("\n")),
          verbatimTextOutput("install.progress"),
          div(id="install.finished"),
          footer = tagList(
            actionButton("install.apkgs", "OK"),
            shinyjs::hidden(actionButton("refresh", "Refresh")),
            span(id="cancel", modalButton("Cancel")),
            shinyjs::showElement("cancel")
          )
        )
      )
    }
  })
  
  # Observer for the 'Ok' button in modal: this will call the .installAnnotPkg fun
  # that will install the pkg or query it in AnnotHub()
  observeEvent(input$install.apkgs, {
    shinyjs::hideElement("install.apkgs")
    withCallingHandlers({
      shinyjs::html(id = "install.text", html = paste("
        <p>Installation of <b>", input$annotPackage, "</b> in progress.<p>
        <p>This may take a while.</p>
        <p>You can check the progress on your R console.</p>"))
      .installAnnotPkg(input$annotPackage)
      shinyjs::html(id = "install.finished", html = "</br><p>Installation Finished</p>", add = TRUE)
      shinyjs::html(id = "install.finished", html = "<p>Please take note: in order to use a new annotation package, 
                      you must refresh this session</p>", add = TRUE)
    },
    message = function(m) {
      shinyjs::html(id = "install.progress", html = m$message, add = TRUE)},
    warning = function(m) {
      shinyjs::html(id = "install.progress", html = m$message, add = TRUE)},
    error = function(m) {
      shinyjs::html(id = "install.progress", html = m$message, add = TRUE)})
    shinyjs::hideElement("cancel")
    shinyjs::show("refresh")
    })
  
  
  # Observer for the "refresh" button: it will reload the session with a js function
  observeEvent(input$refresh, {
    removeModal()
    shinyjs::runjs("{history.go(0)}")
  })
  
  ### render selectInput with population options
  output$pop <- renderUI({
    req(annotPackage())
    selectInput("populations", "Populations", multiple = TRUE,
                choices = populations(annotPackage()))
  })
  

  #### render web parameters and choose between range or individual
  output$webOptions <- renderUI({
    req(input$webOrBed=="web")
    tagList(
      textInput("granges", "Genomic Range", placeholder = "chromosome:initial-final"),
      radioButtons("indOrRange", "Output type",
                   choices = list("Genomic range" = "range", "Individual positions" = "individual"))
    )
  })
  

  
  ### render download buttons
  output$down_btn <- renderUI({
    req(areThereErrors$Yes==TRUE)
    fluidRow(
      column(
        width = 12,
        align = "center",
        downloadButton("dwn_bed", "BED file"),
        downloadButton("dwn_csv", "CSV file")
        )
      )
  })
  
  ################# REACTIVE CORE VALUES #######################
  
  ###### Annotation package object #####
  annotPackage <- reactive({
    req(input$annotPackage)
    name <- input$annotPackage
    if(name=="") return()
    .loadAnnotationPackageObject(input$annotPackage)
  })  %>% bindCache(input$annotPackage)
  
  ##### Uploaded Bed file #####
  uploadedBed<- reactive({
    req(input$annotPackage, input$upload)
    readBed(input$upload$datapath)
  }) # %>% bindCache(input$annotPackage, input$upload)
  
  ### GRanges Object ###
  granges <- reactive({
    req(annotPackage())
    areThereErrors$Yes <- FALSE
    if(input$webOrBed=="web"){
      req(input$granges)
      validate(need(try(GRanges(input$granges), silent=TRUE),
      "Error: The character vector to convert to a GRanges object must contain strings of
      the form 'chr:start-end' or 'chr:start-end:strand', or 'chr:pos' or 'chr:pos:strand'.
      'start', 'end' and 'pos' must be numeric.
      Note that '..' is a valid alternate start/end separator.
      Strand can be '+', '-', '*', or missing."))
      granges <- GRanges(input$granges)
      validate(is_smaller(start(granges), end(granges)))
      validate(is_within_range(granges, annotPackage()))
      if(input$indOrRange=="individual"){
        granges <- GRanges(seqnames=seqnames(granges), IRanges(start(granges):end(granges),width=1))
      }
    } else {
      req(uploadedBed())
      granges <- uploadedBed()
    }
    areThereErrors$Yes <- TRUE
    return(granges)
  })
  

  
  #### GRanges object from the selected annotPkg with added GScores ####
  gsObject <- reactive({
    req(annotPackage(), granges())
    areThereErrors$Yes <- FALSE      # reset the error flag to FALSE
    annot.pkg <- annotPackage()
    granges <- granges()
    if(is.null(input$populations)) {
      population <- defaultPopulation(annot.pkg)
    } else {
      population <- input$populations
    }
    
   gsObject <- gscores(annot.pkg, granges, pop=population)
    
    areThereErrors$Yes <- TRUE      # error flag to TRUE if everything went ok
    
    return(gsObject)
    
  }) %>%  
    bindCache(annotPackage(), input$populations, granges())  %>%
    bindEvent(input$run)

  ################## OUTPUT TEXT AND TABLES ############################## 

  ### Annotation Package
  output$annotPackageInfo <- renderPrint({
    req(annotPackage())
    annotPackage()
  })
  
  ### Citation
  output$citation <- renderPrint({
    req(annotPackage())
    citation(annotPackage())
  })
  
  ### Datatable
  output$printGs <- DT::renderDataTable({
    req(gsObject())
    if( input$annotPackage == "") return ()
    data.table::as.data.table(gsObject())
  })
  

  ### Session Info
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })

  
  ################## DOWNLOAD BUTTONS #########################

  output$dwn_bed <- downloadFile(gsObject(), "bed")
  output$dwn_csv <- downloadFile(gsObject(), "csv")
  
  
  
  ################## QUIT BUTTON ############################## 
  
  observeEvent(input$quit, {
    stopApp()
  })

}
