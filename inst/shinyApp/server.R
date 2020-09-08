server <- function(input, output, session) {

  
  ################# REACTIVE CORE VALUES #######################
  
  ###### Annotation package object #####
  annotPackage <- reactive({
    req(input$annotPackage)
    name <- input$annotPackage
    if(name=="") return()
    .loadAnnotationPackageObject(input$annotPackage)
  })
  
  ##### Uploaded Bed file #####
  uploadedBed<- reactive({
    req(input$annotPackage, input$upload)
    readBed(input$upload$datapath)
  })
  
  #### GRanges object from the selected annotPkg with added GScores ####
  gsObject <- eventReactive(input$run, {
    req(input$annotPackage)
    granges <- tryCatch({
      GRanges(input$granges)
    }, error = function(err){
      stop(print(paste("There is an error with the GRange object:\n", err)))
      return(NULL)
    })
    annot.pkg <- annotPackage()
    if(is.null(input$populations)) {
      population <- defaultPopulation(annot.pkg)
    } else {
      population <- input$populations
      }
    
    switch(input$webOrBed,
           web ={
             req(input$indOrRange)
             validate(is_smaller(start(granges), end(granges)))
             validate(is_within_range(granges, annot.pkg))
             
             tryCatch({
               switch(input$indOrRange,
                      individual = gscores(annot.pkg, 
                                           GRanges(seqnames=seqnames(granges), 
                                                   IRanges(start(granges):end(granges), width=1)),
                                           pop = population),
                      range = gscores(annot.pkg, granges, pop = population)
               )
             }, error = function(err) {
               return(NULL)
             })
           },
           bed = {
             req(uploadedBed())
             tryCatch({
               gscores(annot.pkg, uploadedBed(), pop = population)
             }, error=function(err){
               stop(print(paste("There seems to be a problem with the bed file\n", err)))
               return(NULL)
               }
             )
           })
    
    
  })
  
  
  ################## UI HIDE AND SHOW  ##################
  
  observe({
    if(input$webOrBed == 'web'){
      shinyjs::showElement("webOptions")
      shinyjs::hideElement("upload")
    } else {
      shinyjs::hideElement("webOptions")
      shinyjs::showElement("upload")
    }
  })
  
  observe({
    if(isTruthy(input$granges)){
      shinyjs::enable("run")
    } else {
      shinyjs::disable("run")
    }
  })


  ################## GENERATE INPUTS  ######################## 
  

  output$apkg <- renderUI({
    organism <- input$organism
    category <- input$category
    options <- if(organism=="All") options else options[which(options$Organism==organism),]
    options <- if(category=="All") options else options[which(options$Category==category),]
    tags$div(id="cssref", 
        selectInput("annotPackage", "Select an Annotation Package",
                choices = c("Choose a package" = "", row.names(options))))
  })
  
  # this section generates the necessary css style in order to
  # programmatically change the selectInput() choices' colors
  
  output$css.apkgs <- renderUI({
    req(options)
    options <- availableGScores()
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
            modalButton("Quit")
          )
        )
      )
    }
  })
  
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
    shinyjs::show("refresh")
    })
  
  observeEvent(input$refresh, {
    removeModal()
    shinyjs::runjs("{history.go(0)}")
  })
  

  
  #### render web parameters and choose between range or individual
  output$webOptions <- renderUI({
    req(input$webOrBed=="web")
    tagList(
      textInput("granges", "Genomic Range", value = "chr22:50967020-50967025"),
      radioButtons("indOrRange", "Output type",
                   choices = list("Genomic range" = "range", "Individual positions" = "individual"))
    )
  })
  
  ### render selectinput with population options
  output$pop <- renderUI({
    req(annotPackage())
    selectInput("populations", "Select an available population", multiple = TRUE,
                choices = populations(annotPackage()))
  })
  
  ### render download buttons
  output$down_btn <- renderUI({
    req(gsObject)
    fluidRow(
             downloadButton("dwn_bed", "Download BED"),
             downloadButton("dwn_csv", "Download CSV")
             )
  })
  
  
  
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
