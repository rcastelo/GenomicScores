server <- function(input, output, session) {
  
  # store in a session Object the available Apkgs that user has
  session$userData$apkgs <- availableGScores()
 
  ################## INPUT VALUES ########################
  
  chromo <- reactive({
    req(input$chromo, annotPackage())
    chromo <- input$chromo
    seqlevelsStyle(chromo) <- seqlevelsStyle(annotPackage())[1]
    validate(
      need(chromo %in% seqnames(annotPackage()), 
           "ERROR: Chromosome name is not present in the GScores object")
    )
    chromo
  })
  
  rStart <- reactive({
    req(input$rStart)
    validate( not_empty_or_char(input$rStart) )
    input$rStart
    })
  
  rEnd <- reactive({
    req(input$rEnd)
    validate( not_empty_or_char(input$rEnd) )
    input$rEnd
  })
  
  fullRange <- reactive(
    paste0(chromo(),":", rStart(),"-", rEnd())
  )
  
  ################# REACTIVE CORE VALUES #######################
  
  ###### Annotation package object #####
  annotPackage <- reactive({
    req(input$annotPackage)
    if(input$annotPackage=="" || !input$annotPackage %in% installed.packages())
      return()
    .loadAnnotationPackageObject(input$annotPackage, "GScores")
  })

  
  ##### Uploaded Bed file #####
  uploadedBed<- reactive({
    req(input$annotPackage, input$upload)
    readBed(input$upload$datapath)
  })
  
  #### GRange from web with gscores ####
  gsObject <- eventReactive(input$run, {
    req(input$annotPackage)
    switch(input$webOrBed,
           web ={
             req(input$indOrRange)
             validate(is_smaller(rStart(), rEnd()))
             validate(is_within_range(chromo(), rStart(), rEnd(), annotPackage()))
             
             switch(input$indOrRange,
                    individual = gscores(annotPackage(), 
                                         GRanges(seqnames=chromo(),
                                                 IRanges(rStart():rEnd(), width=1)),
                                         pop = input$populations),
                    range = gscores(annotPackage(), 
                                    GRanges(fullRange()),
                                    pop = input$populations)
             )
           },
           bed = {
             req(uploadedBed())
             tryCatch({
               gscores(annotPackage(), uploadedBed(), pop = input$populations)
             }, error=function(err){
               stop(print(paste("There seems to be a problem with the bed file", err, sep = "\n")))
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
  
  ##### deactivates 'run' btn until there is a loaded Annotation package
  observe({
    shinyjs::toggleState(id="run",
                         condition = input$annotPackage!="" &&
                           class(annotPackage())=="GScores")
  })
  
  ################## GENERATE INPUTS  ######################## 
  
  
  output$org.cat <- renderUI({
    options <- session$userData$apkgs
    a <- list()
    a[["org"]] <- selectInput("organism", "Select an Organism",
                              selected = "All",
                              choices = c("All" = "All", 
                                          unique(options[options$Installed,]$Organism)))
    a[["cat"]] <-  selectInput("category", "Select a Category",
                               selected = "All",
                               choices = c("All" = "All", 
                                       unique(options[options$Installed,]$Category)))
    tagList(a)
    
  })
  

  output$apkg <- renderUI({
    req(input$organism, input$category)
    options <- session$userData$apkgs
    organism <- input$organism
    category <- input$category
    options <- if(organism=="All") options else na.omit(options[options$Organism==organism,])
    options <- if(category=="All") options else na.omit(options[options$Category==category,])

    tags$div(id="cssref", 
        selectInput("annotPackage", "Select an Annotation Package",
                choices = c("Choose a package" = "", options$Name)))
  })
  
  # this section generates the necessary css style in order to
  # programmatically change the selectInput() choices' colors
  
  output$css.apkgs <- renderUI({
    options <- session$userData$apkgs
    names <- session$userData$apkgs$Name
    tags$style(
      HTML(unlist(
        lapply(names, function(x){
          if(options[options$Name==x,]$Installed) {
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
    options <- session$userData$apkgs
    if(!name=="" && !name %in% options[options$Installed,]$Name){
      showModal(
        modalDialog(
          title = "Installed Annotation Packages",
          "The annotation package you chose is not installed and cannot be used yet. 
          Would you like to install it?",
          div(p("\n")),
          verbatimTextOutput("install.text"),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("install.apkgs", "OK")
          )
        )
      )
    }
  })
  
  observeEvent(input$install.apkgs, {
    if(input$annotPackage %in% installed.packages()) {
      removeModal()
      session$reload()
    } else {
      withCallingHandlers({
        shinyjs::html(id = "install.text", html = paste("Installation of ", input$annotPackage, " in progress.\n"))
        shinyjs::html(id = "install.text", html = "Please wait, this may take a while.\n", add = TRUE)
        shinyjs::html(id = "install.text", html = "You can check the progress on your R console.\n\n", add = TRUE)
        BiocManager::install(input$annotPackage, update=FALSE)
        shinyjs::html(id = "install.text", html = "\n", add = TRUE)
        shinyjs::html(id = "install.text", html = "Installation Finished", add = TRUE)
        },
        message = function(m) {
          shinyjs::html(id = "install.text", html = m$message, add = TRUE)},
        warning = function(m) {
          shinyjs::html(id = "install.text", html = m$message, add = TRUE)},
        error = function(m) {
          shinyjs::html(id = "install.text", html = m$message, add = TRUE)})
      }
    })
  

  
  #### render web parameters (chr name, start, end, choose between range or indv.)
  output$webOptions <- renderUI({
    req(input$webOrBed=="web")
    tagList(
      fluidRow(id="granges.inputs",
               column(4,
                      textInput("chromo", "Chr name", value = "chr22")
               ),
               column(4,
                      textInput("rStart", "Start", value = "50967020")
               ),
               column(4,
                      textInput("rEnd", "End", value = "50967025")
               )
      ),
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
    req(gsObject())
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
