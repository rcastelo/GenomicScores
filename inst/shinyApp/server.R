server <- function(input, output, session) {
  
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
    if(input$annotPackage=="")return()
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
  
  ################## GENERATE INPUTS  ######################## 
  
  output$apkg <- renderUI({
    options <- availableGScores(installed=TRUE)
    organism <- input$organism
    category <- input$category
    options <- if(organism=="All") options else options[options$Organism==organism,]
    options <- if(category=="All") options else options[options$Category==category,]
    selectInput("annotPackage", "Select a GScores object",
                choices = c("Choose an installed annotation package" = "",
                            options$Name))
  })
  
  #### render web parameters (chr name, start, end, choose between range or indv.)
  output$webOptions <- renderUI({
    req(input$webOrBed=="web")
    tagList(
      fluidRow(id="algo",
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
    req(input$annotPackage)
    annotPackage()
  })
  
  ### Citation
  output$citation <- renderPrint({
    req(input$annotPackage)
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

}
