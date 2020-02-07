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
    if(input$annotPackage=="")return()
    .loadAnnotationPackageObject(input$annotPackage, "GScores")
  })
  
  ##### Uploaded Bed file #####
  uploadedBed<- reactive({
    req(input$annotPackage, input$upload)
    readBed(input$upload$datapath)
  })
  
  #### GRange from web with gscores ####
  gsObject <- reactive({
    req(input$annotPackage, input$indOrRange)
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
    
  })
  
  #### GRange from bed file with gscores #####
  gsUpload <- reactive({
    gscores(annotPackage(), uploadedBed(), pop = input$populations)
  })

  
  ################## UI HIDE AND SHOW  ##################
  
  observe({
    if(input$webOrBed == 'web'){
      shinyjs::showElement("webOptions")
      shinyjs::hideElement("upload")
      shinyjs::showElement("printGsWeb")
      shinyjs::hideElement("printGsBed")
      shinyjs::showElement("dwn_web_bed")
      shinyjs::showElement("dwn_web_csv")
      shinyjs::hideElement("dwn_bed_bed")
      shinyjs::hideElement("dwn_bed_csv")
    } else {
      shinyjs::hideElement("webOptions")
      shinyjs::showElement("upload")
      shinyjs::hideElement("printGsWeb")
      shinyjs::showElement("printGsBed")
      shinyjs::hideElement("dwn_web_bed")
      shinyjs::hideElement("dwn_web_csv")
      shinyjs::showElement("dwn_bed_bed")
      shinyjs::showElement("dwn_bed_csv")
    }
  })
  
  ################## GENERATE INPUTS  ######################## 
  
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
    req(input$annotPackage)
    selectInput("populations", "Select an available population", multiple = TRUE,
                choices = populations(annotPackage()))
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
  
  ### Web Datatable
  output$printGsWeb <- DT::renderDataTable({
    if(input$webOrBed != "web" || input$annotPackage == "") return ()
    data.table::as.data.table(gsObject())
  })
  
  ### Bed Datatable
  output$printGsBed <- DT::renderDataTable({
    req(input$upload)
    if(input$webOrBed != "bed" ) return ()
    data.table::as.data.table(gsUpload())
  })
  
  ### Session Info
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })

  
  ################## DOWNLOAD BUTTONS #########################

  output$dwn_web_bed <- downloadFile(gsObject(), "bed")
  output$dwn_web_csv <- downloadFile(gsObject(), "csv")
  output$dwn_bed_bed <- downloadFile(gsUpload(), "bed")
  output$dwn_bed_csv <- downloadFile(gsUpload(), "csv")

}
