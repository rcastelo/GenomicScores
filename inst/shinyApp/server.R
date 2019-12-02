server <- function(input, output, session) {
  
  ################## INPUT VALUES ########################
  
  annotPackage <- reactive(input$annotPackage)
  
  chromo <- reactive({
    req(input$chromo, phast())
    validate(
      need(input$chromo %in% seqnames(phast()), 
           "ERROR: Chromosome name is not present in this annotation package")
    )
    input$chromo
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
  phast <- reactive({
    if(annotPackage()=="")return()
    get(annotPackage())
  })
  
  ##### Uploaded Bed file #####
  uploadedBed<- reactive({
    req(input$annotPackage, input$upload)
    rtracklayer::import(input$upload$datapath, format="bed")
  })
  
  #### GRange from web with gscores ####
  gsObject <- reactive({
    req(input$annotPackage, input$indOrRange)
    validate(is_smaller(rStart(), rEnd()))
    validate(is_within_range(chromo(), rStart(), rEnd()))

    switch(input$indOrRange,
           individual = gscores(phast(), 
                                GRanges(seqnames=chromo(),
                                        IRanges(rStart():rEnd(), width=1)),
                                pop = input$populations),
           range = gscores(phast(), 
                           GRanges(fullRange()),
                           pop = input$populations)
    )
    
  })
  
  #### GRange from bed file with gscores #####
  gsUpload <- reactive({
    gscores(phast(), uploadedBed(), pop = input$populations)
  })

  
  ################## UI HIDE AND SHOW  ##################
  
  observe({
    if(input$webOrBed == 'web'){
      shinyjs::show("webOptions")
      shinyjs::hide("upload")
      shinyjs::show("printGsWeb")
      shinyjs::hide("printGsBed")
      shinyjs::show("dwn_web_bed")
      shinyjs::show("dwn_web_csv")
      shinyjs::hide("dwn_bed_bed")
      shinyjs::hide("dwn_bed_csv")
    } else {
      shinyjs::hide("webOptions")
      shinyjs::show("upload")
      shinyjs::hide("printGsWeb")
      shinyjs::show("printGsBed")
      shinyjs::hide("dwn_web_bed")
      shinyjs::hide("dwn_web_csv")
      shinyjs::show("dwn_bed_bed")
      shinyjs::show("dwn_bed_csv")
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
      radioButtons("indOrRange", "Choose between range or individual scores",
                   choices = list("Range" = "range", "Individual" = "individual"))
    )
  })
  
  ### render selectinput with population options
  output$pop <- renderUI({
    req(input$annotPackage)
    selectInput("populations", "Select an available population", multiple = TRUE,
                choices = populations(phast()))
  })
  
  
  ################## OUTPUT TEXT AND TABLES ############################## 

  ### Annotation Package
  output$phastInfo <- renderPrint({
    req(input$annotPackage)
    phast()
  })
  
  ### Citation
  output$citation <- renderPrint({
    req(input$annotPackage)
    citation(phast())
  })
  
  ### Web Datatable
  output$printGsWeb <- DT::renderDataTable({
    if(input$webOrBed != "web" || input$annotPackage == "") return ()
    as.data.table(gsObject())
  })
  
  ### Bed Datatable
  output$printGsBed <- DT::renderDataTable({
    req(input$upload)
    if(input$webOrBed != "bed" ) return ()
    as.data.table(gsUpload())
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
