server <- function(input, output, session) {
  
  ###### values ######
  selectAnn <- reactive(input$selectAnn)
  chromo <- reactive(input$chromo)
  rStart <- reactive(input$rStart)
  rEnd <- reactive(input$rEnd)
  fullRange <- reactive(
    paste0(chromo(),":", rStart(),"-", rEnd())
  )
  
  ###### Annotation package #####
  phast <- reactive({
    if(selectAnn()=="")return()
    get(selectAnn())
  })
  
  ##### Uploaded Bed file #####
  uploadedBed<- reactive({
    req(input$selectAnn, input$upload)
    rtracklayer::import(input$upload$datapath, format="bed")
  })
  
  #### GRange from web with gscores ####
  gsObject <- reactive({
    req(input$selectAnn, input$indOrRange)
    req(rEnd()>=rStart())
    switch(input$indOrRange,
           individual = gscores(phast(), GRanges(seqnames=chromo(), IRanges(rStart():rEnd(), width=1)),
                                pop = input$populations),
           range = gscores(phast(), GRanges(fullRange()), pop = input$populations)
    )
  })
  
  #### GRange from bed file with gscores #####
  gsUpload <- reactive({
    gscores(phast(), uploadedBed(), pop = input$populations)
  })

  
  ###### UI Options  ######
  
  observe({
    if(input$webOrBed == 'web'){
      shinyjs::show("webOptions")
      shinyjs::hide("upload")
      shinyjs::show("printGsWeb")
      shinyjs::hide("text4")
      shinyjs::hide("printGsBed")
      shinyjs::show("downGscoreWeb")
      shinyjs::show("downGscoreWebCsv")
      shinyjs::hide("downGscoreBed")
      shinyjs::hide("downGscoreBedCsv")
    } else {
      shinyjs::hide("webOptions")
      shinyjs::show("upload")
      shinyjs::hide("printGsWeb")
      shinyjs::show("text4")
      shinyjs::show("printGsBed")
      shinyjs::hide("downGscoreWeb")
      shinyjs::hide("downGscoreWebCsv")
      shinyjs::show("downGscoreBed")
      shinyjs::show("downGscoreBedCsv")
    }
  })
  
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
  
  output$pop <- renderUI({
    req(input$selectAnn)
    selectInput("populations", "Select an available population", multiple = TRUE,
                choices = populations(phast()))
  })
  
  ###### Output functions   ###### 

  output$phastInfo <- renderPrint({
    req(input$selectAnn)
    phast()
  })
  
  output$citation <- renderPrint({
    req(input$selectAnn)
    citation(phast())
  })
  
  output$printGsWeb <- DT::renderDataTable({
    if(input$webOrBed != "web" ) return ()
    as.data.table(gsObject())
  })
  
  output$printGsBed <- DT::renderDataTable({
    req(input$upload)
    if(input$webOrBed != "bed" ) return ()
    as.data.table(gsUpload())
  })
  
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })
  
  ###### Download Buttons #######
  
  output$downGscoreWeb <- downloadFile(gsObject(), "bed")
  
  output$downGscoreWebCsv <- downloadFile(gsObject(), "csv")
  
  output$downGscoreBed <- downloadFile(gsUpload(), "bed")
  
  output$downGscoreBedCsv <- downloadFile(gsUpload(), "csv")

}
