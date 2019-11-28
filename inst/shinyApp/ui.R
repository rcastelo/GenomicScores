ui <- fluidPage(
    theme = shinytheme("spacelab"),
    useShinyjs(),
    titlePanel(h2("SHINY GSCORE", align="center"), windowTitle = "SHINY GSCORE"),
    
    sidebarLayout(
        
        sidebarPanel(
            selectInput("selectAnn", "Select an Annotation Package",
                        choices = c("Choose an installed annotation package" = "", 
                                    avAnnotations())),
            uiOutput("pop"),
            radioButtons("webOrBed", "Web parameters or BED file?",
                         choices = list("Web" = "web", "BED file" = "bed")),
            uiOutput("webOptions"),
            fileInput("upload", "Upload your Bed format file"),
            width = 3
            
        ),

        mainPanel(
            tabsetPanel(type="tabs",
                        tabPanel("GScore",
                                 fluidRow(id="info",
                                          column(6,
                                                 verbatimTextOutput("phastInfo")
                                          ),
                                          column(6,
                                                 verbatimTextOutput("citation")
                                          )
                                 ),
                                 
                                 withLoader(DT::dataTableOutput("printGsWeb")),
                                 downloadButton("downGscoreWeb", "Download BED"),
                                 downloadButton("downGscoreWebCsv", "Download CSV"),
                                 withLoader(DT::dataTableOutput("printGsBed")),
                                 downloadButton("downGscoreBed", "Download BED"),
                                 downloadButton("downGscoreBedCsv", "Download CSV")
                                 ),
                        tabPanel("About",
                                 includeMarkdown("about.md")),
                        tabPanel("Session Info",
                                 verbatimTextOutput("sessionInfo"))),
            width = 9
        )
    )
)