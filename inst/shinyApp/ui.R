ui <- fluidPage(
    theme = shinythemes::shinytheme("spacelab"),
    shinyjs::useShinyjs(),
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    titlePanel(div(h2("GenomicScores WebApp", align="left"),
               tags$img(src="GenomicScores.png", align="right", height=75, width=75)),
               windowTitle = "GenomicScores"),
    
    sidebarLayout(
        
        sidebarPanel(
            selectInput("annotPackage", "Select a GScores object",
                        choices = c("Choose an installed annotation package" = "", 
                                    avAnnotations())),
            uiOutput("pop"),
            radioButtons("webOrBed", "Input genomic coordinates",
                         choices = list("Manually" = "web", "Uploading BED file" = "bed")),
            uiOutput("webOptions"),
            fileInput("upload", "Upload your Bed format file"),
            width = 3
            
        ),

        mainPanel(
            tabsetPanel(type="tabs",
                        tabPanel("GScore",
                                 fluidRow(id="info",
                                          column(6,
                                                 verbatimTextOutput("annotPackageInfo")
                                          ),
                                          column(6,
                                                 verbatimTextOutput("citation")
                                          )
                                 ),
                                 shinycustomloader::withLoader(DT::dataTableOutput("printGsWeb")),
                                 downloadButton("dwn_web_bed", "Download BED"),
                                 downloadButton("dwn_web_csv", "Download CSV"),
                                 shinycustomloader::withLoader(DT::dataTableOutput("printGsBed")),
                                 downloadButton("dwn_bed_bed", "Download BED"),
                                 downloadButton("dwn_bed_csv", "Download CSV")
                                 ),
                        tabPanel("About",
                                 includeMarkdown("about.md")),
                        tabPanel("Session Info",
                                 verbatimTextOutput("sessionInfo"))),
            width = 9
        )
    )
)
