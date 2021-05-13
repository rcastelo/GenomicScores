ui <- fluidPage(
    theme = shinythemes::shinytheme("spacelab"),
    shinyjs::useShinyjs(),
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
        uiOutput("css.apkgs")
    ),
    titlePanel(div(h2("GenomicScores WebApp", align="left"),
               tags$img(src="GenomicScores.png", align="right", height=75, width=75)),
               windowTitle = "GenomicScores"),
    
    sidebarLayout(
        
        sidebarPanel(
            selectInput("organism", "Select an Organism",
                        choices = c("All" = "All", unique(options$Organism)),
                        selected = "All"),
            selectInput("category", "Select a Category",
                        choices = c("All" = "All", unique(options$Category)),
                        selected = "All"),
            tags$div(id="cssref", 
                     selectInput("annotPackage", "Select an Annotation Package",
                                 choices = NULL)),
            uiOutput("pop"),
            radioButtons("webOrBed", "Input genomic coordinates",
                         choices = list("Manually" = "web", "Uploading BED file" = "bed")),
            uiOutput("webOptions"),
            fileInput("upload", "Upload your Bed format file"),
            fluidRow(
                actionButton("run", "Run"),
                actionButton("quit", "Quit")),
            width = 3
            
        ),

        mainPanel(
            tabsetPanel(type="tabs",
                        tabPanel("GScore",
                                 fluidRow(id="info",
                                          column(6,
                                                 shinycustomloader::withLoader(verbatimTextOutput("annotPackageInfo"))
                                          ),
                                          column(6,
                                                 shinycustomloader::withLoader(verbatimTextOutput("citation"))
                                          )
                                 ),
                                 shinycustomloader::withLoader(DT::dataTableOutput("printGs")),
                                 uiOutput("down_btn")
                                 ),
                        tabPanel("About",
                                 includeMarkdown("about.md")),
                        tabPanel("Session Info",
                                 verbatimTextOutput("sessionInfo"))),
            width = 9
        )
    )
)
