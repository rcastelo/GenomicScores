ui <- dashboardPage(
    title = "Genomic Scores",
    
    dashboardHeader(
        tags$li(class = "dropdown",
                tags$div(id = "app_title", "GenomicScores")
        ),
        title = tags$img(src="GenomicScores.png", height=75, width=75),
        titleWidth = 300
    ),
    
    dashboardSidebar(
        width = 300,
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
            uiOutput("css.apkgs")
        ),
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
            column(
                width = 12, align = "left", 
                actionButton("run", "Run", icon = icon("cog"), class = "run-btn"),
                actionButton("quit", "Quit", icon = icon("times"), class = "cls-btn")) 
            )
            
    ),
    
    dashboardBody(
        shinyjs::useShinyjs(),
        fluidRow(
            tabBox(
                width = 12,
                tabPanel("GScore",
                         fluidRow(id="info",
                                  column(6,
                                         withLoader(verbatimTextOutput("annotPackageInfo"))
                                  ),
                                  column(6,
                                         withLoader(verbatimTextOutput("citation"))
                                  )
                         ),
                         withLoader(DT::dataTableOutput("printGs")),
                         uiOutput("down_btn")
                ),
                tabPanel("About",
                         includeMarkdown("about.md")
                         ),
                tabPanel("Session Info",
                         verbatimTextOutput("sessionInfo")))
        )
    )
)