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
        br(),
        selectInput("organism", "Organism",
                    choices = c("All" = "All", unique(options$Organism)),
                    selected = "All"),
        selectInput("category", "Category",
                    choices = c("All" = "All", unique(options$Category)),
                    selected = "All"),
        tags$div(id="cssref",
                 selectInput("annotPackage", "Annotation Package",
                             choices = NULL)),
        uiOutput("pop"),
        br(), 
        br(),
        radioButtons("webOrBed", "Genomic Coordinates",
                     choices = list("Manually" = "web", "Uploading BED file" = "bed")),
        uiOutput("webOptions"),
        fileInput("upload", "Select your Bed format file"),
        br(),
        fluidRow(
            column(
                width = 12, align = "center", 
                actionButton("run", "RUN", icon = icon("cog"), 
                             class = "run-btn", width = "60%"),
                actionButton("quit", "QUIT", icon = icon("times"),
                             class = "cls-btn" , width = "60%")) 
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