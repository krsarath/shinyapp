library(shiny)
library(bslib)
ui <- page_sidebar(
  title = "DEnoseq",
  tabsetPanel("one"),
  sidebar = sidebar("Upload"),
  card()
  
  
)
server <- function(input, output) {
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
title <- tags$a(href='https://drive.google.com/file/d/1TN9MYCql5dqtIn74M3mDrEom9UkN7B8q/view?usp=drive_link',
                tags$img(src="denoseq.png", height = '50', width = '50'),
                'DENOSEQ', target="_blank")
title1 <- tags$a(href )
tags$img(srcfile("/home/sarath/R/app/Denoseq/denoseq.png"))
title = title1,titleWidth = 300),
tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom1.css")
fileInput("upload", "Input count file ", accept = c(".csv", ".txt",".xlsx")),
sidebarMenuOutput( tableOutput("contents")),
actionButton("clear","Clear")
sidebarMenu(
  "Differential Expression",
  selectInput("","Select Differential Expression Method", choices = c("DESeq2","EdgeR","NOISeq")),
  actionButton("run_analysis", "Run Analysis")),
output$contents <- renderTable({
  req(input$upload)
  df <- read.csv(input$upload$datapath)
  
  
  sidebarMenu(
    "Upload",
    fileInput("upload", "Input count file ", accept = c(".csv", ".txt",".xlsx")),
    sidebarMenuOutput( tableOutput("contents")),
    actionButton("clear","Clear")),
  
  sidebarMenu(
    "Differential Expression",
    selectInput("","Select Differential Expression Method", choices = c("DESeq2","EdgeR","NOISeq")),
    actionButton("run_analysis", "Run Analysis")),
  shinyDashboardThemes(
    theme = "grey_dark"),
  