library(shiny)
library(shinyBS)
library(shinyjs)
library(DT)
library(ggplot2)
library(shinycssloaders)
library(shinythemes)
library(shinydashboard)
library(bslib)
library(readxl)

ui <- dashboardPage(
  dashboardHeader(title = "Denoseq"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload", tabName = "upload", icon = icon("file-upload")),
      menuItem("Differential Expression", tabName = "differential_expression", icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "upload",
              fileInput("upload", "Input count file", accept = c(".csv", ".txt", ".xlsx")),
              checkboxInput("header", "Header", TRUE),
              radioButtons("sep", "Separator",
                           choices = c(Comma = ",",
                                       Semicolon = ";",
                                       Tab = "\t"),
                           selected = ","),
              radioButtons("quote", "Quote",
                           choices = c(None = "",
                                       "Double Quote" = '"',
                                       "Single Quote" = "'"),
                           selected = '"'),
              DT::dataTableOutput("contents"),
              actionButton("clear", "Clear")
      ),
      tabItem(tabName = "differential_expression",
              selectInput("diff_method", "Select Differential Expression Method", choices = c("DESeq2", "EdgeR", "NOISeq")),
              actionButton("run_analysis", "Run Analysis")
      )
    )
  )
)

server <- function(input, output, session) {
  output$contents <- DT::renderDataTable({
    req(input$upload)
    
    file_ext <- tools::file_ext(input$upload$name)
    
    if (file_ext == "csv" || file_ext == "txt") {
      df <- read.csv(input$upload$datapath,
                     header = input$header,
                     sep = input$sep,
                     quote = input$quote)
    } else if (file_ext == "xlsx") {
      df <- read_excel(input$upload$datapath)
    } else {
      return(NULL)
    }
    
    DT::datatable(df, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  observeEvent(input$clear, {
    session$sendCustomMessage(type = 'resetFileInput', message = 'upload')
  })
  
  # Custom message handler to reset file input
  session$onFlushed(function() {
    session$sendCustomMessage(type = 'resetFileInputHandler', message = 'upload')
  })
}

shinyApp(ui, server)
