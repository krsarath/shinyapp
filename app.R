library(shiny)
library(shinyBS) 
library(shinyjs)  
library(DT)       
library(ggplot2)  
library(shinycssloaders)
library(shinythemes)
library(shinydashboard)
library(bslib)

ui <- fluidPage(
  
  titlePanel("DENOSeq"),
  
  tabsetPanel(
    tabPanel("Upload",
               fileInput("upload", "Input the count file", accept = c(".csv", ".txt",".xlsx")),
             actionButton("reset","Reset"),
             tableOutput("data_preview")),    
             
    tabPanel("Differential expression",
    selectInput("diffrential_expression", "Select Differential Expression Method", choices = c("DESeq2", "EdgeR", "NOISeq")),
    actionButton("run","Run"),
    tableOutput("deresult")),
    
    tabPanel("Gene Ontology",
             fileInput("upload", "Input Annotation file"),
    actionButton("reset","Reset")),
    plotOutput("goresult")),
  
      mainPanel(
       ))
      

server <- function(input, output) {
  output$data_preview <- renderTable({
    req(input$upload)
    df <- read.csv(input$upload$datapath)
    
  })
  
  output$deresult <- renderTable({
    df <- read.delim(deresult)
  })
  output$goresult <- renderPlot({
    plot(goresult)
  })
  
  
}

shinyApp(ui, server)
