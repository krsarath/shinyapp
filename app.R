library(shiny)
library(shinyBS) 
library(shinyjs)  
library(DT)       
library(ggplot2)  
library(shinycssloaders)
library(shinythemes)
library(shinydashboard)
library(bslib)

method <- c('DESeq2', 'EdgeR', 'NOISeq')

ui <- fluidPage(
  themeSelector(),
  theme = shinytheme("superhero"),
  titlePanel("DENOSeq"),
  tabsetPanel(
    tabPanel("Upload"),
    tabPanel("Differential expression",
    selectInput("Diffrential expression", "Select Differential Expression Method", choices = c("DESeq2", "EdgeR", "NOISeq"))),
    tabPanel("Gene Ontology")
  ),
  sidebarLayout(
      sidebarPanel(
       fileInput("upload", "Input the count file", accept = c(".csv", ".txt",".xlsx")),
       fileInput("upload", "Input Annotation file")),
  
      mainPanel(
        tableOutput("contents"))))
      

server <- function(input, output) {
  output$contents <- renderTable({
    req(input$upload)
    df <- read.csv(input$upload$datapath)
    
  })
  
}

shinyApp(ui, server)
