#new
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
  titlePanel("ShinyDE"),
  themeSelector(),
  tabsetPanel(
    tabPanel("Upload",
             fileInput("upload_counts", "Input the count file", accept = c(".csv", ".txt", ".xlsx")),
             actionButton("reset_counts", "Reset"),
             tableOutput("data_preview")
    ),
    tabPanel("Differential expression",
             selectInput("diffrential_expression", "Select Differential Expression Method", choices = c("DESeq2", "EdgeR", "NOISeq")),
             selectInput("controls", "Select Control Samples", choices = NULL),
             selectInput("treatments", "Select Treatment Samples", choices = NULL),
             selectInput("plot_type", "Visualization", choices = c("pca", "tsne", "heatmap", "Functional and Pathway Enrichment", "Volcano Plots"), selected = "heatmap"),
             actionButton("run_diff_exp", "Run"),
             tableOutput("deresult")
    ),
    tabPanel("Gene Ontology",
             fileInput("upload_annotation", "Input Annotation file"),
             actionButton("reset_annotation", "Reset"),
             plotOutput("goresult")
    )
  )
)

server <- function(input, output, session) {
  
  output$data_preview <- renderTable({
    req(input$upload_counts)
    df <- read.csv(input$upload_counts$datapath)
    head(df) # Display only the head for preview
  })
  
  observe({
    req(input$upload_counts)
    df <- read.csv(input$upload_counts$datapath)
    samples <- colnames(df)
    updateSelectInput(session, "controls", choices = samples)
    updateSelectInput(session, "treatments", choices = samples)
  })
  
  output$deresult <- renderTable({
    req(input$upload_counts)
    df <- read.csv(input$upload_counts$datapath)
    # Placeholder for differential expression results, replace with actual DESeq2/EdgeR/NOISeq processing
    head(df) # Display only the head for results preview
  })
  
  output$goresult <- renderPlot({
    req(input$upload_annotation)
    # Placeholder for Gene Ontology results plot, replace with actual GO processing and plotting
    plot(1:10, 1:10) # Example plot
  })
  
  observeEvent(input$reset_counts, {
    updateSelectInput(session, "controls", choices = NULL)
    updateSelectInput(session, "treatments", choices = NULL)
  })
  
  observeEvent(input$reset_annotation, {
    # Placeholder for resetting the Gene Ontology input
  })
}

shinyApp(ui, server)
