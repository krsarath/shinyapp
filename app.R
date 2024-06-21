
library(shiny)
library(shinyBS)
library(shinyjs)
library(DT)
library(ggplot2)
library(shinycssloaders)
library(shinythemes)
library(shinydashboard)
library(bslib)
library(dashboardthemes)

options(shiny.maxRequestSize = 50 * 1024^2)

ui <- dashboardPage(skin = "red",
                    dashboardHeader(title = span(img(src="denose.png", width = 200))),
                    
                    dashboardSidebar(
                      sidebarMenu(
                        id = "tabs",
                        menuItem("Upload", tabName = "upload", icon = icon("upload")),
                        menuItem("Differential Expression", tabName = "differential_expression", icon = icon("bar-chart")),
                        menuItem("Gene Ontology", tabName = "gene_ontology", icon = icon("list-alt")),
                        conditionalPanel(
                          condition = "input.tabs == 'differential_expression'",
                          selectInput("diffrential_expression", "Select Differential Expression Method", choices = c("DESeq2", "EdgeR", "NOISeq")),
                          uiOutput("controls_ui"),
                          uiOutput("treatments_ui"),
                          actionButton("run", "Run")
                        )
                      )
                    ),
                    
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "upload",
                                fluidRow(
                                  box(title = "Upload Counts File", width = 12,
                                      fileInput("upload_counts", "Input the count file", accept = c(".csv", ".txt", ".xlsx")),
                                      actionButton("reset_upload", "Reset"),
                                      checkboxInput("show_full_data", "Show full data", value = FALSE),
                                      tableOutput("data_preview")
                                  )
                                )
                        ),
                        tabItem(tabName = "differential_expression",
                                fluidRow(
                                  box(title = "Differential Expression Results", width = 12,
                                      tableOutput("deresult")
                                  )
                                )
                        ),
                        tabItem(tabName = "gene_ontology",
                                fluidRow(
                                  box(title = "Gene Ontology Analysis", width = 12,
                                      fileInput("upload_annotation", "Input Annotation file"),
                                      actionButton("reset_go", "Reset"),
                                      plotOutput("goresult")
                                  )
                                )
                        )
                      )
                    )
)

server <- function(input, output, session) {
  data <- reactiveVal(NULL)
  
  observeEvent(input$upload_counts, {
    req(input$upload_counts)
    df <- read.csv(input$upload_counts$datapath, row.names = 1)
    data(df)
  })
  
  output$data_preview <- renderTable({
    req(data())
    if (input$show_full_data) {
      data()
    } else {
      head(data())
    }
  })
  
  observe({
    req(data())
    samples <- colnames(data())
    
    output$controls_ui <- renderUI({
      selectInput("Controls", "Select Control Samples", choices = samples, multiple = TRUE)
    })
    
    output$treatments_ui <- renderUI({
      selectInput("Treatments", "Select Treatment Samples", choices = samples, multiple = TRUE)
    })
  })
  
  output$deresult <- renderTable({
    req(input$run, data(), input$Controls, input$Treatments, input$diffrential_expression)
    # Placeholder for differential expression results
    # Replace with actual DE analysis code
    data.frame(Sample = c("Gene1", "Gene2"), Log2FoldChange = c(1.2, -0.5), PValue = c(0.05, 0.01))
  })
  
  output$goresult <- renderPlot({
    req(input$upload_annotation)
    # Placeholder for GO result plot
    # Replace with actual GO analysis and plotting code
    plot(cars)
  })
}

shinyApp(ui, server)
