library(shiny)
library(shinyBS)  # for modal
library(shinyjs)  # easy javascript functionalities with shiny
library(DT)       # for interactive data tables
library(ggplot2)  #
library(shinycssloaders)
library(shinythemes)
library(shinydashboard)
library(bslib)

method <- c('DESeq2', 'EdgeR', 'NOISeq')

ui <- fluidPage(
  themeSelector(),
  theme = shinytheme("superhero"),
  titlePanel("DENOSeq"),
  dashboardPage(
    dashboardHeader(title = "Home"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Data", tabName = "tab1", icon = icon("th"))
      )
    ),
    dashboardBody(
      shinyjs::useShinyjs(),
      tabItems(
        tabItem(
          tabName = "tab1",
          box(
            width = 5,
            title = "Input Data",
            status = "primary",
            solidHeader = TRUE,
            fileInput("file", "Upload Data"),
            shinyjs::hidden(
              div(
                id = "data_b",
                style = "display:inline-block",
                actionButton("data", "View Data", icon = icon('table'))
              )
            ),
            shinyjs::hidden(
              div(
                id = "plot_b",
                style = "display:inline-block",
                actionButton("plot", "View Plot", icon = icon("bar-chart"))
              )
            )
          )
        )
      ),
      bsModal(
        id = "dataset",
        title = "Diamonds Dataset",
        trigger = "data",
        size = "large",
        withSpinner(dataTableOutput("data_set"))
      ),
      bsModal(
        id = "Plot",
        title = "Plot",
        trigger = "plot",
        size = "large",
        sliderInput(
          inputId = "b",
          label = "Select the bin width value",
          min = 50,
          max = 500,
          value = 100
        ),
        br(),
        withSpinner(plotOutput("plot_gg"))
      ),
      navbarPage("R",
      navbarMenu(
        "Method",
        tabPanel("DESeq2"),
        tabPanel("EdgeR"),
        tabPanel("NOISeq")
      ),
      navbarMenu(
        "Visualize",
        tabPanel("Heatmap", plotOutput("Heatmap")),
        tabPanel("MA plot", plotOutput("MA plot")),
        tabPanel("MD plot", plotOutput("MD plot")),
        tabPanel("MDS plot", plotOutput("MDS plot")),
        tabPanel("DE plot", plotOutput("DE plot")),
        tabPanel("Venn diagram", plotOutput("Venn diagram")),
        tabPanel("Bubble plot", plotOutput("Bubble plot")),
        tabPanel("Volcano plot", plotOutput("Volcano plot")),
        tabPanel("Dispersion plot", plotOutput("Dispersion plot")),
        tabPanel("PCA", plotOutput("PCA"))
      ),
      navbarMenu(
        "Results",
        tabPanel("View Results", tableOutput("View Results"))
      ),
      navbarMenu(
        "Download",
        tabPanel("Download", verbatimTextOutput("Download"))
      )
    )
  )
))

server <- function(input, output) {
  observeEvent(input$file, shinyjs::show("data_b"))
  observeEvent(input$file, shinyjs::show("plot_b"))
  
  data_uploaded <- reactive({
    file1 <- input$file
    if (is.null(file1)) {
      return()
    }
    read.table(
      file = file1$datapath,
      sep = ",",
      header = TRUE,
      stringsAsFactors = TRUE
    )
  })
  
  output$data_set <- renderDataTable(data_uploaded(), options = list(scrollX = TRUE))
  
  output$plot_gg <- renderPlot(
    ggplot(data = data_uploaded()) +
      geom_histogram(binwidth = input$b, aes(x = price)) +
      ggtitle("Diamond Price Distribution") +
      xlab(paste("Diamond Price & binwidth as ", input$b)) +
      ylab("Frequency") +
      theme_minimal() +
      xlim(0, 2500)
  )
}

shinyApp(ui = ui, server = server)