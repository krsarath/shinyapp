setwd("/home/sarath/R/app/denoseq")
library(shiny)

method <- c('DESeq2','EdgeR','NOISeq')

ui <- fluidPage(
# Application title
  titlePanel("Denoseq"),
# Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
          fileInput("upload","Input the count file", accept = c(".csv", ".txt")),
          fileInput("upload","Sample info",accept = c(".csv", ".txt")),
          fileInput("upload","Annotation file"))),
        navbarPage("Denoseq",
          navbarMenu("Method",
                     tabPanel("DESeq2"),
                     tabPanel("EdgeR"),
                     tabPanel("NOISeq"))),
                   

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
