
library(shiny)
method <- c('DESeq2','EdgeR','NOISeq')
ui <- fluidPage(
# Application title
  
  
       navbarPage("DENOSeq",
          navbarMenu("Method",
                     tabPanel("DESeq2"),
                     tabPanel("EdgeR"),
                     tabPanel("NOISeq")),
          
        navbarMenu("Visualize",
                   tabPanel("Heatmap"),
                   tabPanel("MA plot"),
                   tabPanel("MD plot"),
                   tabPanel("MDS plot"),
                   tabPanel("DE plot"),
                   tabPanel("Venn diagram"),
                   tabPanel("Bubble plot"),
                   tabPanel("Volcano plot"),
                   tabPanel("Dispersion plot"),
                   tabPanel("PCA")),
       navbarMenu("Results",
                  tabPanel("View Results")),
       navbarMenu("Download",
         tabPanel("Download"))),
       
            sidebarLayout(
              sidebarPanel(
                fileInput("upload","Input the count file", accept = c(".csv", ".txt")),
           fileInput("upload","Sample info",accept = c(".csv", ".txt")),
           fileInput("upload","Annotation file")),
        mainPanel(
  plotOutput("distPlot"))))
  

        # Show a plot of the generated distribution
        

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


