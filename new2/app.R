#new2
library(shiny)
library(shinyBS)
library(shinyjs)
library(DT)
library(ggplot2)
library(shinycssloaders)
library(shinythemes)
library(shinydashboard)
library(bslib)

options(shiny.maxRequestSize = 50 * 1024^2)
ui <- navbarPage(
  title = "ShinyDE",
  theme = shinytheme("simplex"),
  tabPanel("Upload",
           sidebarLayout(
             sidebarPanel(
               fileInput("upload_counts", "Input the count file", accept = c(".csv", ".txt", ".xlsx")),
               actionButton("reset_counts", "Reset"),
               tableOutput("data_preview")
             ),
             mainPanel()
           )
  ),
  tabPanel("Differential expression",
           sidebarLayout(
             sidebarPanel(
               selectInput("diffrential_expression", "Select Differential Expression Method", choices = c("DESeq2", "EdgeR", "NOISeq")),
               selectInput("controls", "Select Control Samples", choices = NULL, multiple = TRUE),
               selectInput("treatments", "Select Treatment Samples", choices = NULL, multiple = TRUE),
               selectInput("plot_type", "Visualization", choices = c("pca", "heatmap", "Functional and Pathway Enrichment", "Volcano Plots", "MA Plot", "Dispersion Plot"), selected = "heatmap"),
               actionButton("run_diff_exp", "Run")
             ),
             mainPanel(
               tableOutput("deresult"),
               plotOutput("diff_exp_plot")
             )
           )
  ),
  tabPanel("Gene Ontology",
           sidebarLayout(
             sidebarPanel(
               fileInput("upload_annotation", "Input Annotation file"),
               actionButton("reset_annotation", "Reset"),
               actionButton("run_go_analysis", "Run GO Analysis")
             ),
             mainPanel(
               plotOutput("goresult")
             )
           )
  )
)
server <- function(input, output, session) {
  
  # Data preview
  output$data_preview <- renderTable({
    req(input$upload_counts)
    df <- read.csv(input$upload_counts$datapath)
    head(df) # Display only the head for preview
  })
  
  # Update sample choices based on uploaded data
  observe({
    req(input$upload_counts)
    df <- read.csv(input$upload_counts$datapath)
    samples <- colnames(df)
    updateSelectInput(session, "controls", choices = samples)
    updateSelectInput(session, "treatments", choices = samples)
  })
  
  # Run differential expression analysis
  observeEvent(input$run_diff_exp, {
    req(input$upload_counts)
    req(input$diffrential_expression == "DESeq2")
    
    # Read count data
    count_data <- read.csv(input$upload_counts$datapath, header = TRUE, row.names = 1)
    
    # Get user-selected control and treatment samples
    controls <- input$controls
    treatments <- input$treatments
    
    # Create sample info dataframe based on user selection
    sample_info <- data.frame(
      sampleID = c(controls, treatments),
      Treatment = c(rep("control", length(controls)), rep("treated", length(treatments)))
    )
    rownames(sample_info) <- sample_info$sampleID
    
    # Filter count data to include only selected samples
    countdata <- count_data[, sample_info$sampleID, drop = FALSE]
    
    # Filtering count data
    countdata <- countdata[which(rowMeans(!is.na(countdata)) > 0.5), ]
    countdata <- countdata %>% mutate(across(where(is.numeric), ~ replace(., is.na(.), median(., na.rm = TRUE))))
    
    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = round(countdata), colData = sample_info, design = ~ Treatment)
    dds <- DESeq(dds)
    deseq_result <- results(dds)
    deseq_result <- as.data.frame(deseq_result)
    
    # Debug: Print a summary of the deseq_result
    print(summary(deseq_result))
    
    deseq_result_ordered <- deseq_result[order(deseq_result$pvalue),]
    write.csv(deseq_result_ordered, 'deresult.csv')
    
    output$deresult <- renderTable({
      head(deseq_result_ordered)
    })
    
    # Plots
    output$diff_exp_plot <- renderPlot({
      plot_type <- input$plot_type
      
      if (plot_type == "Volcano Plots") {
        # Volcano Plot
        resLFC <- lfcShrink(dds, coef = "Treatment_treated_vs_control", type = 'ashr')
        resLFC <- as.data.frame(resLFC)
        resLFC$diffexpressed <- NA
        resLFC$diffexpressed[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05] <- "UP"
        resLFC$diffexpressed[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05] <- "DOWN"
        resLFC$delabel <- NA
        
        ggplot(data = resLFC, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
          geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
          geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
          geom_point(size = 2) +
          geom_text_repel() +
          scale_color_manual(values = c('blue', 'red'), labels = c("Downregulated", "Not significant", "Upregulated")) +
          theme(text = element_text(size = 20))
        
      } else if (plot_type == "pca") {
        vsd <- vst(dds, blind = FALSE)
        plotPCA(vsd, intgroup = "Treatment")
        
      } else if (plot_type == "heatmap") {
        vsd <- vst(dds, blind = FALSE)
        sampleDists <- dist(t(assay(vsd)))
        sampleDistMatrix <- as.matrix(sampleDists)
        colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
        pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
                 clustering_distance_cols = sampleDists, col = colors)
        
      } else if (plot_type == "MA Plot") {
        plotMA(deseq_result, ylim = c(-12, 12))
        
      } else if (plot_type == "Dispersion Plot") {
        plotDispEsts(dds)
      }
    })
  })
  
  # Gene Ontology result plot
  observeEvent(input$run_go_analysis, {
    req(input$upload_annotation)
    req(input$upload_counts)
    
    # Read DESeq results
    deseq_result <- read.csv('deresult.csv', row.names = 1)
    
    # Extract significantly upregulated genes
    sig_genes <- rownames(deseq_result[deseq_result$log2FoldChange > 0.5,])
    
    # Perform GO enrichment analysis
    GO_results <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")
    
    # Render bar plot for GO enrichment results
    output$goresult <- renderPlot({
      barplot(GO_results, showCategory = 15)
    })
    
    # Save bar plot as PNG
    png("GO_results.png", res = 250, width = 1400, height = 1800)
    barplot(GO_results, showCategory = 15)
    dev.off()
  })
  
  # Reset input selections
  observeEvent(input$reset_counts, {
    updateSelectInput(session, "controls", choices = NULL)
    updateSelectInput(session, "treatments", choices = NULL)
  })
  
  observeEvent(input$reset_annotation, {
    # Placeholder for resetting the Gene Ontology input
  })
}

shinyApp(ui, server)
