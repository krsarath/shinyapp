#shiny vista
library(shiny)
library(shinythemes)
library(DESeq2)
library(BiocManager)
library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(NOISeq)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(plotly)
library(VennDiagram)
library(ggVennDiagram)

options(shiny.maxRequestSize = 50 * 1024^2)

ui <- navbarPage(
  title = span(img(src="shiny-vista.png", width = 120)),
  theme = shinytheme("simplex"),
  tabPanel("Upload",
           sidebarLayout(
             sidebarPanel(
               fileInput("upload_counts", "Input the count file", accept = c(".csv",".xlsx"))
             ),
             mainPanel(tableOutput("data_preview"))
           )
  ),
  tabPanel("Differential expression",
           sidebarLayout(
             sidebarPanel(
               selectInput("diffrential_expression", "Select Differential Expression Method", choices = c("DESeq2", "EdgeR", "NOISeq")),
               selectInput("controls", "Select Control Samples", choices = NULL, multiple = TRUE),
               selectInput("treatments", "Select Treatment Samples", choices = NULL, multiple = TRUE),
               selectInput("plot_type", "Visualization", choices = c("PCA", "Dispersion Plot", "MA Plot", "Heatmap", "Volcano Plot", "MDS", "BCV", "MD plot", "DE plot")),
               actionButton("run_diff_exp", "Run")
             ),
             mainPanel(
               tableOutput("deresult"),
               plotOutput("diff_exp_plot"),
               downloadButton("download_deresult", "Download"),
               downloadButton("download_plot",  "Download Plot")
             )
           )
  ),
  tabPanel("GO bubble Plot",
           sidebarLayout(
             sidebarPanel(
               fileInput("deresult", "Upload DE Result", accept = c( ".csv", ".xlsx")),
               fileInput("go_mappings", "Upload GO Terms", accept = c(".csv", ".xlsx")),
               actionButton("run_go", "GO plot")
             ),
             mainPanel(
               plotOutput("go_plot"),
               downloadButton("download_plot", "Download Plot")
             )
           )
  ),
  tabPanel("Venn Diagram",
        sidebarLayout(
          sidebarPanel( fileInput("file1", "Upload csv/xlsx list 1", accept = c(".csv", ".xlsx")),
                        textInput("name1", "listname", value = "list1"),
                        fileInput("file2", "list 2", accept = c(".csv", ".xlsx")),
                        textInput("name2", "listname", value = "list2"),
                        fileInput("file3", "list 3", accept = c(".csv", ".xlsx")),
                        textInput("name3", "listname", value = "list3"),
                        fileInput("file4", "list 4", accept = c(".csv", ".xlsx")),
                        textInput("name4", "listname", value = "list4"),
                        fileInput("file5", "list 5", accept = c(".csv", ".xlsx")),
                        textInput("name5", "listname", value = "list5"),
                        fileInput("file6", "list 6", accept = c(".csv", ".xlsx")),
                        textInput("name6", "listname", value = "list6"),
                        radioButtons("regulationType", "Select Regulation Type:",
                                     choices = list("Upregulated" = "upregulated",
                                                    "Downregulated" = "downregulated"),
                                     selected = "upregulated"),
                        actionButton("plot", "Generate")
          ),
          mainPanel(
            plotOutput("venn_plot"),
            downloadButton("download_plot", "Download Plot")
          )
        )   )
)

server <- function(input, output, session) {
  
  # Data preview
  output$data_preview <- renderTable({
    req(input$upload_counts)
    if (stringr::str_ends(input$upload_counts$datapath, "csv")) {
      df <- read.csv(input$upload_counts$datapath)
    } else if (stringr::str_ends(input$upload_counts$datapath, "xlsx")) {
      df <- readxl::read_excel(input$upload_counts$datapath)
    } 
    
    #df <- read.csv(input$upload_counts$datapath)
    head(df) # Display only the head for preview
  })
  
  # Update sample choices based on uploaded data
  observe({
    req(input$upload_counts)
    if (stringr::str_ends(input$upload_counts$datapath, "csv")) {
      df <- read.csv(input$upload_counts$datapath)
    } else if (stringr::str_ends(input$upload_counts$datapath, "xlsx")) {
      df <- readxl::read_excel(input$upload_counts$datapath)
    } 
    #df <- read.csv(input$upload_counts$datapath)
    samples <- colnames(df)
    updateSelectInput(session, "controls", choices = samples)
    updateSelectInput(session, "treatments", choices = samples)
  })
  
  # Run differential expression analysis for DESeq2
  observeEvent(input$run_diff_exp, {
    req(input$upload_counts)
    if (input$diffrential_expression == "DESeq2") {
      # Read count data
      if (stringr::str_ends(input$upload_counts$datapath, "csv")) {
        count_data <- read.csv(input$upload_counts$datapath)
      } else if (stringr::str_ends(input$upload_counts$datapath, "xlsx")) {
        count_data <- readxl::read_excel(input$upload_counts$datapath)
      } 
      #count_data <- read.csv(input$upload_counts$datapath, header = TRUE, row.names = 1)
      
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
      
      # Filtering count data (optional steps, you can customize as per your data)
      countdata <- countdata[which(rowMeans(!is.na(countdata)) > 0.5), ]
      countdata <- countdata %>% mutate(across(where(is.numeric), ~ replace(., is.na(.), median(., na.rm = TRUE))))
      
      # Create DESeqDataSet
      dds <- DESeqDataSetFromMatrix(countData = round(countdata), colData = sample_info, design = ~ Treatment)
      
      dds <- DESeq(dds)
      deseq_result <- results(dds)
      deseq_result <- as.data.frame(deseq_result)
      deseq_result <- deseq_result[order(deseq_result$pvalue),]
      deseq_result <- deseq_result %>%
        mutate(Regulation = ifelse(log2FoldChange > 1, "upregulated", "downregulated"))
      output$deresult <- renderTable({
        head(deseq_result)
      })
      
      # Downloadable csv of selected dataset
      output$download_deresult <- downloadHandler(
        filename = function() {
          paste("deresult", ".csv", sep = "")
        },
        content = function(file) {
          write.csv(deseq_result, file, row.names = TRUE)
        }
      )
      
      # Plots
      output$diff_exp_plot <- renderPlot({
        plot_type <- input$plot_type
        
        if (plot_type == "Volcano Plot") {
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
          
        } else if (plot_type == "PCA") {
          vsd <- vst(dds, blind = FALSE)
          plotPCA(vsd, intgroup = "Treatment")
          
        } else if (plot_type == "Heatmap") {
          top_hits <- deseq_result[order(deseq_result$padj),][1:30,]
          top_hits <- row.names(top_hits)
          top_hits
          rld <- rlog(dds, blind = F)
          pheatmap(assay(rld)[top_hits,])
        } else if (plot_type == "MA Plot") {
          # Prepare data for ggmaplot
          ma_data <- data.frame(
            name = rownames(deseq_result),
            baseMean = deseq_result$baseMean,
            log2FoldChange = deseq_result$log2FoldChange,
            pvalue = deseq_result$pvalue,
            padj = deseq_result$padj
          )
          
          # Create MA plot with ggmaplot
          ggmaplot(ma_data, 
                   #main = expression("Control" %->% "Treated"),
                   fdr = 0.05, fc = 2, size = 0.4,
                   #palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames = as.vector(ma_data$name),
                   #legend = "top", top = 20,
                   #font.label = c("bold", 11),
                   #font.legend = "bold",
                   #font.main = "bold",
                   ggtheme = ggplot2::theme_linedraw()
          )
          
        } else if (plot_type == "Dispersion Plot") {
          plotDispEsts(dds)
        }
        
        
      })
    }
  })
  # Run differential expression analysis for EdgeR
  observeEvent(input$run_diff_exp, {
    req(input$upload_counts)
    if (input$diffrential_expression == "EdgeR") {
      # Read count data
      if (stringr::str_ends(input$upload_counts$datapath, "csv")) {
        count_data <- read.csv(input$upload_counts$datapath)
      } else if (stringr::str_ends(input$upload_counts$datapath, "xlsx")) {
        count_data <- readxl::read_excel(input$upload_counts$datapath)
      } 
      #count_data <- read.csv(input$upload_counts$datapath, header = TRUE, row.names = 1)
      
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
      
      # Filtering count data (optional steps, you can customize as per your data)
      countdata <- countdata[which(rowMeans(!is.na(countdata)) > 0.5), ]
      countdata <- countdata %>% mutate(across(where(is.numeric), ~ replace(., is.na(.), median(., na.rm = TRUE))))
      group <- sample_info$Treatment
      group <- factor(group)
      table(group)
      dge_list <- DGEList(counts = countdata, genes = row.names(countdata), group = group)  
      dge_list <- dge_list[filterByExpr(dge_list), , keep.lib.sizes = F]
      dge_list <- calcNormFactors(dge_list)
      design <- model.matrix(~0 + group)
      colnames(design) <- levels(group)
      est_disp <- estimateDisp(dge_list, design, robust = TRUE)
      qlfit <- glmQLFit(dge_list, design, robust = TRUE)
      qlftest <- glmQLFTest(qlfit)
      qlftest <- as.data.frame(topTags(qlftest,n = nrow(dge_list)))
      qlftest <- qlftest[order(qlftest$PValue),]
      qlftest <- qlftest %>%
        mutate(Regulation = ifelse(logFC > 1, "upregulated", "downregulated"))
      output$deresult <- renderTable({
        head(qlftest)
      })
      
      # Downloadable csv of selected dataset
      output$download_deresult <- downloadHandler(
        filename = function() {
          paste("deresult", ".csv", sep = "")
        },
        content = function(file) {
          write.csv(qlftest, file, row.names = TRUE)
        }
      )
      
      # Plots
      output$diff_exp_plot <- renderPlot({
        plot_type <- input$plot_type
        
        if (plot_type == "Volcano Plot") {
          et <- exactTest(est_disp)
          et$table$logFC[is.na(et$table$logFC)] <- 0
          et$table$PValue[is.na(et$table$PValue)] <- 1
          plot(et$table$logFC, -10 * log10(et$table$PValue), main = "Volcano plot", xlab = "logFC", ylab = "-10 * log(P-val)")
          points(et$table$logFC[et$table$logFC > 0 & et$table$PValue < 0.01], -10 * log10(et$table$PValue[et$table$logFC > 0 & et$table$PValue < 0.01]), col = "red")
          points(et$table$logFC[et$table$logFC < 0 & et$table$PValue < 0.01], -10 * log10(et$table$PValue[et$table$logFC < 0 & et$table$PValue < 0.01]), col = "green")
          legend("topright", legend = c("UP", "Down", "Not significant"), col = c("red", "green", "black"), pch = 20)
          
        } else if (plot_type == "MDS") {
          dge_list <- DGEList(counts = countdata, genes = row.names(countdata), group = group)  
          dge_list <- dge_list[filterByExpr(dge_list), , keep.lib.sizes = F]
          dge_list <- calcNormFactors(dge_list)
          plotMDS(dge_list)
          
        } else if (plot_type == "BCV") {
          est_disp <- estimateDisp(dge_list, design, robust = TRUE)
          est_disp$common.dispersion
          est_disp$trended.dispersion
          est_disp$tagwise.dispersion
          plotBCV(est_disp)
          
        } else if (plot_type == "MD plot") {
          et <- exactTest(est_disp)
          top_degs <- topTags(et)
          cpm(dge_list)[rownames(top_degs),]
          summary(de <- decideTests(et, lfc = 1))
          detags <- rownames(dge_list)[as.logical(de)]
          plotMD(et, de.tags = detags)
          abline(h = c(-1, 1), col = "blue")
          
        } else if (plot_type == "Dispersion Plot") {
          plotQLDisp(qlfit)
          
        } else if (plot_type == "Heatmap") {
          logcpm <- cpm(dge_list, log = TRUE)
          heatmap(logcpm)
        }
      })
    }
  })
  
  observeEvent(input$run_diff_exp, {
    req(input$upload_counts)
    if (input$diffrential_expression == "NOISeq") {
      if (stringr::str_ends(input$upload_counts$datapath, "csv")) {
        count_data <- read.csv(input$upload_counts$datapath)
      } else if (stringr::str_ends(input$upload_counts$datapath, "xlsx")) {
        count_data <- readxl::read_excel(input$upload_counts$datapath)
      } 
      # Read count data
      #count_data <- read.csv(input$upload_counts$datapath, header = TRUE, row.names = 1)
      
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
      
      group <- factor(sample_info$Treatment)
      countdata <- filtered.data(countdata, factor = sample_info$Treatment, norm = TRUE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")
      mydata <- readData(data = countdata, factors = sample_info)
      myresults <- noiseq(mydata, factor = "Treatment", k = 0, norm = "tmm", 
                          condition = c("control", "treated"), pnr = 0.2, nss = 5, 
                          v = 0.02, lc = 1, replicates = "technical")
      myresults.deg <- degenes(myresults, q = 0.8, M = NULL)
      myresults.deg1 <- degenes(myresults, q = 0.8, M = "up")
      myresults.deg2 <- degenes(myresults, q = 0.8, M = "down")
      
      # Extracting NOISeq results for display
      results_df <- as.data.frame(myresults@results[[1]])
      results_df <- results_df %>%
        mutate(Regulation = ifelse(M > 1, "upregulated", "downregulated"))
      output$deresult <- renderTable({
        head(results_df)
      })
      
      # Downloadable csv of selected dataset
      output$download_deresult <- downloadHandler(
        filename = function() {
          paste("myresults_deg", ".csv", sep = "")
        },
        content = function(file) {
          write.csv(results_df, file, row.names = TRUE)
        }
      )
      
      # Plots
      output$diff_exp_plot <- renderPlot({
        plot_type <- input$plot_type
        
        if (plot_type == "MD plot") {
          DE.plot(myresults, q = 0.8, graphic = "MD")
          
        } else if (plot_type == "DE plot") {
          DE.plot(myresults, q = 0.8, graphic = "expr", log.scale = TRUE)
          
        } else if (plot_type == "Heatmap") {
          N <- 20
          top_genes <- results_df[order(abs(results_df$prob), decreasing = TRUE), ]
          top_genes <- head(top_genes, N)
          top_gene_names <- rownames(top_genes)
          top_gene_exprs <- countdata[top_gene_names, ]
          pheatmap(top_gene_exprs,
                   main = "Heatmap of Top Significant Genes",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
          
        } 
        
      })
    }
  })
  go_plot_data <- reactive({
    req(input$deresult, input$go_mappings)
    
    result_de <- read.csv(input$deresult$datapath)
    term_GO <- read.csv(input$go_mappings$datapath)
    
    # Merge data
    merged_data <- inner_join(term_GO, result_de, by = "gene_id")
    
    # Calculate the count of genes per GO category
    merged_data <- merged_data %>%
      mutate(Regulation = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
      group_by(category_id, Regulation) %>%
      summarise(count = n(),
                log2FoldChange = mean(log2FoldChange),
                pvalue = mean(pvalue)) %>%
      ungroup()
    
    # Select top N categories by count of genes
    top_categories <- merged_data %>%
      group_by(category_id) %>%
      summarise(total_count = sum(count)) %>%
      top_n(n = 50, wt = total_count) %>%
      pull(category_id)
    
    # Filter data to include only the top N categories
    filtered_data <- merged_data %>% filter(category_id %in% top_categories)
    
    return(filtered_data)
  })
  
  output$go_plot <- renderPlot({
    filtered_data <- go_plot_data()
    
    ggplot(filtered_data, aes(x = -log10(pvalue), y = category_id, size = count, color = Regulation)) +
      geom_point(alpha = 0.7) +
      scale_size_continuous(range = c(2, 12)) +  # Adjust size range as needed
      scale_color_manual(values = c("Up" = "orange", "Down" = "purple")) +  # Color for up/down regulation
      labs(x = "-log10(P-Value)", y = "Significant GOs", size = "Count", color = "Regulation") +
      theme_minimal()  # Adjust theme as per your preference
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("go_plot", ".png", sep = "")
    },
    content = function(file) {
      filtered_data <- go_plot_data()
      
      ggsave(file, plot = ggplot(filtered_data, aes(x = -log10(pvalue), y = category_id, size = count, color = Regulation)) +
               geom_point(alpha = 0.7) +
               scale_size_continuous(range = c(2, 12)) +  # Adjust size range as needed
               scale_color_manual(values = c("Up" = "orange", "Down" = "purple")) +  # Color for up/down regulation
               labs(x = "-log10(P-Value)", y = "Significant GOs", size = "Count", color = "Regulation") +
               theme_light(), width = 12, height = 8, units = "in", device = "png")
    }
  )
  data <- reactive({
    files <- list(input$file1, input$file2, input$file3, input$file4, input$file5, input$file6 )
    valid_files <- files[!sapply(files, is.null)]
    
    data_list <- lapply(valid_files, function(file) {
      ext <- tools::file_ext(file$datapath)
      if (ext == "csv") {
        read_csv(file$datapath)
      } else if (ext == "xlsx") {
        read_excel(file$datapath)
      }
    })
    custom_names <- c(input$name1, input$name2, input$name3, input$name4, input$name5, input$name6)
    names(data_list) <- custom_names[seq_along(data_list)]
    data_list
    
  })
  
  output$venn_plot <- renderPlot({
    req(input$plot)
    
    data_list <- data()
    regulation_type <- input$regulationType
    
    gene_lists <- lapply(data_list, function(df) {
      if ("Regulation" %in% colnames(df) && "gene_id" %in% colnames(df)) {
        df %>%
          filter(Regulation == regulation_type) %>%
          pull(gene_id)
      } else {
        character(0)
      }
    })
    
    names(gene_lists) <- names(data_list)
    
    venn <- ggVennDiagram(gene_lists, label_alpha = 0, category.names = names(data_list)) +
      theme(legend.position = "right") +
      scale_fill_gradient(low = "red", high = "green")
    
    print(venn)
  })
  
}

shinyApp(ui, server)

