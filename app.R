#denoseq
library(shiny)
library(shinythemes)
library(DESeq2)
library(BiocManager)
library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(NOISeq)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(RColorBrewer)
library(plotly)
library(VennDiagram)
library(ggVennDiagram)

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
               selectInput("plot_type", "Visualization", choices = c("PCA", "Dispersion Plot", "MA Plot", "Heatmap", "Volcano Plot", "MDS", "BCV", "MD plot", "DE plot", "Venn plot")),
               actionButton("run_diff_exp", "Run")
             ),
             mainPanel(
               tableOutput("deresult"),
               plotOutput("diff_exp_plot"),
               downloadButton("download_deresult", "Download")
             )
           )
  ),
  tabPanel("Gene Ontology",
           sidebarLayout(
             sidebarPanel(
               fileInput("upload_annotation", "Input Annotation file (Optional)", accept = c(".csv", ".txt", ".xlsx", ".gtf", ".gff")),
               textInput("paste_gene_ids", "Or Paste Gene IDs (comma or newline separated)"),
               selectInput("go_plot", "Choose plot type", choices = c("Dot plot", "Bar plot", "Bubble plot", "Cnet plot")),
               actionButton("reset_annotation", "Reset"),
               actionButton("run_go_analysis", "Run GO Analysis")
             ),
             mainPanel(
               plotOutput("goresult"),
               plotOutput("go_plot")
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
  
  # Run differential expression analysis for DESeq2
  observeEvent(input$run_diff_exp, {
    req(input$upload_counts)
    if (input$diffrential_expression == "DESeq2") {
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
      
      # Filtering count data (optional steps, you can customize as per your data)
      countdata <- countdata[which(rowMeans(!is.na(countdata)) > 0.5), ]
      countdata <- countdata %>% mutate(across(where(is.numeric), ~ replace(., is.na(.), median(., na.rm = TRUE))))
      
      # Create DESeqDataSet
      dds <- DESeqDataSetFromMatrix(countData = round(countdata), colData = sample_info, design = ~ Treatment)
      dds <- DESeq(dds)
      deseq_result <- results(dds)
      deseq_result <- as.data.frame(deseq_result)
      deseq_result <- deseq_result[order(deseq_result$pvalue),]
      
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
          top_hits <- deseq_result[order(deseq_result$padj),][1:10,]
          top_hits <- row.names(top_hits)
          top_hits
          rld <- rlog(dds, blind = F)
          pheatmap(assay(rld)[top_hits,])
        } else if (plot_type == "MA Plot") {
          resLFC <- lfcShrink(dds, coef = 2, type = 'ashr')
          
          plotMA(resLFC, ylim = c(-10, 10))
          
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
      qlftest <- as.data.frame(topTags(qlftest))
      qlftest <- qlftest[order(qlftest$PValue),]
      
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
          logcpm <- cpm(dge_list, log=TRUE)
          heatmap(logcpm)
        }
      })
    }
  })
  
  observeEvent(input$run_diff_exp, {
    req(input$upload_counts)
    if (input$diffrential_expression == "NOISeq") {
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
      
      # Filtering count data (optional steps, you can customize as per your data)
      exprs <- as.matrix(countdata)
      exprs <- apply(exprs, 2, function(x) {
        if (any(is.na(x))) {
          x[is.na(x)] <- median(x, na.rm = TRUE)
        }
        return(x)
      })
      myTMM <- tmm(exprs, long = 1000, lc = 0)
      group <- factor(sample_info$Treatment)
      table(group)
      myfilt <- filtered.data(countdata, factor = group, norm = TRUE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")
      mydata2 <- NOISeq::readData(data = myfilt, factors = group)
      myresults <- noiseq(mydata2, factor = "group", k = 0, norm = "tmm", condition = c("control", "treatment"), pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "technical")
      myresults.deg <- degenes(myresults, q = 0.8, M = NULL)
      myresults.deg1 <- degenes(myresults, q = 0.8, M = "up")
      myresults.deg2 <- degenes(myresults, q = 0.8, M = "down")
      write.csv(myresults, 'myresults.csv', row.names = TRUE)
      
      output$deresult <- renderTable({
        head(myresults)
      })
      
      # Downloadable csv of selected dataset
      output$download_deresult <- downloadHandler(
        filename = function() {
          paste("myresults", ".csv", sep = "")
        },
        content = function(file) {
          write.csv(myresults, file, row.names = TRUE)
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
          top_genes <- myresults.deg[order(abs(myresults.deg$ranking), decreasing = TRUE), ]
          top_genes <- head(top_genes, N)
          top_gene_names <- rownames(top_genes)
          top_gene_exprs <- exprs[top_gene_names, ]
          heatmaply(top_gene_exprs, 
                    main = "Heatmap of Top Significant Genes",
                    xlab = "Samples", 
                    ylab = "Genes",
                    scale = "row",  
                    colors = viridis::viridis(256))
          
        } else if (plot_type == "Venn plot") {
          up <- rownames(myresults.deg1)
          down <- rownames(myresults.deg2)
          all <- rownames(exprs)
          no_expression_genes <- setdiff(all, union(up, down))
          gene_sets <- list(
            UP = up,
            DOWN = down,
            ALL = no_expression_genes
          )
          ggVennDiagram(gene_sets, label_alpha = 0)
        }
      })
    }
  })
}

shinyApp(ui, server)
