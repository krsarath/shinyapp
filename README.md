[![R tests](https://github.com/krsarath/shinyapp/actions/workflows/r-test.yml/badge.svg)](https://github.com/krsarath/shinyapp/actions/workflows/r-test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Shiny App](https://img.shields.io/badge/Shiny-vista-blue?logo=RStudio)](https://krsarath.shinyapps.io/shinyvistaApp/)
# Shiny VISTA

# An R Shiny Web-Application for Differential Expression Analysis and Transcriptome data visualization 

This R Shiny web application allows users to perform differential expression analysis and visualize transcriptome data, making it easier to interpret gene expression patterns. The app supports uploading CSV/XLSX files, performing statistical analysis, and generating various plots to aid in the analysis.

- Authors : Sarath kumar R and Dr. J. Sreekumar.
- Website : https://krsarath.shinyapps.io/shinyvistaApp/ shinyvista Ver. 1.0

## Statement of Need

DGE analysis of RNA-seq data requires proficiency in R programming 
and Bioconductor, limiting accessibility for wet-lab researchers 
without bioinformatics support. Most existing web tools support 
only a single statistical method. Shiny Vista provides a unified, 
installation-free interface integrating DESeq2, edgeR, and NOISeq, 
enabling wet-lab biologists, plant scientists, and students to 
perform multi-method DGE analysis and cross-method validation 
without writing any code.

## Features

- **Differential Expression Analysis**: Identify differentially expressed genes between conditions.
- **Data Upload**: Upload transcriptome data in CSV or XLSX format.
- **Interactive Visualizations**:
  - Volcano plots
  - MA plots
  - Heatmaps
  - Venn diagrams
- **Downloadable Results**: Export analysis results and plots.
- **GO Enrichment Analysis**: Bubble plot visualisation of enriched biological processes

## Getting Started

### Prerequisites

- **R** (version 4.x or higher)
- **R Shiny** (version 1.6 or higher)
- Additional R packages: `ggplot2`, `DESeq2`, `edgeR`, `Noiseq`, `dplyr`, `shiny`, `shinyWidgets`, `readxl`, `readr`

### Installation

Clone the repository and navigate to the project directory:

``bash
- git clone https://github.com/krsarath/shinyapp.git
- cd shinyapp

# Install the required packages
```r
install.packages(c("shiny", "ggplot2", "dplyr", "shinyWidgets", 
                   "readxl", "readr", "ggVennDiagram", "pheatmap"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "NOISeq"))
```

#  Running the App

To run the app locally, use the following R command
shiny::runApp('path_to_your_app_directory')



# How to Use

- Upload RNA-seq Count File: Use the 'Upload' tab to select your count file.
- Perform Differential Expression Analysis: Navigate to the 'Differential Expression' tab, choose your settings, and click 'Run'.
- Perform Gene Ontology enrichment visualization: In the 'GO Bubble Plot' tab, upload your DE results and GO files, then generate the plot.
- Create Venn Diagrams: Use the 'Venn Diagram' tab to upload multiple DE result files and compare them visually.
- View Results: Explore the generated plots and tables.
- Download: Export the results and plots for further analysis.


# Example Data

Example datasets are included in the data directory to help you get started.[matrix.counts (copy).zip](url)

## Validation

Validated using the *Drosophila melanogaster* RNA-seq benchmark 
dataset (GEO: GSE18508; Brooks et al., 2011):

- DESeq2: 8,131 DEGs
- edgeR: 7,084 DEGs  
- NOISeq: 2,282 DEGs (conservative high-confidence subset)

## Citation

If you use Shiny Vista, please cite:

> Sarath Kumar R., J. Sreekumar. Shiny Vista: A Web-Based R Shiny 
> Application for Streamlined Differential Expression Analysis and 
> Interactive Visualization of Transcriptome Data. 2026.  
> https://github.com/krsarath/shinyapp

# License

This project is licensed under the MIT License - see the LICENSE file for details.

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/krsarath/shinyapp/blob/main/LICENSE)


# Contact

For any questions or issues, please contact

[![Gmail](https://img.shields.io/badge/Sarath_Kumar_R-DB4437?style=for-the-badge&logo=Gmail&logoColor=white)](mailto:krsarath20@gmail.com) [![LinkedIn](https://img.shields.io/badge/Sarath_Kumar_R-0A66C2?style=for-the-badge&logo=LinkedIn&logoColor=white)](https://www.linkedin.com/in/sarathkr/)



