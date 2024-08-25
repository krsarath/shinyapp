# Shiny VISTA

# A Web-Application provides Differential Expression Analysis and Transcriptome data visualization 

This R Shiny web application allows users to perform differential expression analysis and visualize transcriptome data, making it easier to interpret gene expression patterns. The app supports uploading CSV/XLSX files, performing statistical analysis, and generating various plots to aid in the analysis.

- Authors : Sarath kumar R and Dr. J. Sreekumar.
- Website : https://krsarath.shinyapps.io/shinyvistaApp/ shinyvista Ver. 1.0

## Features

- **Differential Expression Analysis**: Identify differentially expressed genes between conditions.
- **Data Upload**: Upload transcriptome data in CSV or XLSX format.
- **Interactive Visualizations**:
  - Volcano plots
  - MA plots
  - Heatmaps
  - Venn diagrams
- **Downloadable Results**: Export analysis results and plots.

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
install.packages(c("shiny", "ggplot2", "DESeq2", "edgeR", "Noiseq", "dplyr", "shinyWidgets", "readxl", "readr"))

#  Running the App

To run the app locally, use the following R command
shiny::runApp('path_to_your_app_directory')



# How to Use

- Upload RNA-seq Count File: Use the 'Upload' tab to select your count file.
- Perform Differential Expression Analysis: Navigate to the 'Differential Expression' tab, choose your settings, and click 'Run'.
- Visualize GO Terms: In the 'GO Bubble Plot' tab, upload your DE results and GO files, then generate the plot.
- Create Venn Diagrams: Use the 'Venn Diagram' tab to upload multiple DE result files and compare them visually.
- View Results: Explore the generated plots and tables.
- Download: Export the results and plots for further analysis.


# Example Data

Example datasets are included in the data directory to help you get started.[matrix.counts (copy).zip](url)

# License

This project is licensed under the MIT License - see the LICENSE file for details.[LICENSE](url)

# Contact

For any questions or issues, please contact

[![Gmail](https://img.shields.io/badge/Sarath_Kumar_R-DB4437?style=for-the-badge&logo=Gmail&logoColor=white)](mailto:akhilsarath37@gmail.com) [![LinkedIn](https://img.shields.io/badge/Sarath_Kumar_R-0A66C2?style=for-the-badge&logo=LinkedIn&logoColor=white)](https://www.linkedin.com/in/sarathkr/)



