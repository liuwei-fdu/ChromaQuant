# ChromaQuant: Chromatographic Peak Quantification Software

[![Shiny App](https://img.shields.io/badge/launch-shiny%20app-blue.svg)](https://hplcfdu.shinyapps.io/ChromaQuant/)

ChromaQuant is a customized chromatographic data analysis platform developed using the R programming language and the Shiny framework. The software enables automated, batch-oriented processing of chromatographic data, significantly enhancing reproducibility and analytical throughput. It provides an interactive interface for data import, preprocessing, peak detection, quantification, statistical analysis, and advanced visualization.

## Overview & Workflow

The application guides the user from raw data upload to final analysis through a series of dedicated tabs. The workflow includes loading chromatogram and sample information, applying corrections, and generating various analytical plots and data tables.

![ChromaQuant Interface](https://raw.githubusercontent.com/liuwei-fdu/ChromaQuant/main/ChromaQuant%20Interface.png)

## Key Features

* **Batch Processing:** Directly import and process multiple raw chromatographic data files from a user-specified directory.
* **Automated Preprocessing:** Includes numerical conversion, formatting, baseline correction, and noise reduction.
* **Peak Detection & Quantification:**
    * Utilizes the `findpeaks` function from the `pracma` package for robust peak detection.
    * Employs custom algorithms to cluster closely eluting peaks.
    * Quantifies peak areas using spline interpolation and numerical integration.
* **Data Correction & Normalization:**
    * **Peak-area normalization** to correct for variations in sample loading or instrument drift.
    * **Retention time correction** via linear regression to ensure accurate peak alignment across samples.
* **Statistical Analysis:**
    * Performs **Analysis of Variance (ANOVA)** for rigorous inter-group comparisons of peak areas.
    * Calculates **Coefficient of Variation (CV)** to assess data reproducibility.
* **Advanced Visualization:**
    * **Dimensionality Reduction:** Principal Component Analysis (PCA), t-SNE, and UMAP plots for sample clustering.
    * **Data Quality Plots:** Boxplots for peak area distribution, violin plots for CV, and deviation plots for retention time stability.
    * **Correlation Plots:** Assess the correlation of peak intensities between samples within the same group.
* **Interactive Interface:** Dynamically adjust parameters, select data transformations (e.g., Log2, Log10), and view results in real-time.
* **Exporting Results:** All generated tables and plots can be easily downloaded.

## Live Demo

A live version of the application is hosted on shinyapps.io. You can access it here:
[**https://hplcfdu.shinyapps.io/ChromaQuant/**](https://hplcfdu.shinyapps.io/ChromaQuant/)

## Installation and Local Usage

To run ChromaQuant on your local machine, you will need to have R and RStudio installed.

### 1. Prerequisites

* [R](https://cran.r-project.org/) (version 4.0 or higher recommended)
* [RStudio](https://www.rstudio.com/products/rstudio/download/) (recommended for ease of use)

### 2. Dependencies

First, you need to install the required R packages. Open R or RStudio and run the following command in the console to install them all at once:

```R
install.packages(c(
  "shiny", "shinyFiles", "data.table", "dplyr", "stringr", "readr", 
  "pracma", "zoo", "ggplot2", "gridExtra", "ggsci", "umap", "splancs", 
  "shinyjqui", "rstatix", "shinyWidgets", "vioplot", "ggpubr", "tidyr"
))
```

### 3. Running the App

1.  Download or clone this repository.
2.  Save the `ChromaQuant.R` file (and any other repository files) to a directory on your computer.
3.  Open the `ChromaQuant.R` file in RStudio.
4.  Click the "Run App" button at the top of the script editor, or run the following command in the console:

```R
# Make sure your working directory is set to where the file is located
# setwd("path/to/your/folder")
shiny::runApp('ChromaQuant.R')
```

## Input Data Format

ChromaQuant requires two types of input: chromatogram data files and a sample information file.

* **Chromatogram Data:** Raw data should be in `.txt` format, with each file representing a single sample. The software is designed to parse text files containing specific headers (e.g., `[Header]`, `Wavelength(nm)`).
* **Sample Information:** A `.csv` file is required to map the sample file names to their respective experimental groups. This file should contain a `sample_name` column (without the `.txt` extension) and at least one column defining the groups (e.g., `group`).
* **Demo Data:** A demo dataset is available for download directly from the app's UI to help you get started.

