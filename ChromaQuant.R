# ==============================
# Load required packages
# ==============================
library(shiny)
library(shinyFiles)
library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(pracma)
library(zoo)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(umap)
library(splancs)
library(shinyjqui)
library(rstatix)
library(shinyWidgets)
library(vioplot)
library(ggpubr)
library(tidyr)

# ==============================
# Functions
# ==============================

# Function to read and merge raw data
peak_intensity_combined <- function(folder_path, wavelength) {
  file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)
  
  all_intensities <- list()
  
  for (file in file_list) {
    LC_raw <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>% as.data.frame()
    
    sample_name <- gsub(paste0(folder_path), "", file)
    sample_name <- gsub(".txt", "", sample_name)
    
    header_name <- paste0("Wavelength(nm)\t", wavelength, "nm")
    intensity_from <- which(LC_raw$`[Header]` == header_name) + 1
    
    # Locate by wavelength header
    if (wavelength == "280") {
      intensity_to <- which(LC_raw$`[Header]` == "[LC Chromatogram(Detector A-Ch2)]") - 1
    } else if (wavelength == "215") {
      intensity_to <- which(LC_raw$`[Header]` == "[LC Status Trace(Pump A Pressure)]") - 1
    } else {
      stop("Invalid wavelength. Please enter '215' or '280'.")
    }
    
    intensity_data <- LC_raw[intensity_from:intensity_to, ]
    data_split <- strsplit(intensity_data, "\t")
    data_matrix <- matrix(
      unlist(data_split),
      ncol = length(data_split[[1]]),
      byrow = TRUE
    )
    colnames(data_matrix) <- data_split[[1]]
    data_matrix <- data_matrix[-1, ]
    
    R_time <- as.numeric(data_matrix[, 1])
    intensity <- as.numeric(data_matrix[, 2])
    
    all_intensities[[sample_name]] <- list(R_time = R_time, intensity = intensity)
  }
  
  # Merge intensities into a matrix
  result_matrix <- data.frame(R_time = all_intensities[[1]]$R_time)
  for (sample_name in names(all_intensities)) {
    col_name <- sample_name
    result_matrix[[col_name]] <- all_intensities[[sample_name]]$intensity
  }
  
  return(result_matrix)
}

# Function to find peaks
find_peaks <- function(data2, tolerance = 0.5) {
  
  # Internal function: find initial peaks in each column
  find_and_filter_peaks <- function(data2) {
    peak_rt <- c()
    for (i in 2:ncol(data2)) {
      subset_data <- data2[[i]][data2[[1]] >= 25 & data2[[1]] <= 100]
      max_height <- max(subset_data)
      min_peak_height <- 0.1 * max_height
      intensity <- data2[[i]]
      
      peak_indices <- findpeaks(
        intensity,
        nups = 1,
        ndowns = 3,
        zero = "-",
        minpeakheight = min_peak_height
      )
      
      peak_rt0 <- peak_indices[,2] * 0.1
      peak_rt0 <- peak_rt0[peak_rt0 > 25 & peak_rt0 < 100]
      peak_rt <- c(peak_rt, peak_rt0)
    }
    return(peak_rt)
  }
  
  # Internal function: group close RTs
  group_retention_times <- function(sorted_data, tolerance) {
    groups <- list()
    current_group <- c(sorted_data[1])
    
    for (i in 2:length(sorted_data)) {
      if (sorted_data[i] - current_group[1] <= tolerance) {
        current_group <- c(current_group, sorted_data[i])
      } else {
        groups <- c(groups, list(current_group))
        current_group <- c(sorted_data[i])
      }
    }
    groups <- c(groups, list(current_group))
    return(groups)
  }
  
  # Internal function: merge group averages that are close
  merge_groups <- function(groups, tolerance) {
    result <- list()
    current_group <- groups[[1]]
    
    for (i in 2:length(groups)) {
      if (abs(mean(groups[[i]]) - mean(current_group)) <= tolerance) {
        current_group <- c(current_group, groups[[i]])
      } else {
        result <- c(result, list(current_group))
        current_group <- groups[[i]]
      }
    }
    result <- c(result, list(current_group))
    return(result)
  }
  
  peak_rt <- find_and_filter_peaks(data2)
  sorted_data <- sort(unique(peak_rt))
  groups <- group_retention_times(sorted_data, tolerance)
  merged_groups <- merge_groups(groups, tolerance)
  return(merged_groups)
}

# Function to process peak retention time matrix
process_peak_rt <- function(data2, groups) {
  result_rt <- data.frame()
  
  for (i in 2:ncol(data2)) {
    name <- colnames(data2)[i]
    intensity <- data2[[i]]
    subset_data <- data2[[i]][data2[[1]] >= 25 & data2[[1]] <= 100]
    max_height <- max(subset_data)
    min_peak_height <- 0.1 * max_height
    
    peak_indices <- findpeaks(intensity, nups = 1, ndowns = 3, zero="-", minpeakheight = min_peak_height)
    
    peak_area <- numeric(nrow(peak_indices))
    for (j in 1:nrow(peak_indices)) {
      x <- peak_indices[j, 3]:peak_indices[j, 4]
      y <- intensity[x]
      f <- splinefun(x, y, method = "natural")
      peak_area[j] <- integrate(f, lower = x[1], upper = x[length(x)])$value
    }
    peak_indices <- cbind(peak_indices, peak_area)
    
    peak_indices_sub <- peak_indices[peak_indices[, 2] > 250 & peak_indices[, 2] < 1000, ]
    peak_rt <- peak_indices_sub[, 2] * 0.1
    
    group_names <- sapply(peak_rt, function(x) {
      for (m in 1:length(groups)) {
        if (x >= min(groups[[m]]) && x <= max(groups[[m]])) {
          return(paste("peak", m))
        }
      }
    })
    
    result0 <- data.frame(Value = peak_rt, Group = group_names)
    names(result0) <- c(name, "peak")
    
    if (i == 2) {
      result_rt <- result0
    } else {
      result_rt <- merge(result_rt, result0, by = "peak", all = TRUE)
    }
  }
  
  result_rt <- result_rt %>%
    group_by(peak) %>%
    summarize(across(
      everything(),
      ~ if (all(is.na(.))) NA else max(., na.rm = TRUE)
    )) %>%
    as.data.frame()
  return(result_rt)
}

# Function to process peak area data (version 2), returns peak area table
process_peak_area2 <- function(data2, groups) {
  result <- data.frame()
  
  for (i in 2:ncol(data2)) {
    name <- colnames(data2)[i]
    intensity <- data2[[i]]
    
    subset_data <- data2[[i]][data2[[1]] >= 25 & data2[[1]] <= 100]
    max_height <- max(subset_data)
    min_peak_height <- 0.1 * max_height
    
    peak_indices <- findpeaks(intensity, nups = 1, ndowns = 3, zero="-", minpeakheight = min_peak_height)
    peak_indices <- peak_indices[peak_indices[, 2] > 250 & peak_indices[, 2] < 1000, ]
    
    peak_area <- numeric(nrow(peak_indices))
    
    for (j in 1:nrow(peak_indices)) {
      # Approximate area calculation (closed polygon) around peak ±5 points
      x <- (peak_indices[j, 2] - 5):(peak_indices[j, 2] + 5)
      y <- intensity[x]
      x <- c(peak_indices[j, 2] - 5, x, peak_indices[j, 2] + 5)
      y <- c(0, y, 0)
      
      peak_area[j] <- areapl(cbind(x, y))
    }
    
    peak_indices <- cbind(peak_indices, peak_area)
    peak_indices_sub <- peak_indices[peak_indices[, 2] > 250 & peak_indices[, 2] < 1000, ]
    peak_rt <- peak_indices_sub[, 2] * 0.1
    peak_area <- peak_indices_sub[, 5]
    
    # Match to peak groups
    group_names <- sapply(peak_rt, function(x) {
      for (m in 1:length(groups)) {
        if (x >= min(groups[[m]]) && x <= max(groups[[m]])) {
          return(paste("peak", m))
        }
      }
      return("No Group")
    })
    
    result0 <- data.frame(Value = peak_rt, Peak_Area = peak_area, peak = group_names)
    names(result0) <- c(name, "Peak_Area", "peak")
    
    if (i == 2) {
      result <- result0[, c("Peak_Area", "peak")]
      names(result) <- c(name, "peak")
    } else {
      result0 <- result0[, c("Peak_Area", "peak")]
      names(result0) <- c(name, "peak")
      result <- merge(result, result0, by = "peak", all = TRUE)
    }
  }
  
  result <- result %>%
    group_by(peak) %>%
    summarize(across(
      everything(),
      ~ if (all(is.na(.))) NA else max(., na.rm = TRUE)
    )) %>%
    as.data.frame()
  
  return(result)
}

# Function to calculate ANOVA
calculate_anova <- function(data) {
  # Return "no-pvalue" if not enough replicates in any group
  if (length(unique(data$group)) < 2 || any(table(data$group) < 2)) {
    return(tibble(p = "no-pvalue"))
  }
  
  tryCatch({
    anova_result <- anova_test(data, value ~ group)
    return(tibble(p = as.character(anova_result$p)))  # Convert p-value to string
  }, error = function(e) {
    return(tibble(p = "no-pvalue"))
  })
}

# Function to draw t-SNE
drawTSNE <- function(DF, ptColors, rowNormalization=F, colNormalization=F,
                     perplexity=10, strTitle='tSNE') {
  M <- DF[, colnames(DF) != 'label']
  if(rowNormalization){
    M <- data.frame(t(apply(M, 1, function(v){
      (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
    })))
  }
  if(colNormalization){
    M <- apply(M, 2, function(v){
      (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
    })
  }
  M[is.na(M)] <- 0
  
  tsn = tsne(M, perplexity=perplexity)
  tsn <- data.frame(tsn, DF$label)
  colnames(tsn) <- c("X", "Y", "label")
  rownames(tsn) <- rownames(M)
  
  p <- ggplot(tsn, aes(x=X, y=Y, colour=label)) + geom_point(size=4)
  p <- p + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.title   = element_text(size=15),
    axis.line.x  = element_line(color="black", size = 0.5),
    axis.line.y  = element_line(color="black", size = 0.5),
    panel.background = element_blank()
  )
  p <- p + labs(title=strTitle)
  p <- p + scale_colour_manual(values=ptColors)
  p
}

# Function to draw PCA
drawPCA <- function(DF, ptColors, rowNormalization=F, colNormalization=F, strTitle=NULL){
  M <- DF[, colnames(DF) != 'label']
  
  if(rowNormalization){
    M <- data.frame(t(apply(M, 1, function(v){
      (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
    })))
  }
  if(colNormalization){
    M <- apply(M, 2, function(v){
      (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
    })
  }
  M[is.na(M)] <- 0
  
  m1 <- prcomp(M, center=TRUE, scale.=FALSE)
  Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation
  Y  <- Y[, c(1, 2)]
  
  Y <- data.frame(Y, DF$label)
  colnames(Y) <- c("PC1", "PC2", "label")
  
  if(is.null(strTitle)){
    strTitle <- sprintf("PCA: %d features", ncol(M))
  }
  
  eigs <- m1$sdev^2
  percentages <- eigs[1:2] / sum(eigs)
  
  p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=4)
  p <- p + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x  = element_line(color="black", size = 0.25),
    axis.line.y  = element_line(color="black", size = 0.25),
    plot.title   = element_text(size=16),
    panel.background = element_blank()
  )
  
  strLabx <- sprintf("PC1 (%.2f%%)", percentages[1]*100)
  p <- p + labs(
    x = strLabx,
    y = sprintf("PC2 (%.2f%%)", percentages[2]*100),
    title = strTitle
  )
  p <- p + scale_colour_manual(values=ptColors)
  p
}

# Function to draw UMAP
drawUMAP <- function(M1, ptColors, strTitle="UMAP", rowNormalization=TRUE, colNormalization=FALSE) {
  if(! 'label' %in% colnames(M1)){
    stop("A column named 'label' must exist in the data frame.")
  }
  
  tmp <- M1[, colnames(M1) != 'label']
  if(rowNormalization){
    tmp <- data.frame(t(apply(tmp, 1, function(v){
      (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
    })), stringsAsFactors=FALSE)
    rownames(tmp) <- rownames(M1)
  }
  if(colNormalization){
    tmp <- apply(tmp, 2, function(v){
      (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
    })
  }
  tmp[is.na(tmp)] <- 0
  
  set.seed(2020)
  obj <- umap(d=tmp, method='naive')
  df1 <- data.frame(obj$layout)
  df1$label <- M1$label
  colnames(df1) <- c('X', 'Y', 'label')
  
  p <- ggplot(df1, aes(x=X, y=Y, colour=label)) + geom_point(size=4)
  p <- p + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.title   = element_text(size=16),
    axis.line.x  = element_line(color="black", size = 0.5),
    axis.line.y  = element_line(color="black", size = 0.5),
    panel.background = element_blank()
  )
  p <- p + labs(title=strTitle)
  p <- p + scale_colour_manual(values=ptColors)
  p
}

# ==============================
# UI
# ==============================
ui <- fluidPage(
  titlePanel("ChromaQuant: Chromatogram Data Analysis"),
  sidebarLayout(
    sidebarPanel(
      tags$a(
        href = "demo_data.zip",          
        "Download demo data",          
        download = NA,                 
        target = "_blank",               
        class = "btn btn-primary"        
      ),
      
      fileInput("files", "Upload Chromatogram Data", multiple = TRUE, accept = ".txt"),
      fileInput("file1", "Upload Sample Info",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      selectInput("group_column", "Select Group Column:", choices = NULL),
      uiOutput("level_order_ui"),
      selectInput("data_transform", "Data Transformation:",
                  choices = list("None" = "none", "Log2" = "log2", "Log10" = "log10"), selected = "log2"),
      numericInput("wavelength", "Wavelength (nm):", value = 280),
      numericInput("peak_elution_time", "Peak Elution Time (min):", value = 25, min = 0, max = 100),
      numericInput("ylimMax", "Y-axis Max Value:", value = 22000, min = 1000, max = 1000000),
      numericInput("tolerance", "Peak Tolerance:", value = 0.5, min = 0, max = 2),
      materialSwitch(inputId = "peak_area_correction", label = "Peak Area Correction", status = "primary", right = TRUE),
      materialSwitch(inputId = "retention_time_correction", label = "Retention Time Correction", status = "primary", right = TRUE),
      downloadButton("downloadData", "Download Processed Data")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Peak Retention Time Matrix", tableOutput("peakRtsTable"),
                 downloadButton("downloadPeakRts", "Download Peak Retention Time Matrix")),
        tabPanel("Peak Area Preview", dataTableOutput("peakAreaTable"),
                 downloadButton("downloadPeakArea", "Download Peak Area Data")),
        tabPanel("Chromatogram", plotOutput("chromatographicPlots", height = "8000px"),
                 downloadButton("downloadChromatographicPlots", "Download Chromatogram")),
        tabPanel("Peak Area", plotOutput("peakAreaPlots"),
                 downloadButton("downloadPeakAreaPlots", "Download Peak Area Plot")),
        tabPanel("Deviation Boxplot", plotOutput("deviationsBoxplot"),
                 downloadButton("downloadDeviationsBoxplot", "Download Deviation Boxplot")),
        tabPanel("UMAP", plotOutput("umapPlot"),
                 downloadButton("downloadUmapPlot", "Download UMAP Plot")),
        tabPanel("PCA", plotOutput("pcaPlot"),
                 downloadButton("downloadPcaPlot", "Download PCA Plot")),
        tabPanel("t-SNE", plotOutput("tsnePlot"),
                 downloadButton("downloadTsnePlot", "Download t-SNE Plot")),
        tabPanel("CV Violin Plot", plotOutput("cvViolinPlot"),
                 downloadButton("downloadCvViolinPlot", "Download CV Violin Plot")),
        tabPanel("Correlation Plot",
                 plotOutput("correlationMatrixPlots", width = "100%", height = 2000),
                 downloadButton("downloadCorrelationMatrixPlots", "Download Correlation Plots"))
      )
    )
  )
)

# ==============================
# Server
# ==============================
server <- function(input, output, session) {
  shinyDirChoose(input, 'folder', roots = c(home = '~', root = '/'), session = session)
  
  folder_path <- reactive({
    req(input$folder)
    parseDirPath(c(home = '~', root = '/'), input$folder)
  })
  
  data_input <- reactive({
    req(input$files)
    file_paths <- input$files$datapath
    file_names <- input$files$name
    
    all_intensities <- list()
    for(i in seq_along(file_paths)){
      LC_raw <- read_delim(file_paths[i], delim="\t", escape_double=FALSE, trim_ws=TRUE, show_col_types=FALSE) %>% as.data.frame()
      

     
      sample_name <- gsub(".txt", "", file_names[i])
      header_name <- paste0("Wavelength(nm)\t", input$wavelength)
      
      
      intensity_from <- which(grepl(header_name, LC_raw$`[Header]`, fixed=TRUE)) + 1
  

      vec <- LC_raw[[1]]
           is_num_line <- grepl("^\\s*-?\\d+\\.?\\d*\\s+-?\\d+\\.?\\d*", vec[intensity_from:length(vec)])
      first_non_num_index <- which(!is_num_line)[2]
      intensity_to <- intensity_from + first_non_num_index - 2
      
  
      
    
      intensity_data <- LC_raw[intensity_from:intensity_to, ]
      data_split <- strsplit(intensity_data, "\t")
      data_matrix <- matrix(unlist(data_split), ncol = length(data_split[[1]]), byrow = TRUE)
      colnames(data_matrix) <- data_split[[1]]
      data_matrix <- data_matrix[-1, ]
      R_time <- as.numeric(data_matrix[, 1])
      intensity <- as.numeric(data_matrix[, 2])
      
      all_intensities[[sample_name]] <- list(R_time = R_time, intensity = intensity)
    }
    
    result_matrix <- data.frame(R_time = all_intensities[[1]]$R_time)
    for (sample_name in names(all_intensities)) {
      result_matrix[[sample_name]] <- all_intensities[[sample_name]]$intensity
    }
    
    return(result_matrix)
  })
  
  # Data preprocessing function
  preprocessData <- function(data) {
    data <- apply(data, 2, as.numeric) %>% as.data.frame()
    data$RT <- ceiling(data$R_time * 10) / 10
    
    data2 <- data[-1, -1] %>%
      group_by(RT) %>%
      summarise(across(everything(), mean)) %>%
      replace(., . < 0, 0)
    return(data2)
  }
  
  # Peak area correction function
  correctPeakArea <- function(data) {
    if (input$peak_area_correction) {
      peak.elution.time <- input$peak_elution_time
      data.tmp1 <- filter(data, RT <= peak.elution.time)
      data.tmp2 <- filter(data, RT > peak.elution.time)
      
      data.tmp3 <- data.tmp2[, -1]
      meanall <- mean(apply(data.tmp3, 2, mean))
      df_normalized <- apply(data.tmp3, 2, function(x) x / mean(x) * meanall)
      
      data.tmp4 <- cbind(data.tmp2[1], df_normalized) %>% as.data.frame()
      names(data.tmp4)[1] <- "RT"
      data3 <- rbind(data.tmp1, data.tmp4)
    } else {
      data3 <- data
    }
    return(data3)
  }
  
  # Retention time correction function
  correctRetentionTime <- function(data) {
    if (input$retention_time_correction) {
      peaks <- find_peaks(data, tolerance = input$tolerance)
      result_rt <- process_peak_rt(data, peaks)
      
      peak.RT.cal <- result_rt$peak[apply(
        result_rt, 1,
        function(x) any(sum(!is.na(x)) > 0.9 * ncol(result_rt))
      )]
      
      result_rt.cal <- result_rt[result_rt$peak %in% peak.RT.cal, ]
      result_rt.cal1 <- result_rt.cal[, 2:ncol(result_rt.cal)]
      
      result_rt.cal1$exp.rt <- apply(result_rt.cal1, 1, function(x) mean(x, na.rm = TRUE))
      result_rt.cal1$exp.rt <- round(result_rt.cal1$exp.rt, 1)
      
      exp_rt <- result_rt.cal1$exp.rt
      columns_to_analyze <- setdiff(names(result_rt.cal1), "exp.rt")
      
      col_name <- columns_to_analyze[1]
      model <- lm(result_rt.cal1[[col_name]] ~ exp_rt)
      coefficients <- coef(model)
      intercept <- round(coefficients[1], 2)
      slope <- round(coefficients[2], 2)
      
      data.cal <- data[, names(data) %in% c("RT", col_name)]
      data.cal$RT.cal <- round((data.cal$RT - intercept) / slope, 1)
      data.cal <- data.cal[, -1]
      
      for (col_name in columns_to_analyze[-1]) {
        model <- lm(result_rt.cal1[[col_name]] ~ exp_rt)
        coefficients <- coef(model)
        intercept <- round(coefficients[1], 2)
        slope <- round(coefficients[2], 2)
        
        data.cal0 <- data[, names(data) %in% c("RT", col_name)]
        data.cal0$RT.cal <- round((data.cal0$RT - intercept) / slope, 1)
        data.cal0 <- data.cal0[, -1]
        data.cal <- merge(data.cal, data.cal0, by = "RT.cal", all = TRUE)
      }
      
      data.cal.aggregated <- aggregate(data.cal[, -1],
                                       by = list(data.cal$RT.cal), mean, na.rm = TRUE
      )
      names(data.cal.aggregated)[1] <- "RT"
      
      data4 <- data.cal.aggregated[data.cal.aggregated$RT > 0, ]
      data4[is.na(data4)] <- NA
      data4 <- apply(data4, 2, function(x) na.approx(x, rule = 2))
      data4 <- as.data.frame(data4)
    } else {
      data4 <- data
      data4 <- as.data.frame(data4)
    }
    return(data4)
  }
  
  # Read raw data from files
  # processed_data <- reactive({
  #   data2 <- preprocessData(data_input())
  #   data3 <- correctPeakArea(data2)
  #   data4 <- correctRetentionTime(data3)
  #   return(data4)
  # })
  # 
  processed_data <- reactive({
    req(data_input())
    input$peak_elution_time
    input$peak_area_correction
    data2 <- preprocessData(data_input())
    data3 <- correctPeakArea(data2)
    data4 <- correctRetentionTime(data3)
    data4
  })
  
  
  # Calculate peak groups
  peaks <- reactive({
    find_peaks(processed_data(), tolerance = input$tolerance)
  })
  
  # Calculate peak areas
  peak_areas <- reactive({
    peak_areas_raw <- process_peak_area2(processed_data(), peaks())
    transformed_data <- switch(
      input$data_transform,
      "none"  = peak_areas_raw[, -1],
      "log2"  = log2(peak_areas_raw[, -1] + 1),
      "log10" = log10(peak_areas_raw[, -1] + 1),
      peak_areas_raw[, -1]
    )
    peak_areas_transformed <- cbind(peak_areas_raw[, 1], transformed_data)
    colnames(peak_areas_transformed)[1] <- colnames(peak_areas_raw)[1]
    return(peak_areas_transformed)
  })
  
  # Calculate peak retention time matrix
  peak_rts <- reactive({
    process_peak_rt(processed_data(), peaks())
  })
  
  # Calculate deviations from peak means
  deviations <- reactive({
    lapply(peaks(), function(group) {
      mean_val <- mean(group)
      group - mean_val
    })
  })
  
  # Process intensity matrix for plots (boxplot, UMAP, PCA, etc.)
  processed_intensity_mat <- reactive({
    req(input$file1)
    intensity_mat <- peak_areas()[, 2:ncol(peak_areas())]
    row.names(intensity_mat) <- peak_areas()[, 1]
    
    intensity_mat1 <- data.frame(t(intensity_mat))
    sample_info <- read.csv(input$file1$datapath)
    intensity_mat1$sample_name <- row.names(intensity_mat1)
    
    intensity_mat1 <- merge(sample_info, intensity_mat1, by = "sample_name", all.x = TRUE)
    intensity_mat1[[input$group_column]] <- factor(intensity_mat1[[input$group_column]], levels = input$level_order)
    
    return(intensity_mat1)
  })
  
  # Calculate CV data for violin plot
  cv_data <- reactive({
    req(input$file1)
    
    peak_data_long <- processed_intensity_mat() %>%
      select(starts_with("peak"), input$group_column) %>%
      pivot_longer(cols = starts_with("peak"), names_to = "peak", values_to = "area")
    
    cv_data <- peak_data_long %>%
      group_by(!!sym(input$group_column), peak) %>%
      summarise(cv = sd(area, na.rm = TRUE) / mean(area, na.rm = TRUE) * 100, .groups = "drop") %>%
      group_by(!!sym(input$group_column)) %>%
      summarise(cv_values = list(cv), median_cv = median(cv, na.rm = TRUE))
    
    return(cv_data)
  })
  
  # Update group column and level order when sample_info.csv is uploaded
  observe({
    req(input$file1)
    sample_info <- read.csv(input$file1$datapath)
    
    group_col <- names(sample_info)[grepl("group", names(sample_info), ignore.case = TRUE)]
    if (length(group_col) > 0) {
      selected_col <- group_col[1]
    } else {
      selected_col <- names(sample_info)[1]
    }
    updateSelectInput(inputId = "group_column", choices = names(sample_info), selected = selected_col)
    
    output$level_order_ui <- renderUI({
      group_levels <- unique(sample_info[[input$group_column]])
      orderInput(inputId = "level_order", label = "Group Levels Order:", items = group_levels)
    })
  })
  
  # Render CV violin plot
  output$cvViolinPlot <- renderPlot({
    req(cv_data())
    
    plot(0, 0, type = "n",
         xlim = c(0.5, nrow(cv_data()) + 0.5),
         ylim = c(0, max(unlist(cv_data()$cv_values), na.rm = TRUE) * 1.1),
         xaxt = "n", xlab = "Group", ylab = "CV (%)", main = "Coefficient of Variation of Peak Area by Group")
    axis(1, at = 1:nrow(cv_data()), labels = cv_data()[[input$group_column]])
    
    for (i in 1:nrow(cv_data())) {
      vioplot(cv_data()$cv_values[[i]], at = i, add = TRUE, col = pal_npg("nrc")(7)[i])
      text(i, cv_data()$median_cv[i], labels = round(cv_data()$median_cv[i], 2), pos = 3, cex = 0.8)
    }
  })
  
  # Render correlation matrix plots
  output$correlationMatrixPlots <- renderPlot({
    req(input$file1, input$group_column)
    intensity_mat1 <- processed_intensity_mat()
    groups <- unique(intensity_mat1[[input$group_column]])
    
    all_plots <- list()
    
    for(group_val in groups) {
      group_data <- intensity_mat1 %>%
        filter(!!sym(input$group_column) == group_val) %>%
        select(starts_with("peak"))
      
      if(ncol(group_data) > 0) {
        samples <- rownames(group_data)
        if(length(samples) > 6) {
          samples <- samples[1:6]
        }
        
        combinations <- combn(samples, 2, simplify = FALSE)
        group_plots <- lapply(combinations, function(comb) {
          sample1 <- group_data[comb[1], ]
          sample2 <- group_data[comb[2], ]
          
          plot_data <- data.frame(
            peak = names(sample1),
            intensity1 = as.numeric(sample1),
            intensity2 = as.numeric(sample2)
          )
          
          ggplot(plot_data, aes(x = intensity1, y = intensity2)) +
            geom_point() +
            geom_smooth(method = "lm", se = FALSE) +
            stat_cor(method = "pearson",  label.x = min(plot_data$intensity1, na.rm = TRUE),
                     label.y = max(plot_data$intensity2, na.rm = TRUE)) +
            labs(
              title = paste("Group:", group_val, "-", comb[1], "VS", comb[2]),
              x = paste("Intensity of", comb[1]),
              y = paste("Intensity of", comb[2])
            ) +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_x_continuous(limits = range(plot_data$intensity1, na.rm = TRUE)) +   # X axis auto-range
            scale_y_continuous(limits = range(plot_data$intensity2, na.rm = TRUE))     # Y axis auto-range
        })
        
        if (length(group_plots) > 0) {
          num_plots <- length(group_plots)
          ncol_dynamic <- min(6, num_plots)
          nrow_dynamic <- ceiling(num_plots / ncol_dynamic)
          all_plots[[group_val]] <- do.call(grid.arrange, c(group_plots, ncol = ncol_dynamic, nrow = nrow_dynamic))
        }
      }
    }
    all_plots <- all_plots[!sapply(all_plots, is.null)]
    if(length(all_plots) > 0) {
      do.call(grid.arrange, c(all_plots, nrow = length(all_plots)))
    }
  })
  
  # Function to draw chromatogram plots
  drawChromatographicPlots <- function(data4) {
    par(mfrow = c(ncol(data4)-1, 1))
    for (i in 2:ncol(data4)) {
      title <- colnames(data4)[i]
      plot(data4[, 1], data4[, title],
           xlab = "R.Time (min)",
           ylab = "Intensity (µV)",
           type = "l",
           ylim = c(0, input$ylimMax),
           main = title)
      
      peak_rt_select <- as.numeric(c(as.matrix(peak_rts()[, title])))
      label <- peak_rts()$peak
      label <- label[!is.na(peak_rt_select)]
      label <- gsub("peak ", "", label)
      label <- as.numeric(label)
      label <- label[order(label)]
      
      peak_rt_select <- peak_rt_select[!is.na(peak_rt_select)]
      peak_rt_select <- round(peak_rt_select, 1)
      peak_rt_select <- peak_rt_select[order(peak_rt_select)]
      
      points(x = peak_rt_select, y = data4[data4[, 1] %in% c(peak_rt_select), title], col = label, pch = 1)
      text(x = peak_rt_select, y = data4[data4[, 1] %in% c(peak_rt_select), title], labels = label, pos = 3, cex = 0.8)
      legend("topright", legend = unique(label), col = unique(label), pch = 1)
    }
  }
  
  # Output peak retention time matrix table
  output$peakRtsTable <- renderTable({
    peak_rts()
  })
  
  # Output peak area data table
  output$peakAreaTable <- renderDataTable({
    peak_areas()
  })
  
  # Output chromatogram plots
  output$chromatographicPlots <- renderPlot({
    drawChromatographicPlots(processed_data())
  }, height = 8000)
  
  # Function to generate peak area boxplots
  generatePeakAreaPlots <- function() {
    req(input$file1)
    intensity_mat1 <- processed_intensity_mat()
    
    anova_results <- intensity_mat1 %>%
      pivot_longer(cols = starts_with("peak"), names_to = "variable", values_to = "value") %>%
      group_by(variable) %>%
      group_modify(~ calculate_anova(.)) %>%
      ungroup()
    
    df_long <- intensity_mat1 %>%
      pivot_longer(cols = starts_with("peak"), names_to = "variable", values_to = "value") %>%
      group_by(!!sym(input$group_column), variable) %>%
      summarise(value = value, .groups = "drop")
    
    num_plots <- length(unique(df_long$variable))
    plot_list <- list()
    
    for (i in 1:num_plots) {
      df_subset <- df_long[df_long$variable == unique(df_long$variable)[i], ]
      peak_name <- gsub("_mean$", "", unique(df_long$variable)[i])
      
      p_value <- anova_results$p[anova_results$variable == peak_name]
      title <- paste0(peak_name, "\n")
      if (!is.na(p_value) && p_value != "no-pvalue") {
        title <- paste0(title, "ANOVA p = ", format.pval(as.numeric(p_value), digits = 3))
      } else {
        title <- paste0(title, "Insufficient data for ANOVA")
      }
      
      plot_list[[i]] <- ggplot(df_subset, aes(x = !!sym(input$group_column), y = value, fill = !!sym(input$group_column))) +
        geom_boxplot() +
        ggtitle(title) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.title.x = element_text(margin = margin(t = 15)),
          axis.title.y = element_text(margin = margin(r = 15)),
          axis.text = element_text(size = 10),
          axis.ticks = element_line(color = "black"),
          axis.ticks.length = unit(0.3, "cm")
        ) +
        scale_fill_npg()
    }
    return(plot_list)
  }
  
  # Output peak area boxplots
  output$peakAreaPlots <- renderPlot({
    plot_list <- generatePeakAreaPlots()
    num_plots <- length(plot_list)
    num_columns <- 3
    num_rows <- ceiling(num_plots / num_columns)
    
    grid.arrange(grobs = plot_list, nrow = num_rows, ncol = num_columns)
  }, height = 5000, width = 1200)
  
  # Output deviations boxplot
  output$deviationsBoxplot <- renderPlot({
    boxplot(deviations(),
            main = "Deviations from Peaks Means",
            ylab = "Deviations of Peaks RT (min)",
            xlab = "Peaks",
            outline = TRUE)
    abline(h = 0, col = "red", lty = 2)
  })
  
  # Output UMAP plot
  output$umapPlot <- renderPlot({
    intensity_mat1 <- processed_intensity_mat()
    anova_results <- intensity_mat1 %>%
      pivot_longer(cols = starts_with("peak"), names_to = "variable", values_to = "value") %>%
      group_by(variable) %>%
      group_modify(~ calculate_anova(.)) %>%
      ungroup()
    
    anova_results$p <- as.numeric(anova_results$p)
    select_peaks <- c(as.matrix(anova_results[anova_results$p < 0.05, ]$variable))
    select_peaks <- select_peaks[!is.na(select_peaks)]
    
    pca.mat <- intensity_mat1[, c(select_peaks, "group")]
    names(pca.mat) <- gsub("group", "label", names(pca.mat))
    npg_colors <- ggsci::pal_npg("nrc")(7)
    
    drawUMAP(pca.mat, rowNormalization = FALSE, colNormalization = TRUE, ptColors = npg_colors)
  })
  
  # Output PCA plot
  output$pcaPlot <- renderPlot({
    intensity_mat1 <- processed_intensity_mat()
    anova_results <- intensity_mat1 %>%
      pivot_longer(cols = starts_with("peak"), names_to = "variable", values_to = "value") %>%
      group_by(variable) %>%
      group_modify(~ calculate_anova(.)) %>%
      ungroup()
    
    anova_results$p <- as.numeric(anova_results$p)
    select_peaks <- c(as.matrix(anova_results[anova_results$p < 0.05, ]$variable))
    select_peaks <- select_peaks[!is.na(select_peaks)]
    
    pca.mat <- intensity_mat1[, c(select_peaks, "group")]
    names(pca.mat) <- gsub("group", "label", names(pca.mat))
    npg_colors <- ggsci::pal_npg("nrc")(7)
    
    drawPCA(pca.mat, rowNormalization = FALSE, colNormalization = TRUE, ptColors = npg_colors)
  })
  
  # Output t-SNE plot
  output$tsnePlot <- renderPlot({
    intensity_mat1 <- processed_intensity_mat()
    anova_results <- intensity_mat1 %>%
      pivot_longer(cols = starts_with("peak"), names_to = "variable", values_to = "value") %>%
      group_by(variable) %>%
      group_modify(~ calculate_anova(.)) %>%
      ungroup()
    
    anova_results$p <- as.numeric(anova_results$p)
    select_peaks <- c(as.matrix(anova_results[anova_results$p < 0.05, ]$variable))
    select_peaks <- select_peaks[!is.na(select_peaks)]
    
    pca.mat <- intensity_mat1[, c(select_peaks, "group")]
    names(pca.mat) <- gsub("group", "label", names(pca.mat))
    npg_colors <- ggsci::pal_npg("nrc")(7)
    
    drawTSNE(pca.mat, rowNormalization = FALSE, colNormalization = TRUE, ptColors = npg_colors)
  })
  
  # Download processed data
  output$downloadData <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("processed_data_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(processed_data(), file)
    }
  )
  
  # Download peak retention times
  output$downloadPeakRts <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("peak_retention_times_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(peak_rts(), file)
    }
  )
  
  # Download peak area data
  output$downloadPeakArea <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      data_transform_str <- input$data_transform
      paste0("peak_area_data_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", data_transform_str, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(peak_areas(), file)
    }
  )
  
  # Download chromatogram plots
  output$downloadChromatographicPlots <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("chromatographic_plots_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 15, height = 200)
      drawChromatographicPlots(processed_data())
      dev.off()
    }
  )
  
  # Download peak area boxplots
  output$downloadPeakAreaPlots <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("peak_area_plots_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      plot_list <- generatePeakAreaPlots()
      num_plots <- length(plot_list)
      num_columns <- 3
      num_rows <- ceiling(num_plots / num_columns)
      pdf(file, width = 15, height = 50)
      grid.arrange(grobs = plot_list, nrow = num_rows, ncol = num_columns)
      dev.off()
    }
  )
  
  # Download deviations boxplot
  output$downloadDeviationsBoxplot <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("deviations_boxplot_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      boxplot(deviations(),
              main = "Deviations from Peaks Means",
              ylab = "Deviations of Peaks RT (min)",
              xlab = "Peaks",
              outline = TRUE)
      abline(h = 0, col = "red", lty = 2)
      dev.off()
    }
  )
  
  # Download UMAP plot
  output$downloadUmapPlot <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("umap_plot_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 10, height = 8)
      intensity_mat1 <- processed_intensity_mat()
      anova_results <- intensity_mat1 %>%
        pivot_longer(cols = starts_with("peak"), names_to = "variable", values_to = "value") %>%
        group_by(variable) %>%
        group_modify(~ calculate_anova(.)) %>%
        ungroup()
      anova_results$p <- as.numeric(anova_results$p)
      select_peaks <- c(as.matrix(anova_results[anova_results$p < 0.05, ]$variable))
      select_peaks <- select_peaks[!is.na(select_peaks)]
      
      pca.mat <- intensity_mat1[, c(select_peaks, "group")]
      names(pca.mat) <- gsub("group", "label", names(pca.mat))
      npg_colors <- ggsci::pal_npg("nrc")(7)
      print(drawUMAP(pca.mat, rowNormalization = FALSE, colNormalization = TRUE, ptColors = npg_colors))
      dev.off()
    }
  )
  
  # Download PCA plot
  output$downloadPcaPlot <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("pca_plot_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 10, height = 8)
      intensity_mat1 <- processed_intensity_mat()
      anova_results <- intensity_mat1 %>%
        pivot_longer(cols = starts_with("peak"), names_to = "variable", values_to = "value") %>%
        group_by(variable) %>%
        group_modify(~ calculate_anova(.)) %>%
        ungroup()
      anova_results$p <- as.numeric(anova_results$p)
      select_peaks <- c(as.matrix(anova_results[anova_results$p < 0.05, ]$variable))
      select_peaks <- select_peaks[!is.na(select_peaks)]
      
      pca.mat <- intensity_mat1[, c(select_peaks, "group")]
      names(pca.mat) <- gsub("group", "label", names(pca.mat))
      npg_colors <- ggsci::pal_npg("nrc")(7)
      print(drawPCA(pca.mat, rowNormalization = FALSE, colNormalization = TRUE, ptColors = npg_colors))
      dev.off()
    }
  )
  
  # Download t-SNE plot
  output$downloadTsnePlot <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("tsne_plot_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 10, height = 8)
      intensity_mat1 <- processed_intensity_mat()
      anova_results <- intensity_mat1 %>%
        pivot_longer(cols = starts_with("peak"), names_to = "variable", values_to = "value") %>%
        group_by(variable) %>%
        group_modify(~ calculate_anova(.)) %>%
        ungroup()
      anova_results$p <- as.numeric(anova_results$p)
      select_peaks <- c(as.matrix(anova_results[anova_results$p < 0.05, ]$variable))
      select_peaks <- select_peaks[!is.na(select_peaks)]
      pca.mat <- intensity_mat1[, c(select_peaks, "group")]
      names(pca.mat) <- gsub("group", "label", names(pca.mat))
      npg_colors <- ggsci::pal_npg("nrc")(7)
      print(drawTSNE(pca.mat, rowNormalization = FALSE, colNormalization = TRUE, ptColors = npg_colors))
      dev.off()
    }
  )
  
  # Download CV violin plot
  output$downloadCvViolinPlot <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      data_transform_str <- input$data_transform
      paste("cv_violin_plot_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", data_transform_str, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 12, height = 8)
      req(cv_data())
      
      plot(0, 0, type = "n",
           xlim = c(0.5, nrow(cv_data()) + 0.5),
           ylim = c(0, max(unlist(cv_data()$cv_values), na.rm = TRUE) * 1.1),
           xaxt = "n", xlab = "Group", ylab = "CV (%)",
           main = "Coefficient of Variation of Peak Area by Group")
      axis(1, at = 1:nrow(cv_data()), labels = cv_data()[[input$group_column]])
      
      for (i in 1:nrow(cv_data())) {
        vioplot(cv_data()$cv_values[[i]], at = i, add = TRUE, col = pal_npg("nrc")(7)[i])
        text(i, cv_data()$median_cv[i], labels = round(cv_data()$median_cv[i], 2), pos = 3, cex = 0.8)
      }
      
      dev.off()
    }
  )
  
  # Download correlation matrix plots
  output$downloadCorrelationMatrixPlots <- downloadHandler(
    filename = function() {
      peak_area_str <- ifelse(input$peak_area_correction, "T", "F")
      retention_time_str <- ifelse(input$retention_time_correction, "T", "F")
      paste("correlation_matrix_plots_PeakareaCal_", peak_area_str, "_RtCal_", retention_time_str, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      groups <- unique(processed_intensity_mat()[[input$group_column]])
      # For each group, count how many plots will be generated
      plots_per_group <- lapply(groups, function(group_val) {
        group_data <- processed_intensity_mat() %>%
          filter(!!sym(input$group_column) == group_val) %>%
          select(starts_with("peak"))
        samples <- rownames(group_data)
        if (length(samples) > 6) samples <- samples[1:6]
        n_comb <- choose(length(samples), 2)
        return(n_comb)
      })
      n_total_plots <- sum(unlist(plots_per_group))
      plots_per_row <- 3
      nrow_total <- ceiling(n_total_plots / plots_per_row)
      pdf_width <- plots_per_row * 5
      pdf_height <- max(5, nrow_total * 5)
      pdf(file, width = pdf_width, height = pdf_height)
      for (group_val in groups) {
        group_data <- processed_intensity_mat() %>%
          filter(!!sym(input$group_column) == group_val) %>%
          select(starts_with("peak"))
        samples <- rownames(group_data)
        if (length(samples) > 6) samples <- samples[1:6]
        combinations <- combn(samples, 2, simplify = FALSE)
        group_plots <- lapply(combinations, function(comb) {
          sample1 <- group_data[comb[1], ]
          sample2 <- group_data[comb[2], ]
          plot_data <- data.frame(
            peak = names(sample1),
            intensity1 = as.numeric(sample1),
            intensity2 = as.numeric(sample2)
          )
          ggplot(plot_data, aes(x = intensity1, y = intensity2)) +
            geom_point() +
            geom_smooth(method = "lm", se = FALSE) +
            stat_cor(method = "pearson",   label.x = min(plot_data$intensity1, na.rm = TRUE),
                     label.y = max(plot_data$intensity2, na.rm = TRUE)) +
            labs(title = paste("Group:", group_val, "-", comb[1], "VS", comb[2]),
                 x = paste("Intensity of", comb[1]),
                 y = paste("Intensity of", comb[2])) +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_x_continuous(limits = range(plot_data$intensity1, na.rm = TRUE)) +
            scale_y_continuous(limits = range(plot_data$intensity2, na.rm = TRUE))
        })
        if (length(group_plots) > 0) {
          num_plots <- length(group_plots)
          ncol_dynamic <- min(plots_per_row, num_plots)
          nrow_dynamic <- ceiling(num_plots / ncol_dynamic)
          do.call(grid.arrange, c(group_plots, ncol = ncol_dynamic, nrow = nrow_dynamic))
        }
      }
      dev.off()
    }
  )
}

# ==============================
# Run the app
# ==============================
shinyApp(ui = ui, server = server)