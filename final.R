# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(readr)

# UPLOAD SIZE LIMIT to 100MB
options(shiny.maxRequestSize = 100 * 1024^2)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "AMR Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Data Overview", tabName = "overview", icon = icon("database")),
      menuItem("Antibiotic Resistance", tabName = "resistance", icon = icon("chart-bar")),
      menuItem("MIC Analysis", tabName = "mic", icon = icon("chart-line")),
      menuItem("Testing Methods", tabName = "methods", icon = icon("flask")),
      menuItem("Raw Data", tabName = "rawdata", icon = icon("table"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .box {
          border-top: 3px solid #3c8dbc;
        }
        .shiny-file-input-progress {
          display: none;
        }
      "))
    ),
    
    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
              fluidRow(
                box(title = "Upload Your AMR Data", status = "primary", solidHeader = TRUE,
                    width = 12,
                    h4("Step 1: Upload your CSV file"),
                    fileInput("file_upload", "Choose CSV File",
                              accept = c(".csv"),
                              buttonLabel = "Browse...",
                              placeholder = "No file selected"),
                    
                    h4("Step 2: Check file info"),
                    verbatimTextOutput("file_info"),
                    
                    h4("Step 3: Load data"),
                    actionButton("load_data", "Load Data", class = "btn-success"),
                    actionButton("load_sample", "Use Sample Data", class = "btn-warning"),
                    
                    br(), br(),
                    h4("Upload Status:"),
                    verbatimTextOutput("upload_status"),
                    
                    helpText("Note: File upload limit is 100MB"),
                    helpText("If upload fails, try the sample data first to test the app")
                )
              )
      ),
      
      # Data Overview Tab
      tabItem(tabName = "overview",
              uiOutput("overview_ui")
      ),
      
      # Antibiotic Resistance Tab
      tabItem(tabName = "resistance",
              uiOutput("resistance_ui")
      ),
      
      # MIC Analysis Tab
      tabItem(tabName = "mic",
              uiOutput("mic_ui")
      ),
      
      # Testing Methods Tab
      tabItem(tabName = "methods",
              uiOutput("methods_ui")
      ),
      
      # Raw Data Tab
      tabItem(tabName = "rawdata",
              uiOutput("rawdata_ui")
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive value to store the data
  amr_data <- reactiveVal(NULL)
  
  # Show file information when selected
  output$file_info <- renderText({
    file <- input$file_upload
    if (is.null(file)) {
      return("No file selected yet.\nPlease choose a CSV file")
    }
    
    paste("File selected:", file$name, "\n",
          "Size:", round(file$size/1024/1024, 2), "MB\n",
          "Type:", file$type, "\n",
          "Click 'Load Data' to process the file")
  })
  
  # Observe load data button
  observeEvent(input$load_data, {
    file <- input$file_upload
    if (is.null(file)) {
      output$upload_status <- renderText({
        "ERROR: No file selected. Please choose a file first."
      })
      return()
    }
    
    # Show loading message
    output$upload_status <- renderText({
      paste("Loading data... This may take a moment for large files.\n",
            "File size:", round(file$size/1024/1024, 2), "MB")
    })
    
    tryCatch({
      # Read the uploaded CSV file
      data <- read_csv(file$datapath, show_col_types = FALSE)
      
      # Clean the data
      data_clean <- clean_amr_data(data)
      
      # Store the cleaned data
      amr_data(data_clean)
      
      # Success message
      output$upload_status <- renderText({
        paste("SUCCESS: Data loaded and processed!\n",
              "Records:", format(nrow(data_clean), big.mark = ","), "\n",
              "Columns:", ncol(data_clean), "\n",
              "Unique antibiotics:", length(unique(data_clean$Antibiotic)), "\n",
              "You can now navigate to other tabs.")
      })
      
    }, error = function(e) {
      output$upload_status <- renderText({
        paste("ERROR loading file:\n", e$message, "\n",
              "Try using the sample data or check your file format.")
      })
    })
  })
  
  # Data cleaning function
  clean_amr_data <- function(data) {
    # Create a copy
    clean_data <- data
    
    # Standardize column names - handle spaces
    colnames(clean_data) <- gsub(" ", "_", colnames(clean_data))
    
    # Clean each column
    clean_data <- clean_data %>%
      mutate(
        # Convert all columns to character first
        across(everything(), as.character),
        
        # Clean resistant phenotype
        Resistant_Phenotype = case_when(
          is.na(Resistant_Phenotype) | Resistant_Phenotype == "" ~ "Not Specified",
          TRUE ~ Resistant_Phenotype
        ),
        
        # Clean laboratory typing method
        Laboratory_Typing_Method = case_when(
          is.na(Laboratory_Typing_Method) | Laboratory_Typing_Method == "" ~ "Not Specified",
          TRUE ~ Laboratory_Typing_Method
        ),
        
        # Clean testing standard
        Testing_Standard = case_when(
          is.na(Testing_Standard) | Testing_Standard == "" ~ "Not Specified",
          TRUE ~ Testing_Standard
        ),
        
        # Convert measurement value to numeric safely
        Measurement_Value = suppressWarnings(as.numeric(Measurement_Value)),
        
        # Handle testing standard year
        Testing_Standard_Year = suppressWarnings(as.numeric(Testing_Standard_Year))
      )
    
    return(clean_data)
  }
  
  # Sample data for testing
  observeEvent(input$load_sample, {
    # Create sample data based on original structure
    set.seed(123)
    n_records <- 1000
    
    sample_data <- data.frame(
      Taxon_ID = rep(562, n_records),
      Genome_ID = paste0("562.", sample(100000:999999, n_records)),
      Genome_Name = paste("Escherichia coli", sample(c("strain", "isolate"), n_records, replace = TRUE), sample(1000:9999, n_records)),
      Antibiotic = sample(c("ceftriaxone", "trimethoprim/sulfamethoxazole", "meropenem", 
                            "cefepime", "piperacillin/tazobactam", "ampicillin"), 
                          n_records, replace = TRUE),
      Resistant_Phenotype = sample(c("Resistant", "Susceptible", "Intermediate", "Not Specified"), 
                                   n_records, replace = TRUE, prob = c(0.3, 0.6, 0.05, 0.05)),
      Measurement = sample(c("<=0.25", "<=0.12", "<=4", ">=32", "16", "8", "2", "1"), 
                           n_records, replace = TRUE),
      Measurement_Sign = sample(c("<=", ">=", "=", ""), n_records, replace = TRUE),
      Measurement_Value = sample(c(0.12, 0.25, 1, 2, 4, 8, 16, 32), n_records, replace = TRUE),
      Measurement_Unit = sample(c("mg/L", "mm"), n_records, replace = TRUE),
      Laboratory_Typing_Method = sample(c("MIC", "Disk diffusion", "Broth dilution"), 
                                        n_records, replace = TRUE),
      Testing_Standard = sample(c("CLSI", "EUCAST"), n_records, replace = TRUE),
      Testing_Standard_Year = sample(2018:2023, n_records, replace = TRUE),
      Computational_Method = "",
      Computational_Method_Performance = "",
      Evidence = "",
      Source = "Laboratory Method",
      PubMed = sample(26142410:38052776, n_records, replace = TRUE)
    )
    
    amr_data(sample_data)
    output$upload_status <- renderText({
      paste("SAMPLE DATA LOADED!\n",
            "Records:", nrow(sample_data), "\n",
            "This is sample data for testing. You can navigate to other tabs now.")
    })
  })
  
  # UI for Overview Tab (only shows when data is loaded)
  output$overview_ui <- renderUI({
    if (is.null(amr_data())) {
      return(
        box(title = "No Data Loaded", status = "warning", solidHeader = TRUE, width = 12,
            h4("Please upload your data first:"),
            p("1. Go to the 'Data Upload' tab"),
            p("2. Upload your CSV file or use sample data"),
            p("3. Click 'Load Data'")
        )
      )
    }
    
    tagList(
      fluidRow(
        box(title = "Dataset Summary", status = "primary", solidHeader = TRUE,
            width = 12,
            valueBoxOutput("totalRecords", width = 3),
            valueBoxOutput("uniqueAntibiotics", width = 3),
            valueBoxOutput("resistantPhenotypes", width = 3),
            valueBoxOutput("testingMethods", width = 3)
        )
      ),
      fluidRow(
        box(title = "Data Structure", status = "info", solidHeader = TRUE,
            width = 6,
            tableOutput("data_structure")
        ),
        box(title = "Data Quality Check", status = "warning", solidHeader = TRUE,
            width = 6,
            plotlyOutput("missingDataPlot", height = "300px")
        )
      ),
      fluidRow(
        box(title = "First 10 Rows of Data", status = "info", solidHeader = TRUE,
            width = 12,
            DTOutput("data_preview")
        )
      )
    )
  })
  
  # UI for Resistance Tab
  output$resistance_ui <- renderUI({
    if (is.null(amr_data())) {
      return(
        box(title = "No Data Loaded", status = "warning", solidHeader = TRUE, width = 12,
            h4("Please upload your data first in the 'Data Upload' tab")
        )
      )
    }
    
    data <- amr_data()
    antibiotics <- unique(data$Antibiotic)
    
    tagList(
      fluidRow(
        box(title = "Filters", status = "primary", solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(4,
                     selectizeInput("antibioticFilter", "Select Antibiotics:",
                                    choices = antibiotics, 
                                    multiple = TRUE,
                                    selected = if(length(antibiotics) > 0) antibiotics[1:min(3, length(antibiotics))] else NULL)
              ),
              column(4,
                     selectInput("methodFilter", "Testing Method:",
                                 choices = c("All", unique(data$Laboratory_Typing_Method)))
              ),
              column(4,
                     selectInput("standardFilter", "Testing Standard:",
                                 choices = c("All", unique(data$Testing_Standard)))
              )
            )
        )
      ),
      fluidRow(
        box(title = "Resistance Profile", status = "info", solidHeader = TRUE,
            width = 6,
            plotlyOutput("resistanceBarPlot", height = "400px")
        ),
        box(title = "Resistance by Testing Method", status = "info", solidHeader = TRUE,
            width = 6,
            plotlyOutput("methodComparisonPlot", height = "400px")
        )
      ),
      fluidRow(
        box(title = "Top Antibiotics - Resistance Rates", status = "info", solidHeader = TRUE,
            width = 12,
            plotlyOutput("resistanceRatesPlot", height = "400px")
        )
      )
    )
  })
  
  # UI for MIC Tab
  output$mic_ui <- renderUI({
    if (is.null(amr_data())) {
      return(
        box(title = "No Data Loaded", status = "warning", solidHeader = TRUE, width = 12,
            h4("Please upload your data first in the 'Data Upload' tab")
        )
      )
    }
    
    data <- amr_data()
    
    tagList(
      fluidRow(
        box(title = "MIC Analysis Controls", status = "primary", solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(6,
                     selectInput("micAntibiotic", "Select Antibiotic:",
                                 choices = unique(data$Antibiotic))
              ),
              column(6,
                     checkboxInput("logTransform", "Log Transform MIC Values", 
                                   value = TRUE)
              )
            )
        )
      ),
      fluidRow(
        box(title = "MIC Distribution Plot", status = "info", solidHeader = TRUE,
            width = 6,
            plotlyOutput("micDistributionPlot", height = "400px")
        ),
        box(title = "MIC Statistics", status = "info", solidHeader = TRUE,
            width = 6,
            tableOutput("micStatsTable")
        )
      ),
      fluidRow(
        box(title = "MIC vs Resistance Correlation", status = "info", solidHeader = TRUE,
            width = 12,
            plotlyOutput("micCorrelationPlot", height = "400px")
        )
      )
    )
  })
  
  # UI for Methods Tab
  output$methods_ui <- renderUI({
    if (is.null(amr_data())) {
      return(
        box(title = "No Data Loaded", status = "warning", solidHeader = TRUE, width = 12,
            h4("Please upload your data first in the 'Data Upload' tab")
        )
      )
    }
    
    tagList(
      fluidRow(
        box(title = "Testing Methods Overview", status = "primary", solidHeader = TRUE,
            width = 6,
            plotlyOutput("methodDistributionPlot", height = "400px")
        ),
        box(title = "Testing Standards", status = "primary", solidHeader = TRUE,
            width = 6,
            plotlyOutput("standardDistributionPlot", height = "400px")
        )
      ),
      fluidRow(
        box(title = "Resistance by Testing Standard", status = "info", solidHeader = TRUE,
            width = 12,
            plotlyOutput("resistanceByStandardPlot", height = "400px")
        )
      )
    )
  })
  
  # UI for Raw Data Tab
  output$rawdata_ui <- renderUI({
    if (is.null(amr_data())) {
      return(
        box(title = "No Data Loaded", status = "warning", solidHeader = TRUE, width = 12,
            h4("Please upload your data first in the 'Data Upload' tab")
        )
      )
    }
    
    box(title = "Raw AMR Data", status = "primary", solidHeader = TRUE,
        width = 12,
        DTOutput("rawDataTable")
    )
  })
  
  # Data structure table
  output$data_structure <- renderTable({
    data <- amr_data()
    req(data)
    
    data.frame(
      Column = colnames(data),
      Type = sapply(data, class),
      Non_NA = sapply(data, function(x) sum(!is.na(x))),
      Unique_Values = sapply(data, function(x) length(unique(x))),
      stringsAsFactors = FALSE
    )
  })
  
  # Data preview
  output$data_preview <- renderDT({
    data <- amr_data()
    req(data)
    
    datatable(head(data, 10), 
              options = list(
                scrollX = TRUE,
                dom = 't'
              ))
  })
  
  # Filtered data based on user selections
  filtered_data <- reactive({
    data <- amr_data()
    req(data)
    
    filtered <- data
    
    if (!is.null(input$antibioticFilter) && length(input$antibioticFilter) > 0) {
      filtered <- filtered %>% filter(Antibiotic %in% input$antibioticFilter)
    }
    
    if (!is.null(input$methodFilter) && input$methodFilter != "All") {
      filtered <- filtered %>% filter(Laboratory_Typing_Method == input$methodFilter)
    }
    
    if (!is.null(input$standardFilter) && input$standardFilter != "All") {
      filtered <- filtered %>% filter(Testing_Standard == input$standardFilter)
    }
    
    return(filtered)
  })
  
  # Value boxes for overview
  output$totalRecords <- renderValueBox({
    data <- amr_data()
    req(data)
    valueBox(
      format(nrow(data), big.mark = ","), "Total Records", icon = icon("database"),
      color = "purple"
    )
  })
  
  output$uniqueAntibiotics <- renderValueBox({
    data <- amr_data()
    req(data)
    valueBox(
      length(unique(data$Antibiotic)), "Unique Antibiotics", icon = icon("pills"),
      color = "green"
    )
  })
  
  output$resistantPhenotypes <- renderValueBox({
    data <- amr_data()
    req(data)
    resistant_count <- sum(data$Resistant_Phenotype == "Resistant", na.rm = TRUE)
    valueBox(
      format(resistant_count, big.mark = ","), "Resistant Phenotypes", icon = icon("exclamation-triangle"),
      color = "red"
    )
  })
  
  output$testingMethods <- renderValueBox({
    data <- amr_data()
    req(data)
    valueBox(
      length(unique(data$Laboratory_Typing_Method)), "Testing Methods", icon = icon("flask"),
      color = "yellow"
    )
  })
  
  # Missing data plot
  output$missingDataPlot <- renderPlotly({
    data <- amr_data()
    req(data)
    
    missing_summary <- data %>%
      summarise(across(everything(), ~sum(is.na(.) | . == "" | . == "Not Specified"))) %>%
      pivot_longer(everything(), names_to = "Variable", values_to = "Missing_Count") %>%
      mutate(Missing_Percentage = round(Missing_Count / nrow(data) * 100, 1))
    
    p <- ggplot(missing_summary, aes(x = reorder(Variable, -Missing_Percentage), y = Missing_Percentage)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      labs(title = "Missing Data by Variable", x = "Variable", y = "Missing Data (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p)
  })
  
  # Resistance bar plot
  output$resistanceBarPlot <- renderPlotly({
    data <- filtered_data()
    req(data)
    
    resistance_summary <- data %>%
      filter(Resistant_Phenotype %in% c("Resistant", "Susceptible", "Intermediate")) %>%
      count(Antibiotic, Resistant_Phenotype) %>%
      group_by(Antibiotic) %>%
      mutate(Percentage = n / sum(n) * 100) %>%
      ungroup()
    
    if(nrow(resistance_summary) == 0) {
      return(plotly_empty(type = "scatter") %>% 
               layout(title = "No data available for selected filters"))
    }
    
    p <- ggplot(resistance_summary, aes(x = Antibiotic, y = Percentage, fill = Resistant_Phenotype)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("Resistant" = "red", "Susceptible" = "green", "Intermediate" = "orange")) +
      labs(title = "Antibiotic Resistance Profile", x = "Antibiotic", y = "Percentage") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p)
  })
  
  # Resistance rates plot
  output$resistanceRatesPlot <- renderPlotly({
    data <- amr_data()
    req(data)
    
    resistance_rates <- data %>%
      filter(Resistant_Phenotype %in% c("Resistant", "Susceptible")) %>%
      count(Antibiotic, Resistant_Phenotype) %>%
      group_by(Antibiotic) %>%
      mutate(Total = sum(n),
             Resistance_Rate = ifelse(Resistant_Phenotype == "Resistant", n / Total * 100, 0)) %>%
      filter(Resistant_Phenotype == "Resistant") %>%
      ungroup() %>%
      arrange(desc(Resistance_Rate)) %>%
      head(10)
    
    if(nrow(resistance_rates) == 0) {
      return(plotly_empty(type = "scatter") %>% 
               layout(title = "No resistance data available"))
    }
    
    p <- ggplot(resistance_rates, aes(x = reorder(Antibiotic, Resistance_Rate), y = Resistance_Rate)) +
      geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
      coord_flip() +
      labs(title = "Top Antibiotics by Resistance Rate", 
           x = "Antibiotic", 
           y = "Resistance Rate (%)") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Method comparison plot
  output$methodComparisonPlot <- renderPlotly({
    data <- filtered_data()
    req(data)
    
    method_summary <- data %>%
      filter(Resistant_Phenotype == "Resistant") %>%
      count(Laboratory_Typing_Method, Antibiotic)
    
    if(nrow(method_summary) == 0) {
      return(plotly_empty(type = "scatter") %>% 
               layout(title = "No resistant cases for selected filters"))
    }
    
    p <- ggplot(method_summary, aes(x = Antibiotic, y = n, fill = Laboratory_Typing_Method)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "Resistance Cases by Testing Method", 
           x = "Antibiotic", 
           y = "Number of Resistant Cases") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p)
  })
  
  # MIC distribution plot
  output$micDistributionPlot <- renderPlotly({
    data <- filtered_data()
    req(data)
    antibiotic <- input$micAntibiotic
    
    if (!is.null(antibiotic) && antibiotic != "") {
      mic_data <- data %>%
        filter(Antibiotic == antibiotic & !is.na(Measurement_Value))
      
      if (nrow(mic_data) > 0) {
        # Calculate resistance rate for this antibiotic
        resistance_rate <- mic_data %>%
          filter(Resistant_Phenotype %in% c("Resistant", "Susceptible")) %>%
          summarise(Resistance_Rate = sum(Resistant_Phenotype == "Resistant") / n() * 100)
        
        if (input$logTransform) {
          mic_data$Log_Value <- log10(mic_data$Measurement_Value + 1)
          y_lab <- "Log10(MIC Value + 1)"
          values <- mic_data$Log_Value
        } else {
          y_lab <- "MIC Value (mg/L)"
          values <- mic_data$Measurement_Value
        }
        
        p <- plot_ly(mic_data, y = ~values, type = "box",
                     boxpoints = "all", jitter = 0.3, pointpos = -1.8,
                     marker = list(size = 3, opacity = 0.5),
                     line = list(color = 'rgb(107,174,214)'),
                     hoverinfo = "y") %>%
          layout(title = paste("MIC Distribution for", antibiotic, 
                               "<br>Resistance Rate:", round(resistance_rate$Resistance_Rate, 1), "%"),
                 yaxis = list(title = y_lab),
                 xaxis = list(title = ""))
        
        p
      } else {
        plotly_empty(type = "scatter") %>% 
          layout(title = paste("No MIC data available for", antibiotic))
      }
    }
  })
  
  # MIC statistics table
  output$micStatsTable <- renderTable({
    data <- filtered_data()
    req(data)
    antibiotic <- input$micAntibiotic
    
    if (!is.null(antibiotic) && antibiotic != "") {
      stats <- data %>%
        filter(Antibiotic == antibiotic & !is.na(Measurement_Value)) %>%
        summarise(
          Antibiotic = antibiotic,
          Count = n(),
          Mean_MIC = round(mean(Measurement_Value, na.rm = TRUE), 3),
          Median_MIC = round(median(Measurement_Value, na.rm = TRUE), 3),
          SD = round(sd(Measurement_Value, na.rm = TRUE), 3),
          Min_MIC = round(min(Measurement_Value, na.rm = TRUE), 3),
          Max_MIC = round(max(Measurement_Value, na.rm = TRUE), 3),
          Resistance_Rate = round(sum(Resistant_Phenotype == "Resistant", na.rm = TRUE) / 
                                    sum(Resistant_Phenotype %in% c("Resistant", "Susceptible"), na.rm = TRUE) * 100, 1)
        )
      
      stats
    }
  }, bordered = TRUE)
  
  # MIC Correlation Plot
  output$micCorrelationPlot <- renderPlotly({
    data <- amr_data()
    req(data)
    correlation_data<- data %>%
      filter(!is.na(Measurement_Value) & 
               Resistant_Phenotype %in% c("Resistant", "Susceptible", "Intermediate")) %>%
      group_by(Antibiotic) %>%
      summarise(
        Avg_MIC = mean(Measurement_Value, na.rm = TRUE),
        Total_Isolates = n(),
        Resistant_Count = sum(Resistant_Phenotype == "Resistant", na.rm = TRUE),
        Resistance_Rate = (Resistant_Count / Total_Isolates) * 100,
        .groups = 'drop'
      ) %>%
      filter(Total_Isolates >= 10) 
    if(nrow(correlation_data) == 0) {
      return(plotly_empty(type = "scatter") %>% 
               layout(title = "No sufficient data for correlation analysis"))
    }
    
    # Calculating correlation coefficient
    cor_coef <- cor(correlation_data$Avg_MIC, correlation_data$Resistance_Rate, 
                    method = "pearson")
    
    p <- ggplot(correlation_data, 
                aes(x = Avg_MIC, y = Resistance_Rate, 
                    text = paste("Antibiotic:", Antibiotic,
                                 "\nAvg MIC:", round(Avg_MIC, 2),
                                 "\nResistance Rate:", round(Resistance_Rate, 1), "%",
                                 "\nIsolates:", Total_Isolates))) +
      geom_point(aes(size = Total_Isolates, color = Antibiotic), alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "darkred", linetype = "dashed", alpha = 0.3) +
      labs(title = paste("MIC Values vs Resistance Rates (Correlation =", round(cor_coef, 2), ")"),
           x = "Average MIC Value (mg/L)",
           y = "Resistance Rate (%)",
           size = "Number of\nIsolates") +
      theme_minimal() +
      theme(legend.position = "right")
    
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "v", y = 0.5))
  })
  
  # Method distribution plot
  output$methodDistributionPlot <- renderPlotly({
    data <- amr_data()
    req(data)
    
    method_summary <- data %>%
      count(Laboratory_Typing_Method) %>%
      mutate(Percentage = n / sum(n) * 100)
    
    p <- plot_ly(method_summary, labels = ~Laboratory_Typing_Method, values = ~n, type = 'pie') %>%
      layout(title = "Distribution of Testing Methods")
    
    p
  })
  
  # Standard distribution plot
  output$standardDistributionPlot <- renderPlotly({
    data <- amr_data()
    req(data)
    
    standard_summary <- data %>%
      count(Testing_Standard) %>%
      mutate(Percentage = n / sum(n) * 100)
    
    p <- plot_ly(standard_summary, labels = ~Testing_Standard, values = ~n, type = 'pie') %>%
      layout(title = "Distribution of Testing Standards")
    
    p
  })
  
  # Resistance by standard plot
  output$resistanceByStandardPlot <- renderPlotly({
    data <- amr_data()
    req(data)
    
    standard_resistance <- data %>%
      filter(Resistant_Phenotype %in% c("Resistant", "Susceptible")) %>%
      count(Testing_Standard, Resistant_Phenotype) %>%
      group_by(Testing_Standard) %>%
      mutate(Percentage = n / sum(n) * 100) %>%
      ungroup()
    
    p <- ggplot(standard_resistance, aes(x = Testing_Standard, y = Percentage, fill = Resistant_Phenotype)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("Resistant" = "red", "Susceptible" = "green")) +
      labs(title = "Resistance Profile by Testing Standard", 
           x = "Testing Standard", 
           y = "Percentage") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Raw data table
  output$rawDataTable <- renderDT({
    data <- amr_data()
    req(data)
    
    datatable(data, 
              options = list(
                scrollX = TRUE,
                pageLength = 10,
                lengthMenu = c(10, 25, 50, 100)
              ),
              filter = 'top')
  })
}

# Run the application
shinyApp(ui = ui, server = server)