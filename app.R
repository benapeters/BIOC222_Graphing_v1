library(shiny)
library(ggplot2)
library(rhandsontable)
library(minpack.lm)  # For non-linear least squares fitting

#UI section
ui <- fluidPage(
  
  # Application title
  titlePanel("BIOC222 Enzyme Kinetics Calculator"),
  
  tabsetPanel(
    
    tabPanel(
      "Initial rate calculation",
      fluidRow(
        column(6,rHandsontableOutput("table1")),
        column(6,
               h4("Initial Rate (Slope)"),
               tableOutput("initial_rates_table"),
               br(),
               h4("Linear Model Statistics"),
               verbatimTextOutput("model_stats")
        )
      ),
      fluidRow(
        column(12,"Move the slider so that only the initial linear portion of the graph is used for the line of best fit")
      ),
      sliderInput("slider_id", "", min = 0, max = 120, step = 15, value =c(15,60)),
      plotOutput("progressCurve")
    ),
    
    tabPanel(
      "Km and Vmax Calculation",
      fluidRow(
        column(6,
               h4("Substrate Concentration and Initial Rates"),
               rHandsontableOutput("kinetics_table"),
               br(),
               actionButton("add_row", "Add Row", class = "btn-primary"),
               actionButton("remove_row", "Remove Row", class = "btn-warning")
        ),
        column(6,
               h4("Michaelis-Menten Parameters"),
               tableOutput("kinetics_params_table"),
               br(),
               h4("Model Statistics"),
               verbatimTextOutput("kinetics_model_stats")
        )
      ),
      br(),
      fluidRow(
        column(12,
               plotOutput("michaelis_menten_plot", height = "500px")
        )
      )
    ),
    
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Define the initial data frame for Tab 1 (Initial rate calculation)
  data1 <- reactiveValues(df = data.frame(
    Time = seq(0, 120, by = 15),
    Assay = rep(NA, 9),
    stringsAsFactors = FALSE
  ))
  
  # Define the initial data frame for Tab 2 (Km and Vmax calculation)
  kinetics_data <- reactiveValues(df = data.frame(
    `Substrate Concentration` = c(1, 0.5, 0.25, 0.125, 0.063, 0.031, 0.016,0.008),
    `Initial Rate` = rep(NA, 8),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  
  #The slider text
  output$slider_label <- renderText({
    "Move the slider so that only the initial linear portion of the graph is used for the line of best fit"
  })
  
  # Convert Assay to numeric (if possible)
  observe({
    data1$df$Assay <- as.numeric(data1$df$Assay)
  })
  
  # Convert kinetics data to numeric
  observe({
    kinetics_data$df$`Substrate Concentration` <- as.numeric(kinetics_data$df$`Substrate Concentration`)
    kinetics_data$df$`Initial Rate` <- as.numeric(kinetics_data$df$`Initial Rate`)
  })
  
  # Render the data table for Tab 1
  output$table1 <- renderRHandsontable({
    rhandsontable(data1$df, rowHeaders = FALSE) %>% 
      hot_col("Time", type = "numeric", strict = TRUE, allowInvalid = FALSE, readOnly = TRUE, format = 0) %>%
      hot_col("Assay", type = "numeric", strict = TRUE, allowInvalid = FALSE, format = 0.000)
  })
  
  # Render the kinetics table for Tab 2
  output$kinetics_table <- renderRHandsontable({
    rhandsontable(kinetics_data$df, rowHeaders = FALSE) %>% 
      hot_col("Substrate Concentration", type = "numeric", strict = TRUE, allowInvalid = FALSE, format = 0.000) %>%
      hot_col("Initial Rate", type = "numeric", strict = TRUE, allowInvalid = FALSE, format = 0.000)
  })
  
  # Update dataframe after inputs for Tab 1
  observe({
    if (!is.null(input$table1)) {
      data1$df <- hot_to_r(input$table1)
    }
  })
  
  # Update kinetics dataframe after inputs for Tab 2
  observe({
    if (!is.null(input$kinetics_table)) {
      kinetics_data$df <- hot_to_r(input$kinetics_table)
    }
  })
  
  # Add row functionality for kinetics table
  observeEvent(input$add_row, {
    new_row <- data.frame(
      `Substrate Concentration` = NA,
      `Initial Rate` = NA,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    kinetics_data$df <- rbind(kinetics_data$df, new_row)
  })
  
  # Remove row functionality for kinetics table
  observeEvent(input$remove_row, {
    if (nrow(kinetics_data$df) > 1) {
      kinetics_data$df <- kinetics_data$df[-nrow(kinetics_data$df), ]
    }
  })
  
  # Calculate initial rates and model statistics for Tab 1
  initial_rates_data <- reactive({
    req(nrow(data1$df) > 0)
    
    slider_value_min <- input$slider_id[1]
    slider_value_max <- input$slider_id[2]
    
    # Subset the data based on the slider value
    subset_data <- data1$df[data1$df$Time >= slider_value_min & data1$df$Time <= slider_value_max, ]
    
    # Check if we have valid data
    if (nrow(subset_data) < 2 || all(is.na(subset_data$Assay))) {
      return(NULL)
    }
    
    results <- list()
    
    tryCatch({
      # Remove rows with NA values
      subset_assay <- subset_data[!is.na(subset_data$Assay), ]
      
      if (nrow(subset_assay) >= 2) {
        fit <- lm(Assay ~ Time, data = subset_assay)
        results$fit <- fit
        # Convert from abs/sec to abs/min by multiplying by 60
        results$rate <- coef(fit)[2] * 60  # slope * 60 seconds/minute
        results$r_squared <- summary(fit)$r.squared
        results$p_value <- summary(fit)$coefficients[2, 4]  # p-value for slope
      }
      
      return(results)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  # Calculate Michaelis-Menten parameters
  michaelis_menten_fit <- reactive({
    req(nrow(kinetics_data$df) > 0)
    
    # Get data and remove NA values
    df <- kinetics_data$df
    df_clean <- df[complete.cases(df), ]
    
    if (nrow(df_clean) < 3) {
      return(NULL)
    }
    
    S <- df_clean$`Substrate Concentration`
    V <- df_clean$`Initial Rate`
    
    if (length(S) < 3 || any(S <= 0) || any(V <= 0)) {
      return(NULL)
    }
    
    results <- list()
    
    tryCatch({
      # Estimate initial parameters
      Vmax_est <- max(V) * 1.2
      Km_est <- S[which.min(abs(V - Vmax_est/2))]
      
      # Fit Michaelis-Menten equation: V = (Vmax * S) / (Km + S)
      fit <- nlsLM(V ~ (Vmax * S) / (Km + S), 
                   start = list(Vmax = Vmax_est, Km = Km_est),
                   data = data.frame(S = S, V = V),
                   control = nls.lm.control(maxiter = 100))
      
      results$fit <- fit
      results$Vmax <- coef(fit)["Vmax"]
      results$Km <- coef(fit)["Km"]
      results$r_squared <- 1 - sum(residuals(fit)^2) / sum((V - mean(V))^2)
      results$data <- data.frame(S = S, V = V)
      
      # Generate prediction data for smooth curve
      S_pred <- seq(min(S), max(S) * 1.2, length.out = 100)
      V_pred <- predict(fit, newdata = list(S = S_pred))
      results$pred_data <- data.frame(S = S_pred, V = V_pred)
      
      return(results)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  # Display initial rates table for Tab 1
  output$initial_rates_table <- renderTable({
    rates_data <- initial_rates_data()
    
    if (is.null(rates_data)) {
      return(data.frame(Message = "Enter data to calculate initial rate"))
    }
    
    rate_table <- data.frame(
      Assay = "Single Assay",
      `Initial Rate (abs/min)` = round(rates_data$rate, 6),
      `R-squared` = round(rates_data$r_squared, 4),
      `P-value` = ifelse(rates_data$p_value < 0.001, "<0.001", round(rates_data$p_value, 4)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    return(rate_table)
  }, digits = 6)
  
  # Display Michaelis-Menten parameters table
  output$kinetics_params_table <- renderTable({
    mm_data <- michaelis_menten_fit()
    
    if (is.null(mm_data)) {
      return(data.frame(Message = "Enter at least 3 complete data points to calculate Km and Vmax"))
    }
    
    params_table <- data.frame(
      Parameter = c("Vmax", "Km", "R-squared"),
      Value = c(round(mm_data$Vmax, 4), round(mm_data$Km, 4), round(mm_data$r_squared, 4)),
      Units = c("Initial Rate Units", "Substrate Concentration Units", ""),
      stringsAsFactors = FALSE
    )
    
    return(params_table)
  }, digits = 4)
  
  # Display model statistics for Tab 1
  output$model_stats <- renderText({
    rates_data <- initial_rates_data()
    
    if (is.null(rates_data)) {
      return("Enter data to see model statistics")
    }
    
    stats_text <- ""
    
    if (!is.null(rates_data$fit)) {
      stats_text <- paste0("Linear Model:\n")
      stats_text <- paste0(stats_text, "Equation: Absorbance = ", 
                           round(coef(rates_data$fit)[1], 6), " + ", 
                           round(coef(rates_data$fit)[2], 6), " × Time\n")
      stats_text <- paste0(stats_text, "Initial Rate: ", round(rates_data$rate, 6), " abs/min\n")
    }
    
    return(stats_text)
  })
  
  # Display kinetics model statistics
  output$kinetics_model_stats <- renderText({
    mm_data <- michaelis_menten_fit()
    
    if (is.null(mm_data)) {
      return("Enter at least 3 complete data points to see model statistics")
    }
    
    stats_text <- ""
    stats_text <- paste0("Michaelis-Menten Model:\n")
    stats_text <- paste0(stats_text, "Equation: V = (", round(mm_data$Vmax, 4), " × [S]) / (", 
                         round(mm_data$Km, 4), " + [S])\n")
    stats_text <- paste0(stats_text, "Vmax = ", round(mm_data$Vmax, 4), " (maximum velocity)\n")
    stats_text <- paste0(stats_text, "Km = ", round(mm_data$Km, 4), " (substrate concentration at Vmax/2)\n")
    stats_text <- paste0(stats_text, "R² = ", round(mm_data$r_squared, 4))
    
    return(stats_text)
  })
  
  # Render the progress Curve scatter plot for Tab 1
  output$progressCurve <- renderPlot({
    req(nrow(data1$df) > 0)
    
    #writes the slider input to the slider value
    slider_value_min <- input$slider_id[1]
    slider_value_max <- input$slider_id[2]
    
    # Subset the data based on the slider value
    subset_data <- data1$df[data1$df$Time >= slider_value_min & data1$df$Time <= slider_value_max, ]
    
    tryCatch({
      ggplot(data1$df, aes(x = Time)) +
        geom_line(aes(y = Assay), color = "grey", size = 2) +
        geom_smooth(data = subset_data, aes(y = Assay), color = "black", method = "lm", se = FALSE, fullrange = TRUE, linetype = "dotted") +
        geom_point(aes(y = Assay), color = "black", size = 3) +
        labs(x = "Time (seconds)", y = "Absorbance", title = "Initial rate of reaction") +
        theme_minimal() +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max(data1$df$Assay, na.rm = TRUE) * 1.1)) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 200), breaks = seq(0, max(data1$df$Time, na.rm = TRUE), by = 60), 
                           minor_breaks = seq(0, max(data1$df$Time, na.rm = TRUE), by = 20)) +
        theme(axis.line.x = element_line(color = "black", size = 1),
              axis.line.y = element_line(color = "black", size = 1),
              plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_line(colour = "grey", linetype = "solid"),
              axis.line = element_line(colour = "black", size = 1, linetype = "solid"))
      
      
    }, error = function(e){"The plot will appear once you have entered your data"})
  })
  
  # Render Michaelis-Menten plot for Tab 2
  output$michaelis_menten_plot <- renderPlot({
    mm_data <- michaelis_menten_fit()
    
    if (is.null(mm_data)) {
      return(ggplot() + 
               geom_text(aes(x = 0.5, y = 0.5, label = "Enter at least 3 complete data points to generate plot"), 
                         size = 6) +
               theme_void() +
               xlim(0, 1) + ylim(0, 1))
    }
    
    tryCatch({
      # Create the plot
      p <- ggplot(mm_data$data, aes(x = S, y = V)) +
        geom_point(color = "black", size = 4) +
        geom_line(data = mm_data$pred_data, aes(x = S, y = V), color = "blue", size = 1.2) +
        labs(x = "Substrate Concentration [S]", 
             y = "Initial Velocity (V)", 
             title = "Michaelis-Menten Plot") +
        theme_minimal() +
        theme(axis.line.x = element_line(color = "black", size = 1),
              axis.line.y = element_line(color = "black", size = 1),
              plot.title = element_text(hjust = 0.5, size = 16),
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12)) +
        theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
              axis.line = element_line(colour = "black", size = 1, linetype = "solid"))
      
      # Add Vmax and Km reference lines
      p <- p + 
        geom_hline(yintercept = mm_data$Vmax, linetype = "dashed", color = "red", alpha = 0.7) +
        geom_vline(xintercept = mm_data$Km, linetype = "dashed", color = "red", alpha = 0.7) +
        geom_hline(yintercept = mm_data$Vmax/2, linetype = "dashed", color = "orange", alpha = 0.7) +
        annotate("text", x = max(mm_data$data$S) * 0.7, y = mm_data$Vmax * 1.05, 
                 label = paste("Vmax =", round(mm_data$Vmax, 3)), color = "red", size = 4) +
        annotate("text", x = mm_data$Km * 1.1, y = max(mm_data$data$V) * 0.9, 
                 label = paste("Km =", round(mm_data$Km, 3)), color = "red", size = 4) +
        annotate("text", x = max(mm_data$data$S) * 0.7, y = (mm_data$Vmax/2) * 1.1, 
                 label = paste("Vmax/2 =", round(mm_data$Vmax/2, 3)), color = "orange", size = 4)
      
      return(p)
      
    }, error = function(e) {
      return(ggplot() + 
               geom_text(aes(x = 0.5, y = 0.5, label = "Error generating plot. Check your data."), 
                         size = 6, color = "red") +
               theme_void() +
               xlim(0, 1) + ylim(0, 1))
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)