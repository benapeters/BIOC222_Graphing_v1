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
  
  # UNITS CONFIGURATION - Change units here to update throughout the app
  conc_units <- "µmol.mL⁻¹"
  rate_units <- "µmol.mL⁻¹.min⁻¹"
  
  # Create display labels with units
  conc_label <- paste0("Substrate Concentration (", conc_units, ")")
  rate_label <- paste0("Initial Rate (", rate_units, ")")
  
  # Define the initial data frame for Tab 2 (Km and Vmax calculation)
  kinetics_data <- reactiveValues(df = data.frame(
    `Substrate Concentration` = c(1, 0.5, 0.25, 0.125, 0.063, 0.031, 0.016,0.008),
    `Initial Rate` = rep(NA, 8),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  
  # Convert kinetics data to numeric
  observe({
    kinetics_data$df$`Substrate Concentration` <- as.numeric(kinetics_data$df$`Substrate Concentration`)
    kinetics_data$df$`Initial Rate` <- as.numeric(kinetics_data$df$`Initial Rate`)
  })
  
  # Render the kinetics table for Tab 2
  output$kinetics_table <- renderRHandsontable({
    # Create a copy of the data with proper column names for display
    display_df <- kinetics_data$df
    names(display_df) <- c(conc_label, rate_label)
    
    rhandsontable(display_df, rowHeaders = FALSE) %>% 
      hot_col(conc_label, type = "numeric", strict = TRUE, allowInvalid = FALSE, format = 0.000) %>%
      hot_col(rate_label, type = "numeric", strict = TRUE, allowInvalid = FALSE, format = 0.000)
  })
  
  # Update kinetics dataframe after inputs for Tab 2
  observe({
    if (!is.null(input$kinetics_table)) {
      temp_df <- hot_to_r(input$kinetics_table)
      # Restore original column names for internal use
      names(temp_df) <- c("Substrate Concentration", "Initial Rate")
      kinetics_data$df <- temp_df
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
  
  # Display Michaelis-Menten parameters table
  output$kinetics_params_table <- renderTable({
    mm_data <- michaelis_menten_fit()
    
    if (is.null(mm_data)) {
      return(data.frame(Message = "Enter at least 3 complete data points to calculate Km and Vmax"))
    }
    
    params_table <- data.frame(
      Parameter = c("Vmax", "Km", "R-squared"),
      Value = c(round(mm_data$Vmax, 4), round(mm_data$Km, 4), round(mm_data$r_squared, 4)),
      Units = c(rate_units, conc_units, ""),
      stringsAsFactors = FALSE
    )
    
    return(params_table)
  }, digits = 4)
  
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
    stats_text <- paste0(stats_text, "Vmax = ", round(mm_data$Vmax, 4), " ", rate_units, " (maximum velocity)\n")
    stats_text <- paste0(stats_text, "Km = ", round(mm_data$Km, 4), " ", conc_units, " (substrate concentration at Vmax/2)\n")
    stats_text <- paste0(stats_text, "R² = ", round(mm_data$r_squared, 4))
    
    return(stats_text)
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
      # Create the plot with units in axis labels
      p <- ggplot(mm_data$data, aes(x = S, y = V)) +
        geom_point(color = "black", size = 4) +
        geom_line(data = mm_data$pred_data, aes(x = S, y = V), color = "blue", size = 1.2) +
        labs(x = paste0("Substrate Concentration [S] (", conc_units, ")"), 
             y = paste0("Initial Velocity (V) (", rate_units, ")"), 
             title = "Michaelis-Menten Plot") +
        theme_minimal() +
        theme(axis.line.x = element_line(color = "black", size = 1),
              axis.line.y = element_line(color = "black", size = 1),
              plot.title = element_text(hjust = 0.5, size = 16),
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12)) +
        theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
              axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
        scale_y_continuous(limits = c(0, NA))
      
      # Add Vmax and Km reference lines with units in annotations
      p <- p + 
        geom_hline(yintercept = mm_data$Vmax, linetype = "dashed", color = "red", alpha = 0.7) +
        geom_vline(xintercept = mm_data$Km, linetype = "dashed", color = "red", alpha = 0.7) +
        geom_hline(yintercept = mm_data$Vmax/2, linetype = "dashed", color = "orange", alpha = 0.7) +
        annotate("text", x = max(mm_data$data$S) * 0.7, y = mm_data$Vmax * 1.05, 
                 label = paste("Vmax =", round(mm_data$Vmax, 3), rate_units), color = "red", size = 4) +
        annotate("text", x = mm_data$Km * 1.1, y = max(mm_data$data$V) * 0.9, 
                 label = paste("Km =", round(mm_data$Km, 3), conc_units), color = "red", size = 4) +
        annotate("text", x = max(mm_data$data$S) * 0.7, y = (mm_data$Vmax/2) * 1.1, 
                 label = paste("Vmax/2 =", round(mm_data$Vmax/2, 3), rate_units), color = "orange", size = 4)
      
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