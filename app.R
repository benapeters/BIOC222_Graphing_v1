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
        column(6,
               h4("Michaelis-Menten Plot"),
               plotOutput("michaelis_menten_plot", height = "400px")
        ),
        column(6,
               h4("Lineweaver-Burk Plot"),
               plotOutput("lineweaver_burk_plot", height = "400px")
        )
      ),
      br(),
      fluidRow(
        column(12,
               h4("Lineweaver-Burk Analysis"),
               tableOutput("lineweaver_params_table"),
               br(),
               verbatimTextOutput("lineweaver_model_stats")
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
  
  # Calculate Lineweaver-Burk parameters
  lineweaver_burk_fit <- reactive({
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
      # Calculate reciprocals
      inv_S <- 1/S
      inv_V <- 1/V
      
      # Linear regression: 1/V = (Km/Vmax) * (1/S) + (1/Vmax)
      lm_fit <- lm(inv_V ~ inv_S)
      
      # Extract parameters
      slope <- coef(lm_fit)[2]  # Km/Vmax
      intercept <- coef(lm_fit)[1]  # 1/Vmax
      
      results$fit <- lm_fit
      results$Vmax <- 1/intercept
      results$Km <- slope * results$Vmax
      results$r_squared <- summary(lm_fit)$r.squared
      results$data <- data.frame(inv_S = inv_S, inv_V = inv_V, S = S, V = V)
      
      # Generate prediction data for line
      inv_S_pred <- seq(min(inv_S), max(inv_S), length.out = 100)
      inv_V_pred <- predict(lm_fit, newdata = data.frame(inv_S = inv_S_pred))
      results$pred_data <- data.frame(inv_S = inv_S_pred, inv_V = inv_V_pred)
      
      # Calculate x and y intercepts
      results$x_intercept <- -intercept/slope  # -1/Km
      results$y_intercept <- intercept  # 1/Vmax
      
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
  
  # Display Lineweaver-Burk parameters table
  output$lineweaver_params_table <- renderTable({
    lb_data <- lineweaver_burk_fit()
    
    if (is.null(lb_data)) {
      return(data.frame(Message = "Enter at least 3 complete data points to calculate Lineweaver-Burk parameters"))
    }
    
    params_table <- data.frame(
      Parameter = c("Vmax (L-B)", "Km (L-B)", "R-squared (L-B)", "Slope", "Y-intercept", "X-intercept"),
      Value = c(round(lb_data$Vmax, 4), round(lb_data$Km, 4), round(lb_data$r_squared, 4),
                round(coef(lb_data$fit)[2], 4), round(lb_data$y_intercept, 4), round(lb_data$x_intercept, 4)),
      Units = c(rate_units, conc_units, "", paste0(conc_units, "/", rate_units), 
                paste0("1/", rate_units), paste0("-1/", conc_units)),
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
  
  # Display Lineweaver-Burk model statistics
  output$lineweaver_model_stats <- renderText({
    lb_data <- lineweaver_burk_fit()
    
    if (is.null(lb_data)) {
      return("Enter at least 3 complete data points to see Lineweaver-Burk statistics")
    }
    
    slope <- coef(lb_data$fit)[2]
    intercept <- coef(lb_data$fit)[1]
    
    stats_text <- ""
    stats_text <- paste0("Lineweaver-Burk Model (Double Reciprocal Plot):\n")
    stats_text <- paste0(stats_text, "Equation: 1/V = ", round(slope, 4), " × (1/[S]) + ", round(intercept, 4), "\n")
    stats_text <- paste0(stats_text, "Linear form: y = ", round(slope, 4), "x + ", round(intercept, 4), "\n")
    stats_text <- paste0(stats_text, "Slope = Km/Vmax = ", round(slope, 4), "\n")
    stats_text <- paste0(stats_text, "Y-intercept = 1/Vmax = ", round(intercept, 4), "\n")
    stats_text <- paste0(stats_text, "X-intercept = -1/Km = ", round(lb_data$x_intercept, 4), "\n")
    stats_text <- paste0(stats_text, "Vmax (from L-B) = ", round(lb_data$Vmax, 4), " ", rate_units, "\n")
    stats_text <- paste0(stats_text, "Km (from L-B) = ", round(lb_data$Km, 4), " ", conc_units, "\n")
    stats_text <- paste0(stats_text, "R² = ", round(lb_data$r_squared, 4))
    
    return(stats_text)
  })
  
  # Render Michaelis-Menten plot for Tab 2
  output$michaelis_menten_plot <- renderPlot({
    mm_data <- michaelis_menten_fit()
    
    if (is.null(mm_data)) {
      return(ggplot() + 
               geom_text(aes(x = 0.5, y = 0.5, label = "Enter at least 3 complete\ndata points to generate plot"), 
                         size = 5) +
               theme_void() +
               xlim(0, 1) + ylim(0, 1))
    }
    
    tryCatch({
      # Create the plot with units in axis labels
      p <- ggplot(mm_data$data, aes(x = S, y = V)) +
        geom_point(color = "black", size = 3) +
        geom_line(data = mm_data$pred_data, aes(x = S, y = V), color = "blue", size = 1) +
        labs(x = paste0("Substrate Concentration [S] (", conc_units, ")"), 
             y = paste0("Initial Velocity (V) (", rate_units, ")"), 
             title = "Michaelis-Menten Plot") +
        theme_minimal() +
        theme(axis.line.x = element_line(color = "black", size = 0.5),
              axis.line.y = element_line(color = "black", size = 0.5),
              plot.title = element_text(hjust = 0.5, size = 14),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10)) +
        theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
              axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
        scale_y_continuous(limits = c(0, NA))
      
      # Add Vmax and Km reference lines with units in annotations
      p <- p + 
        geom_hline(yintercept = mm_data$Vmax, linetype = "dashed", color = "red", alpha = 0.7) +
        geom_vline(xintercept = mm_data$Km, linetype = "dashed", color = "red", alpha = 0.7) +
        geom_hline(yintercept = mm_data$Vmax/2, linetype = "dashed", color = "orange", alpha = 0.7) +
        annotate("text", x = max(mm_data$data$S) * 0.7, y = mm_data$Vmax * 1.05, 
                 label = paste("Vmax =", round(mm_data$Vmax, 3)), color = "red", size = 3) +
        annotate("text", x = mm_data$Km * 1.1, y = max(mm_data$data$V) * 0.9, 
                 label = paste("Km =", round(mm_data$Km, 3)), color = "red", size = 3) +
        annotate("text", x = max(mm_data$data$S) * 0.7, y = (mm_data$Vmax/2) * 1.1, 
                 label = paste("Vmax/2 =", round(mm_data$Vmax/2, 3)), color = "orange", size = 3)
      
      return(p)
      
    }, error = function(e) {
      return(ggplot() + 
               geom_text(aes(x = 0.5, y = 0.5, label = "Error generating plot.\nCheck your data."), 
                         size = 5, color = "red") +
               theme_void() +
               xlim(0, 1) + ylim(0, 1))
    })
  })
  
  # Render Lineweaver-Burk plot
  output$lineweaver_burk_plot <- renderPlot({
    lb_data <- lineweaver_burk_fit()
    
    if (is.null(lb_data)) {
      return(ggplot() + 
               geom_text(aes(x = 0.5, y = 0.5, label = "Enter at least 3 complete\ndata points to generate plot"), 
                         size = 5) +
               theme_void() +
               xlim(0, 1) + ylim(0, 1))
    }
    
    tryCatch({
      # Calculate appropriate axis limits to show both intercepts
      x_min_data <- min(lb_data$data$inv_S)
      x_max_data <- max(lb_data$data$inv_S)
      x_intercept <- lb_data$x_intercept
      y_intercept <- lb_data$y_intercept
      
      # Extend x-axis to show x-intercept with padding
      x_min_plot <- min(x_min_data, x_intercept * 1.2)
      x_max_plot <- max(x_max_data, x_max_data * 1.1)
      
      # Extend y-axis to show y=0 and provide good view of y-intercept
      y_min_plot <- min(0, min(lb_data$data$inv_V) * 0.1)
      y_max_plot <- max(lb_data$data$inv_V) * 1.1
      
      # Create the Lineweaver-Burk plot
      p <- ggplot(lb_data$data, aes(x = inv_S, y = inv_V)) +
        geom_point(color = "black", size = 3) +
        geom_line(data = lb_data$pred_data, aes(x = inv_S, y = inv_V), color = "blue", size = 1) +
        labs(x = paste0("1/[S] (1/", conc_units, ")"), 
             y = paste0("1/V (1/", rate_units, ")"), 
             title = "Lineweaver-Burk Plot") +
        theme_minimal() +
        theme(axis.line.x = element_line(color = "black", size = 0.5),
              axis.line.y = element_line(color = "black", size = 0.5),
              plot.title = element_text(hjust = 0.5, size = 14),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10)) +
        theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
              axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
        scale_x_continuous(limits = c(x_min_plot, x_max_plot)) +
        scale_y_continuous(limits = c(y_min_plot, y_max_plot))
      
      # Add y-intercept line and annotation
      if (lb_data$y_intercept > y_min_plot && lb_data$y_intercept < y_max_plot) {
        p <- p + geom_hline(yintercept = lb_data$y_intercept, linetype = "dashed", color = "red", alpha = 0.7)
        p <- p + annotate("text", x = x_max_plot * 0.7, y = lb_data$y_intercept * 1.1, 
                          label = paste("1/Vmax =", round(lb_data$y_intercept, 3)), color = "red", size = 3)
      }
      
      # Add x-intercept line and annotation (now always visible)
      p <- p + geom_vline(xintercept = lb_data$x_intercept, linetype = "dashed", color = "red", alpha = 0.7)
      p <- p + annotate("text", x = lb_data$x_intercept, y = y_max_plot * 0.9, 
                        label = paste("-1/Km =", round(lb_data$x_intercept, 3)), color = "red", size = 3,
                        hjust = ifelse(lb_data$x_intercept < 0, 1.1, -0.1))
      
      return(p)
      
    }, error = function(e) {
      return(ggplot() + 
               geom_text(aes(x = 0.5, y = 0.5, label = "Error generating plot.\nCheck your data."), 
                         size = 5, color = "red") +
               theme_void() +
               xlim(0, 1) + ylim(0, 1))
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)