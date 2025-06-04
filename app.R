library(shiny)
library(ggplot2)
library(rhandsontable)


#UI section
ui <- fluidPage(
  
  # Application title
  titlePanel("BIOC222 Enzyme Kinetics Calculator"),
  
  tabsetPanel(
    
    tabPanel(
      "Initial rate calculation",
      fluidRow(
        column(6,rHandsontableOutput("table2")),
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
      sliderInput("slider_id", "", min = 20, max = 180, step = 20, value =c(20,60)),
      plotOutput("progressCurve")
    ),
    
    tabPanel(
      "Km and Vmax Calculation"
    ),
    tabPanel(
      "Enzyme Inhibition"
    )
  ),
  
  
  
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Define the initial data frame for Tab 2
  data2 <- reactiveValues(df = data.frame(
    Time = seq(20, 180, by = 20),
    Assay = rep(NA, 9),
    stringsAsFactors = FALSE
  ))
  
  #The slider text
  output$slider_label <- renderText({
    "Move the slider so that only the initial linear portion of the graph is used for the line of best fit"
  })
  
  # Convert Assay to numeric (if possible)
  observe({
    data2$df$Assay <- as.numeric(data2$df$Assay)
  })
  
  # Render the data table for Tab 2
  output$table2 <- renderRHandsontable({
    rhandsontable(data2$df, rowHeaders = FALSE) %>% 
      hot_col("Time", type = "numeric", strict = TRUE, allowInvalid = FALSE, readOnly = TRUE, format = 0) %>%
      hot_col("Assay", type = "numeric", strict = TRUE, allowInvalid = FALSE, format = 0.000)
  })
  
  # Update dataframe after inputs
  observe({
    if (!is.null(input$table2)) {
      data2$df <- hot_to_r(input$table2)
    }
  })
  
  # Calculate initial rates and model statistics
  initial_rates_data <- reactive({
    req(nrow(data2$df) > 0)
    
    slider_value_min <- input$slider_id[1]
    slider_value_max <- input$slider_id[2]
    
    # Subset the data based on the slider value
    subset_data <- data2$df[data2$df$Time >= slider_value_min & data2$df$Time <= slider_value_max, ]
    
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
  
  # Display initial rates table
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
  
  # Display model statistics
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
  
  # Render the progress Curve scatter plot for Tab 2
  output$progressCurve <- renderPlot({
    req(nrow(data2$df) > 0)
    
    #writes the slider input to the slider value
    slider_value_min <- input$slider_id[1]
    slider_value_max <- input$slider_id[2]
    
    # Subset the data based on the slider value
    subset_data <- data2$df[data2$df$Time >= slider_value_min & data2$df$Time <= slider_value_max, ]
    
    tryCatch({
      ggplot(data2$df, aes(x = Time)) +
        geom_line(aes(y = Assay), color = "grey", size = 2) +
        geom_smooth(data = subset_data, aes(y = Assay), color = "black", method = "lm", se = FALSE, fullrange = TRUE, linetype = "dotted") +
        geom_point(aes(y = Assay), color = "black", size = 3) +
        labs(x = "Time (seconds)", y = "Absorbance", title = "Initial rate of reaction") +
        theme_minimal() +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max(data2$df$Assay, na.rm = TRUE) * 1.1)) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 200), breaks = seq(0, max(data2$df$Time, na.rm = TRUE), by = 60), 
                           minor_breaks = seq(0, max(data2$df$Time, na.rm = TRUE), by = 20)) +
        theme(axis.line.x = element_line(color = "black", size = 1),
              axis.line.y = element_line(color = "black", size = 1),
              plot.title = element_text(hjust = 0.5)) +
        theme(panel.grid.major = element_line(colour = "grey", linetype = "solid"),
              axis.line = element_line(colour = "black", size = 1, linetype = "solid"))
      
      
    }, error = function(e){"The plot will appear once you have entered your data"})
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)