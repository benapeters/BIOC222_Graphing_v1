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
    ),
    fluidRow(
      column(6,"Move the slider so that only the initial linear portion of the graph is used for the line of best fit")
    ),
    sliderInput("slider_id", "", min = 20, max = 180, step = 20, value = 20),
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
    Assay1 = rep(NA, 9),
    Assay2 = rep(NA, 9),
    stringsAsFactors = FALSE
  ))
  
  #The slider text
  output$slider_label <- renderText({
    "Move the slider so that only the initial linear portion of the graph is used for the line of best fit"
  })
  
  # Convert Assay1 and Assay2 to numeric (if possible)
  observe({
    data2$df$Assay1 <- as.numeric(data2$df$Assay1)
    data2$df$Assay2 <- as.numeric(data2$df$Assay2)
  })
  # Render the data table for Tab 2
  output$table2 <- renderRHandsontable({
    rhandsontable(data2$df, rowHeaders = FALSE) %>% 
      hot_col("Time", type = "numeric", strict = TRUE, allowInvalid = FALSE, readOnly = TRUE, format = 0) %>%
      hot_col("Assay1", type = "numeric", strict = TRUE, allowInvalid = FALSE, format = 0.000) %>%
      hot_col("Assay2", type = "numeric", strict = TRUE, allowInvalid = FALSE, format = 0.000)
  })
  # Update dataframe after inputs
  observe({
    if (!is.null(input$table2)) {
      data2$df <- hot_to_r(input$table2)
    }
  })
  
  # Render the progress Curve scatter plot for Tab 2
  output$progressCurve <- renderPlot({
    req(nrow(data2$df) > 0)
    
    #writes the slider input to the slider value
    slider_value <- input$slider_id
    
    # Subset the data based on the slider value
    subset_data <- data2$df[data2$df$Time <= slider_value, ]
    
    tryCatch({
      output$annotated_points <- renderTable({
        subset_data <- subset_data
        subset_data <- as.data.frame(subset_data)
        
        if (!any(is.na(subset_data))) {
          # Fit the Michaelis-Menten equation to the data
          fit1 <- lm(Assay1 ~ Time, data = subset_data)
          fit2 <- lm(Assay2 ~ Time, data = subset_data)
          
          # Create a data frame with the annotated points
          annotated_points <- data.frame(
            Assay = rep(c("Assay1", "Assay2"), each = 2),
            Time = c(20, 80, 20, 80),
            Value = c(predict(fit1, newdata = data.frame(Time = 20)),
                      predict(fit1, newdata = data.frame(Time = 80)),
                      predict(fit2, newdata = data.frame(Time = 20)),
                      predict(fit2, newdata = data.frame(Time = 80)))
          )
          
          
          # Return the annotated points
          return(annotated_points)}
      }, error = function(e){"The table will work once you have entered data"})
    })
    
    
    
    tryCatch({
      ggplot(data2$df, aes(x = Time)) +
        geom_line(aes(y = Assay1), color = "grey", size = 2) +
        geom_line(aes(y = Assay2), color = "pink", size = 2) +
        geom_smooth(data = subset_data, aes(y = Assay1), color = "black", method = "lm", se = FALSE, fullrange = TRUE, linetype = "dotted") +
        geom_smooth(data = subset_data, aes(y = Assay2), color = "red", method = "lm", se = FALSE, fullrange = TRUE, linetype = "dotted") +
        geom_point(aes(y = Assay1), color = "black", size = 3) +
        geom_point(aes(y = Assay2), color = "red", size = 3) +
        labs(x = "Time (seconds)", y = "Absorbance", title = "Alcohol Dehydrogenase Assay") +
        theme_minimal() +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max(data2$df$Assay1, data2$df$Assay2, na.rm = TRUE) * 1.1)) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 200), breaks = seq(0, max(data2$df$Time, na.rm = TRUE), by = 60), 
                           minor_breaks = seq(0, max(data2$df$Time, na.rm = TRUE), by = 20)) +
        theme(legend.position = "bottomright") +
        theme(axis.line.x = element_line(color = "black", size = 1),
              axis.line.y = element_line(color = "black", size = 1),
              plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values = c("black", "red"), labels = c("Assay1", "Assay2")) +
        theme(panel.grid.major = element_line(colour = "grey", linetype = "solid"),
              axis.line = element_line(colour = "black", size = 1, linetype = "solid"))
      
      
    }, error = function(e){"the plot will appear once you have entered your data"})
  })
   
    
}

# Run the application 
shinyApp(ui = ui, server = server)
