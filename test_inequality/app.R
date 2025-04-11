#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggplot2)
library(scales)

my_sld <- function(...) {
  sld <- sliderInput(...)
  sld$children[[2]]$attribs$`data-immediate` <- "true"
  sld
}

R1 <- seq(from = .3, to = 1, length.out = 50)
R2 <- seq(from = .01, to = 1, length.out = 50)

R_grid <- expand.grid(R1, R2)
names(R_grid) <- c("R1", "R2")



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Visualising inequality"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        my_sld("mu_MD", "mu MD:", min = 0, max = 2.5, value = 1),
        my_sld("tau_MD", "tau MD:", min = 0, max = 1, value = .2),
        my_sld("mu_sigma2x", "mu sigma2_X:", min = 0, max = 2.5, value = 2),
        my_sld("tau_sigma2x", "tau sigma2_X:", min = 0, max = 1, value = .3)
      ),
      

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
      
      ineq_func <- function(R_grid, tau_MD, mu_MD, mu_sigma2x, tau_sigma2x){
        
        val <- ((tau_MD^2 )* ((1/R_grid$R1) - 1)) + (.25 * (mu_MD^2) * ((tau_sigma2x/mu_sigma2x)^2) * ((R_grid$R2 / (R_grid$R1^3)) - 1))
        
        data.frame(R1 = R_grid$R1,
                   R2 = R_grid$R2,
                   val = val > 0)
      }
      
      input_df <- ineq_func(R_grid, 
                            tau_MD = input$tau_MD, 
                            mu_MD = input$mu_MD, 
                            mu_sigma2x = input$mu_sigma2x, 
                            tau_sigma2x = input$tau_sigma2x)
      
      ggplot(input_df) +
        geom_raster(aes(x = R1, y = R2, fill = val)) +
        labs(title = "Does the inequality hold?",
             fill = "Inequality\nholds") 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
