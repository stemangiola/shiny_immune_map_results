library(shiny)
library(gganatogram)

prop_data <- readRDS("cell_proportion_for_shiny_app.rds")

ui <- navbarPage(
    "Immune Proportion Map",
    tabPanel("Full Data",
        sidebarLayout(
            sidebarPanel(
                h3("Settings"),
                p("placeholder")
            ),
            mainPanel(
                h3("placeholder"),
                p("content")
            )
        )
    ),
    tabPanel("Immune Proportions",
        sidebarLayout(
            sidebarPanel(
                h3("Settings"),
                p("placeholder")
            ),
            mainPanel(
                h3("placeholder"),
                p("content")
            )
        )
    ),
    tabPanel("Gene Expression",
        sidebarLayout(
            sidebarPanel(
                h3("Settings"),
                p("placeholder")
            ),
            mainPanel(
                h3("placeholder"),
                p("content")
            )
        )
    )
)

server <- function(input, output) {
    
}

shinyApp(ui = ui, server = server)