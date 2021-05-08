#' gfp UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import shiny
mod_gfp_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Quantify GFP",

    sidebarLayout(
      sidebarPanel = sidebarPanel(
        width = 5,
        shinyMatrix::matrixInput(
          inputId = ns("mat_std_curve"),
          label = "Values for standard curve",
          value = matrix(0, nrow = 3, ncol = 2) %>%
            `colnames<-`(c("GFP Level", "Fluorescence")),
          inputClass = "numeric",
          rows = list(
            extend = TRUE
          ),
          cols = list(
            names = TRUE
          )
        ),

        shinyMatrix::matrixInput(
          inputId = ns("mat_sample_fluorescence"),
          label = "Samples' fluorescence values",
          value = matrix(0, nrow = 3, ncol = 3) %>%
            `colnames<-`(paste("Sample", 1:3)),
          inputClass = "numeric",
          rows = list(
            extend = TRUE
          ),
          cols = list(
            names = TRUE,
            editableNames = TRUE,
            extend = TRUE
          )
        )

      ),
      mainPanel = mainPanel(
        width = 7,

        h4("Results"),
        verbatimTextOutput(ns("input_panel_output"))
      )
    )
  )
}

#' gfp Server Functions
#'
#' @noRd
mod_gfp_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$input_panel_output <- renderPrint({
      utils::str(list(
        "Matrix A" = input$mat_std_curve,
        "Matrix B" = input$mat_sample_fluorescence
      ))
    })
  })
}

## To be copied in the UI
# mod_gfp_ui("gfp_ui_1")

## To be copied in the server
# mod_gfp_server("gfp_ui_1")
