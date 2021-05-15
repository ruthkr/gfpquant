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
          value = matrix(
            # 0,
            c(100, 200, 400, 600, 800, 26442, 40792, 92383, 135403, 177405),
            # byrow = TRUE,
            nrow = 5,
            ncol = 2,
            dimnames = list(NULL, c("GFP Level", "Fluorescence"))
          ),
          inputClass = "numeric",
          rows = list(
            extend = TRUE,
            names = FALSE
          ),
          cols = list(
            names = TRUE
          )
        ),

        shinyMatrix::matrixInput(
          inputId = ns("mat_sample_fluorescence"),
          label = "Samples' fluorescence values",
          value = matrix(
            # 0,
            c(3239, 3351, 3305, 94613, 93828, 93380, 26388, 26840, 27044, 33545, 34215, 34566),
            # byrow = TRUE,
            nrow = 3,
            ncol = 4,
            dimnames = list(NULL, c("Uninfiltrated Leaves", paste("Sample", 1:3)))
          ),
          inputClass = "numeric",
          rows = list(
            extend = TRUE,
            names = FALSE
          ),
          cols = list(
            names = TRUE,
            editableNames = TRUE,
            extend = TRUE
          )
        ),

        actionButton(
          inputId = ns("get_gfp_level"),
          label = "Get GFP",
          class = "btn-primary",
          icon = icon("calculator")
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
      if (input$get_gfp_level == 0) {
        return(print("Input your values"))
      }

      input$get_gfp_level
      isolate({
        mat_std_curve <- input$mat_std_curve
        mat_sample_fluorescence <- input$mat_sample_fluorescence
      })

      vec_a <- mat_std_curve[,1]
      vec_a <- vec_a[1:(length(vec_a) - 1)]
      vec_b <- mat_std_curve[,2]
      vec_b <- vec_b[1:(length(vec_b) - 1)]

      utils::str(list(
        "Matrix A" = mat_std_curve,
        "Matrix B" = mat_sample_fluorescence,
        "Vector A" = vec_a,
        "Vector B" = vec_b,
        "First function" = get_std_curve_value(vec_a, vec_b)
      ))
    })
  })
}

## To be copied in the UI
# mod_gfp_ui("gfp_ui_1")

## To be copied in the server
# mod_gfp_server("gfp_ui_1")
