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
      # sidebarPanel ----
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
          label = "Calculate",
          class = "btn-primary",
          icon = icon("calculator")
        )
      ),
      # mainPanel ----
      mainPanel = mainPanel(
        width = 7,

        h4("Results"),
        # verbatimTextOutput(ns("input_panel_output")),

        plotOutput(ns("plot_std_curve")),
        plotOutput(ns("plot_bar_fluorescence")),
        plotOutput(ns("plot_bar_gfp")),
        DT::DTOutput(ns("table_gfp_kg"))
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

    react_vals <- reactiveValues()

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

    output$plot_std_curve <- renderPlot({
      if (input$get_gfp_level == 0) {
        return(NULL)
      }

      input$get_gfp_level
      isolate({
        mat_std_curve <- input$mat_std_curve
        mat_sample_fluorescence <- input$mat_sample_fluorescence
      })

      # Process data
      std_gfp <- mat_std_curve[,1]
      std_gfp <- std_gfp[1:(length(std_gfp) - 1)]
      std_fluorescence <- mat_std_curve[,2]
      std_fluorescence <- std_fluorescence[1:(length(std_fluorescence) - 1)]

      extra_row <- nrow(mat_sample_fluorescence)
      extra_col <- ncol(mat_sample_fluorescence)
      mat_sample_fluorescence <- mat_sample_fluorescence[-extra_row, -extra_col]
      # browser()

      # Calculations
      list_std_curve <- get_std_curve_value(std_fluorescence, std_gfp)
      df_tidied <- get_fluorescence_input(mat_sample_fluorescence)
      df_with_pred_gfp <- predict_gfp_from_fluorescence(df_tidied, list_std_curve$std_curve_fit)

      # browser()

      gg_plot <- plot_std_curve_and_pred(list_std_curve$std_curve_df, df_with_pred_gfp, list_std_curve$std_curve_fit)

      # Reactive values
      react_vals$df_with_pred_gfp <- df_with_pred_gfp

      return(gg_plot)
    }, res = 96)

    output$plot_bar_fluorescence <- renderPlot({
      if (input$get_gfp_level == 0) {
        return(NULL)
      }

      input$get_gfp_level
      isolate({
        mat_sample_fluorescence <- input$mat_sample_fluorescence
      })

      # Process data
      extra_row <- nrow(mat_sample_fluorescence)
      extra_col <- ncol(mat_sample_fluorescence)
      mat_sample_fluorescence <- mat_sample_fluorescence[-extra_row, -extra_col]

      # Plot
      gg_plot <- plot_bar_fluorescence(mat_sample_fluorescence)

      return(gg_plot)
    }, res = 96)

    output$plot_bar_gfp <- renderPlot({
      if (input$get_gfp_level == 0) {
        return(NULL)
      }

      input$get_gfp_level
      isolate({
        df_with_pred_gfp <- react_vals$df_with_pred_gfp
      })

      # Process data
      df_with_pred_gfp_kg <- df_with_pred_gfp %>%
        dplyr::mutate(
          gfp_final = (gfp / 200) / 1000
        ) %>%
        `colnames<-`(c("Sample", "Fluorescence", "GFP (ng)", "GFP (g/kg)"))

      # Reactive values
      react_vals$df_with_pred_gfp_kg <- df_with_pred_gfp_kg

      # Plot
      gg_plot <- plot_bar_gfp(df_with_pred_gfp_kg)

      return(gg_plot)
    }, res = 96)

    output$table_gfp_kg <- DT::renderDT({
      if (input$get_gfp_level == 0) {
        return(DT::datatable(NULL, style = "bootstrap4"))
      }

      input$get_gfp_level
      isolate({
        df_with_pred_gfp_kg <- react_vals$df_with_pred_gfp_kg
      })

      table <- DT::datatable(
        data = df_with_pred_gfp_kg,
        style = "bootstrap4",
        rownames = FALSE,
        selection = "none",
        extensions = "Buttons",
        options = list(
          pageLength = 10,
          filter = FALSE,
          lengthChange = FALSE,
          scrollX = TRUE,
          dom = "tB",
      # <'row'<'col-sm-12'tr>>
      # <'row'<'col-sm-12 col-md-7'pB><'col-sm-12 col-md-5 text-right'i>>
      # ",
          buttons = list(
            list(
              extend = "csv",
              filename = "data"
            ),
            list(
              extend = "excel",
              filename = "data"
            )
          )
        )
      )

      return(table)
    })

  })
}

## To be copied in the UI
# mod_gfp_ui("gfp_ui_1")

## To be copied in the server
# mod_gfp_server("gfp_ui_1")
