#' bca UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import shiny
mod_bca_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "BCA Protein Assay",

    sidebarLayout(
      # sidebarPanel ----
      sidebarPanel = sidebarPanel(
        width = 5,

        shinyMatrix::matrixInput(
          inputId = ns("mat_std_curve"),
          label = "Values for standard curve",
          value = matrix(
            # c(100, 200, 400, 600, 800, 0, 0, 0, 0, 0),
            c(100, 200, 400, 600, 800, 26442, 40792, 92383, 135403, 177405),
            nrow = 5,
            ncol = 2,
            dimnames = list(NULL, c("BSA Level", "Signal"))
          ),
          inputClass = "numeric",
          rows = list(
            extend = TRUE,
            names = FALSE
          ),
          cols = list(
            extend = TRUE,
            names = TRUE,
            editableNames = TRUE
          )
        ),

        shinyMatrix::matrixInput(
          inputId = ns("mat_sample_signal"),
          label = "Signal absorbance values",
          value = matrix(
            # rep("0", 3),
            # nrow = 1,
            # ncol = 3,
            # dimnames = list(NULL, c("Subtrahend sample", paste("Sample", 1:2)))
            c(3239, 3351, 3305, 94613, 93828, 93380, 26388, 26840, 27044, 33545, 34215, 34566),
            nrow = 3,
            ncol = 4,
            dimnames = list(NULL, c("Subtrahend sample", paste("Sample", 1:3)))
          ),
          inputClass = "numeric",
          rows = list(
            extend = TRUE,
            names = FALSE
          ),
          cols = list(
            extend = TRUE,
            names = TRUE,
            editableNames = TRUE
          )
        ),

        selectizeInput(
          ns("wildtype_selection"),
          choices = NULL,
          label = "Control sample",
          options = list(placeholder = "Select control sample")
        ),

        actionButton(
          inputId = ns("get_prot_level"),
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

        tabsetPanel(
          type = "tabs",
          tabPanel(
            title = "Signal absorbance plot",
            plotOutput(ns("plot_bar_signal"))
          ),
          tabPanel(
            title = "Fitting results",
            plotOutput(ns("plot_std_curve")),
            verbatimTextOutput(ns("std_curve_fit_summary"))
          ),
          tabPanel(
            title = "Protein concentration level plot",
            plotOutput(ns("plot_bar_prot"))
          ),
          tabPanel(
            title = "Result table",
            DT::DTOutput(ns("table_prot"))
          )
        )
      )
    )
  )
}

#' bca Server Functions
#'
#' @noRd
#' @importFrom rlang .data
mod_bca_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    react_vals <- reactiveValues()

    # Calculations ----
    observeEvent(
      input$get_prot_level,
      {
        # Read inputs
        mat_std_curve <- input$mat_std_curve
        mat_sample_signal <- input$mat_sample_signal

        # Process data
        std_gfp <- mat_std_curve[, 1]
        std_fluorescence <- mat_std_curve[, 2:ncol(mat_std_curve)]

        if (length(std_fluorescence) > length(std_gfp)) {
          std_fluorescence <- apply(std_fluorescence, 2 , as.numeric)
          std_fluorescence <- rowMeans(std_fluorescence)
        }

        # Calculations
        list_std_curve <- get_std_curve_value(std_fluorescence, std_gfp)
        df_tidied <- get_fluorescence_input(mat_sample_signal)
        df_with_pred_gfp <- predict_gfp_from_fluorescence(df_tidied, list_std_curve$std_curve_fit)
        df_with_pred_gfp_kg <- df_with_pred_gfp %>%
          `colnames<-`(c("Sample", "Signal absorbance", "Protein conc. (ng)"))

        # Reactive values
        react_vals$list_std_curve <- list_std_curve
        react_vals$df_tidied <- df_tidied
        react_vals$df_with_pred_gfp <- df_with_pred_gfp
        react_vals$df_with_pred_gfp_kg <- df_with_pred_gfp_kg
      }
    )

    # Update wildtype_selection ----
    observeEvent(
      input$mat_sample_signal,
      {
        updateSelectizeInput(
          session = session,
          inputId = "wildtype_selection",
          choices = colnames(input$mat_sample_signal),
          server = TRUE
        )
      }
    )

    # Debugging panel ----
    output$input_panel_output <- renderPrint({
      if (input$get_prot_level == 0) {
        return(print("Input your values"))
      }

      input$get_prot_level
      isolate({
        mat_std_curve <- input$mat_std_curve
        mat_sample_signal <- input$mat_sample_signal
      })

      vec_a <- mat_std_curve[, 1]
      # vec_a <- vec_a[1:(length(vec_a) - 1)]
      vec_b <- mat_std_curve[, 2]
      # vec_b <- vec_b[1:(length(vec_b) - 1)]

      # utils::str(list(
      #   "Matrix A" = mat_std_curve,
      #   "Matrix B" = mat_sample_signal,
      #   "Fit summary" = summary(get_std_curve_value(vec_a, vec_b)$std_curve_fit)
      #   # "Vector A" = vec_a,
      #   # "Vector B" = vec_b,
      #   # "First function" = get_std_curve_value(vec_a, vec_b)
      # ))
      summary(get_std_curve_value(vec_a, vec_b)$std_curve_fit)
    })

    # Plots ----
    output$plot_std_curve <- renderPlot(
      {
        if (input$get_prot_level == 0) {
          return(NULL)
        }

        input$get_prot_level
        isolate({
          list_std_curve <- react_vals$list_std_curve
          df_tidied <- react_vals$df_tidied
          df_with_pred_gfp <- react_vals$df_with_pred_gfp
        })

        # Calculations
        gg_plot <- plot_std_curve_and_pred(
          list_std_curve$std_curve_df,
          df_with_pred_gfp,
          list_std_curve$std_curve_fit,
          xlab = "Protein concentration (ng)",
          ylab = "Signal absorbance"
        )

        return(gg_plot)
      },
      res = 96
    )

    output$plot_bar_signal <- renderPlot(
      {
        if (input$get_prot_level == 0) {
          return(NULL)
        }

        input$get_prot_level
        isolate({
          mat_sample_signal <- input$mat_sample_signal
        })

        # Plot
        gg_plot <- plot_bar_fluorescence(
          mat_sample_signal,
          ylab = "Signal absorbance"
        )

        return(gg_plot)
      },
      res = 96
    )

    output$plot_bar_prot <- renderPlot(
      {
        if (input$get_prot_level == 0) {
          return(NULL)
        }

        input$get_prot_level
        isolate({
          df_with_pred_gfp_kg <- react_vals$df_with_pred_gfp_kg
          wildtype_sample <- input$wildtype_selection
        })

        # Plot
        gg_plot <- plot_bar_prot(
          df_with_pred_gfp_kg,
          wildtype_sample,
          var_to_plot = "Protein conc. (ng)",
          ylab = "Protein concentration (ng)",
          percent_lab = "Percentage of protein concentration"
        )

        return(gg_plot)
      },
      res = 96
    )

    # Table ----
    output$table_prot <- DT::renderDT({
      if (input$get_prot_level == 0) {
        return(DT::datatable(NULL, style = "bootstrap4"))
      }

      input$get_prot_level
      isolate({
        df_with_pred_gfp_kg <- react_vals$df_with_pred_gfp_kg
      })

      table <- datatable_bca(df_with_pred_gfp_kg)

      return(table)
    })

    # Fit summary ----
    output$std_curve_fit_summary <- renderPrint({
      if (input$get_prot_level == 0) {
        return(print(""))
      }

      input$get_prot_level
      isolate({
        std_curve_fit <- react_vals$list_std_curve$std_curve_fit
      })

      return(summary(std_curve_fit))
    })
  })
}

## To be copied in the UI
# mod_bca_ui("bca_ui_1")

## To be copied in the server
# mod_bca_server("bca_ui_1")
