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
      sidebarPanel = shiny::sidebarPanel(),
      mainPanel = mainPanel()
    )
  )
}

#' gfp Server Functions
#'
#' @noRd
mod_gfp_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  })
}

## To be copied in the UI
# mod_gfp_ui("gfp_ui_1")

## To be copied in the server
# mod_gfp_server("gfp_ui_1")
