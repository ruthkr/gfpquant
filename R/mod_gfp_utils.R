#' Plot GFP
#'
#' @param df_with_pred_gfp Data frame with predicted GFP.
#' @param wildtype_sample Name of wildtype/control sample.
#' @param var_to_plot Variable name to plot.
#'
#' @importFrom rlang .data
plot_bar_gfp <- function(df_with_pred_gfp, wildtype_sample, var_to_plot = "GFP (g/kg)") {
  df_with_pred_gfp <- df_with_pred_gfp %>%
    dplyr::rename(wanted_var = var_to_plot)

  total_gfp <- df_with_pred_gfp %>%
    dplyr::filter(.data$Sample == wildtype_sample) %>%
    dplyr::pull(.data$wanted_var)

  gg_gfp <- df_with_pred_gfp %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$Sample, y = .data$wanted_var) +
    ggplot2::geom_bar(
      stat = "identity",
      fill = "mediumseagreen",
      width = 0.5,
      alpha = 1
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 0.01, 0.001),
      sec.axis = ggplot2::sec_axis(
        ~ (. / total_gfp),
        name = "Percentage of GFP",
        labels = scales::percent
      )
    ) +
    ggplot2::labs(y = "GFP (g/kg)")
  # ggplot2::theme(
  #   axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust=1),
  #   axis.title.x= ggplot2::element_blank()
  # )

  return(gg_gfp)
}

#' @noRd
datatable_gfp <- function(data) {
  table <- DT::datatable(
    data = data,
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
  ) %>%
    DT::formatRound(
      columns = c("Fluorescence", "GFP (ng)"),
      digits = 2
    ) %>%
    DT::formatSignif(
      columns = c("GFP (g/kg)"),
      digits = 3
    )

  return(table)
}
