#' Function to calculate standard curve
#'
#' @param std_fluorescence Fluorescence for the standard curve.
#' @param std_gfp GFP for the standard curve.
get_std_curve_value <- function(std_fluorescence, std_gfp) {
  df <- data.frame(
    fluorescence = as.numeric(std_fluorescence),
    gfp = as.numeric(std_gfp)
  )

  fit <- stats::lm(fluorescence ~ gfp, data = df)

  list_std_curve <- list(
    std_curve_df = df,
    std_curve_fit = fit
  )

  return(list_std_curve)
}

#' Function to get the value, given input data
#'
#' @param mdata Sample fluorescence matrix.
#' @importFrom rlang .data
get_fluorescence_input <- function(mdata) {
  df_input <- as.data.frame(mdata, row.names = FALSE) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = as.numeric
      )
    )

  colnames(df_input)[1] <- "uninfiltrated_leaves"

  # Get new df with the fluorescence value subtracted from infiltrated leave
  df_tidied <- df_input %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = ~ mean(.x, na.rm = TRUE)
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = -.data$uninfiltrated_leaves,
        .fns = ~ .x - .data$uninfiltrated_leaves
      )
    ) %>%
    dplyr::select(-c("uninfiltrated_leaves")) %>%
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "sample",
      values_to = "fluorescence"
    )

  return(df_tidied)
}

#' Function to get GFP value if fluorescence known
#'
#' @param df_tidied Tidied fata frame.
#' @param fit Fit object.
#'
#' @importFrom rlang .data
predict_gfp_from_fluorescence <- function(df_tidied, fit) {
  intercept <- fit$coefficients[[1]]
  slope <- fit$coefficients[[2]]

  df_with_pred_gfp <- df_tidied %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      gfp = (.data$fluorescence - intercept) / slope
    )

  return(df_with_pred_gfp)
}

#' Plot standard curve and predictions
#'
#' @param df_std_curve Standard curve data frame.
#' @param df_with_pred_gfp Data frame with predicted GFP.
#' @param fit Fit object.
#'
#' @importFrom rlang .data
plot_std_curve_and_pred <- function(df_std_curve, df_with_pred_gfp, fit) {
  # Plot std curve
  gg <- df_std_curve %>%
    ggplot2::ggplot() +
    ggplot2::aes(
      x = .data$gfp,
      y = .data$fluorescence
    ) +
    # geom_smooth(method = "lm") +
    ggplot2::geom_point(
      color = "#8856a7"
    ) +
    ggplot2::geom_abline(
      ggplot2::aes(
        slope = stats::coef(fit)[[2]],
        intercept = stats::coef(fit)[[1]]
      ),
      linetype = "solid",
      color = "#9ebcda"
    ) +
    #  Vertical lines
    ggplot2::geom_segment(
      data = df_with_pred_gfp,
      ggplot2::aes(
        x = .data$gfp, xend = .data$gfp,
        y = 0, yend = .data$fluorescence
      ),
      linetype = "dashed"
    ) +
    # Horizontal lines
    ggplot2::geom_segment(
      data = df_with_pred_gfp,
      ggplot2::aes(
        x = 0, xend = .data$gfp,
        y = .data$fluorescence, yend = .data$fluorescence
      ),
      linetype = "dashed"
    ) +
    # Predictions
    ggplot2::geom_point(data = df_with_pred_gfp, shape = 2, size = 2) +
    # Labels
    ggrepel::geom_label_repel(
      data = df_with_pred_gfp,
      ggplot2::aes(label = .data$sample),
      seed = 1,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = "GFP (ng)",
      y = "Flouresence"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.01), add = c(2, 0))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.01), add = c(2, 0)))

  return(gg)
}

#' Plot fluorescence
#'
#' @param mat_sample_fluorescence Input data.
#'
#' @importFrom rlang .data
plot_bar_fluorescence <- function(mat_sample_fluorescence) {
  df_tidied <- get_fluorescence_input(mat_sample_fluorescence)

  df <- as.data.frame(mat_sample_fluorescence) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = as.numeric
      )
    )

  colnames(df)[1] <- "uninfiltrated_leaves"
  df <- df %>%
    dplyr::mutate(
      dplyr::across(
        .cols = -.data$uninfiltrated_leaves,
        .fns = ~ .x - .data$uninfiltrated_leaves
      )
    ) %>%
    dplyr::select(-c(.data$uninfiltrated_leaves))

  df_melt <- df %>%
    dplyr::mutate(replicate = as.factor(1:dplyr::n())) %>%
    tidyr::pivot_longer(
      cols = -.data$replicate,
      names_to = "sample"
    )

  gg_fluor <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = df_tidied,
      ggplot2::aes(x = .data$sample, y = .data$fluorescence),
      fill = "#80b1d3",
      width = 0.5,
      alpha = 1
    ) +
    ggbeeswarm::geom_beeswarm(
      data = df_melt,
      ggplot2::aes(x = .data$sample, y = .data$value, color = .data$replicate)
    ) +
    ggplot2::labs(
      x = "Construct",
      y = "Flouresence",
      color = "Replicate"
    ) +
    ggplot2::scale_color_manual(
      values = c("#fdb462", "#bebada", "#fb8072")
    )
  # ggplot2::theme(
  #   axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust=1),
  #   axis.title.x= ggplot2::element_blank()
  # )

  return(gg_fluor)
}

#' Plot GFP
#'
#' @param df_with_pred_gfp Data frame with predicted GFP.
#' @param wildtype_sample Name of wildtype/control sample.
#'
#' @importFrom rlang .data
plot_bar_gfp <- function(df_with_pred_gfp, wildtype_sample) {
  # total_gfp <- sum(df_with_pred_gfp[["GFP (g/kg)"]])
  total_gfp <- df_with_pred_gfp %>%
    dplyr::filter(.data$Sample == wildtype_sample) %>%
    dplyr::pull(.data$`GFP (g/kg)`)

  gg_gfp <- df_with_pred_gfp %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$Sample, y = .data$`GFP (g/kg)`) +
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
