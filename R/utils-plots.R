#' Plot standard curve and predictions
#'
#' @param df_std_curve Standard curve data frame.
#' @param df_with_pred_gfp Data frame with predicted GFP.
#' @param xlab X-axis name.
#' @param ylab Y-axis name.
#' @param fit Fit object.
#'
#' @importFrom rlang .data
plot_std_curve_and_pred <- function(df_std_curve, df_with_pred_gfp, fit, xlab = "GFP (ng)", ylab = "Flouresence") {
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
        y = -1e3, yend = .data$fluorescence
      ),
      linetype = "dashed"
    ) +
    # Horizontal lines
    ggplot2::geom_segment(
      data = df_with_pred_gfp,
      ggplot2::aes(
        x = -1e3, xend = .data$gfp,
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
      x = xlab,
      y = ylab
    ) +
    ggplot2::coord_cartesian(xlim = c(0, NA), ylim = c(0, NA))
  # ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.01), add = c(2, 0))) +
  # ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.01), add = c(2, 0)))

  return(gg)
}

#' Plot fluorescence
#'
#' @param mat_sample_fluorescence Input data.
#' @param ylab Y-axis name.
#'
#' @importFrom rlang .data
plot_bar_fluorescence <- function(mat_sample_fluorescence, ylab = "Flouresence") {
  df_tidied <- get_fluorescence_input(mat_sample_fluorescence)

  df <- as.data.frame(mat_sample_fluorescence) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = as.numeric
      )
    )

  colnames(df)[1] <- "background_sample_column"
  df <- df %>%
    dplyr::mutate(
      dplyr::across(
        .cols = -.data$background_sample_column,
        .fns = ~ .x - .data$background_sample_column
      )
    ) %>%
    dplyr::select(-c(.data$background_sample_column))

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
      y = ylab,
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
