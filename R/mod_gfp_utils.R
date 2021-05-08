#' Title
#'
#' @param std_fluorescence
#' @param std_gfp
get_std_curve_value <- function(std_fluorescence, std_gfp) {
  df <- data.frame(
    fluorescence = std_fluorescence,
    gfp = std_gfp
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
#' @param mdata
#' @importFrom rlang .data
get_fluorescence_input <- function(mdata) {
  df_input <- as.data.frame(mdata, row.names = FALSE)

  # Get new df with the fluorescence value subtracted from infiltrated leave
  df_tidied <- df_input %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = mean
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = -.data$uninfiltrated_leave,
        .fns = ~ .x - .data$uninfiltrated_leave
      )
    ) %>%
    dplyr::select(-c("uninfiltrated_leave")) %>%
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "sample",
      values_to = "fluorescence"
    )

  return(df_tidied)
}

#' Function to get GFP value if fluorescence known
#'
#' @param df_tidied
#' @param fit Fit object.
#'
#' @importFrom rlang .data
predict_gfp_from_fluorescence <- function(df_tidied, fit) {
  intercept <- fit$coefficients[[1]]
  slope <- fit$coefficients[[2]]

  df_with_pred_gfp <- df_tidied %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      gfp = (.data$fluorescence - .data$intercept) / .data$slope
    )

  return(df_with_pred_gfp)
}



#' Title
#'
#' @param df_std_curve
#' @param df_with_pred_gfp
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
  # ggplot2::theme_bw()

  return(gg)
}
