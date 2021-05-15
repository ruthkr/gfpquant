#' Title
#'
#' @param std_fluorescence
#' @param std_gfp
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
        .fns = mean
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
      gfp = (.data$fluorescence - intercept) / slope
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

  return(gg)
}

#' Title
#'
#' @param mat_sample_fluorescence Data
#'
#' @importFrom rlang .data
plot_bar_fluorescence <- function(mat_sample_fluorescence) {
  # browser()
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
    # ggplot2::scale_y_continuous(
    #   breaks = seq(0, 90500, 10000)
    # )

  return(gg_fluor)
}

#' Title
#'
#' @param df_with_pred_gfp Data
#'
#' @importFrom rlang .data
plot_bar_gfp <- function(df_with_pred_gfp) {
  gg_gfp <- df_with_pred_gfp %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = .data$Sample, y = .data$`GFP (g/kg)`) +
    ggplot2::geom_bar(
      stat = "identity",
      fill = "mediumseagreen",
      width = 0.5,
      alpha = 1
    ) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1250, 250))

  return(gg_gfp)
}
