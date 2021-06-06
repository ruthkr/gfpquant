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

  colnames(df_input)[1] <- "subtrahend_column"

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
        .cols = -.data$subtrahend_column,
        .fns = ~ .x - .data$subtrahend_column
      )
    ) %>%
    dplyr::select(-c("subtrahend_column")) %>%
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
