#' @noRd
init_mat_std_curve <- function(is_prod = FALSE, assay = c("gfp", "bca")) {
  assay <- match.arg(assay)

  if (assay == "gfp") {
    levels <- c(100L, 200L, 400L, 600L, 800L)
    col_names <- c("GFP Level", "Fluorescence")
  } else if (assay == "bca") {
    levels <- c(25L, 125L, 250L, 500L, 750L)
    col_names <- c("BSA Level", "Signal")
  }

  if (is_prod) {
    # Production matrix
    if (assay == "gfp") {
      signal <- c(0L, 0L, 0L, 0L, 0L)
    } else if (assay == "bca") {
      signal <- c(0L, 0L, 0L, 0L, 0L)
    }
  } else {
    # Testing matrix
    if (assay == "gfp") {
      signal <- c(26442L, 40792L, 92383L, 135403L, 177405L)
    } else if (assay == "bca") {
      signal <- c(26442L, 40792L, 92383L, 135403L, 177405L)
    }
  }

  mat <- matrix(
    c(levels, signal),
    nrow = 5,
    ncol = 2,
    dimnames = list(NULL, col_names)
  )

  return(mat)
}

#' @noRd
init_mat_sample <- function(is_prod = FALSE) {
  if (is_prod) {
    # Production matrix
    mat <- matrix(
      rep("0", 3),
      nrow = 1,
      ncol = 3,
      dimnames = list(NULL, c("background_sample sample", paste("Sample", 1:2)))
    )
  } else {
    # Testing matrix

    mat <- matrix(
      c(3239, 3351, 3305, 94613, 93828, 93380, 26388, 26840, 27044, 33545, 34215, 34566),
      nrow = 3,
      ncol = 4,
      dimnames = list(NULL, c("background_sample sample", paste("Sample", 1:3)))
    )
  }
  return(mat)
}

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

  colnames(df_input)[1] <- "background_sample_column"

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
        .cols = -.data$background_sample_column,
        .fns = ~ .x - .data$background_sample_column
      )
    ) %>%
    dplyr::select(-c("background_sample_column")) %>%
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
