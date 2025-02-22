#' @import bslib
#' @noRd
theme <- function() {
  theme <- bs_theme(
    version = 4,
    bootswatch = "united",
    primary = "#606BD6",
    secondary = "#606BD6",
    # secondary = "#d0d4d8",
    # success = "#a4c689",
    # warning = "#fdbe4b",
    # info = "#cbd4dd",
    base_font = font_google("Roboto", wght = c(400, 500, 700)),
    heading_font = font_google("Roboto Condensed"),
    code_font = font_google("Fira Code"),
    "font-size-base" = "1rem"
  )

  return(theme)
}

#' @noRd
app_title <- function(title) {
  div(
    img(class = "navbar-custom-logo", src = "www/icon.svg"),
    div(
      class = "title-container",
      span(title),
      span(class = "title-version", as.character(utils::packageVersion("protAssay")))
    )
  )
}
