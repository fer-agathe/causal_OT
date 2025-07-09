wongBlack     <- "#000000"
wongGold      <- "#E69F00"
wongLightBlue <- "#56B4E9"
wongGreen     <- "#009E73"
wongYellow    <- "#F0E442"
wongBlue      <- "#0072B2"
wongOrange    <- "#D55E00"
wongPurple    <- "#CC79A7"


colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    if (grep(x = color, "^#")) {
      color <- deparse(substitute(color))
    }
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}