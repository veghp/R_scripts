plotCircles <- function(circle.colours, do.border = F) {
  # Plot circles. 'circle.colours' is a vector of colours.
  len.col <- length(circle.colours)
  if(do.border) {
    border = circle.colours
   } else {border = NA }
  symbols(x = 2.5*1:len.col, y = rep(1, len.col), circles = rep(1, len.col),
  inches = F, bg = circle.colours, fg = border,
  xlab = NA, ylab = NA, xaxt = "n", yaxt = "n")
}
