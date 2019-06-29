plotData <- function(data, x, y, z = NULL, line = "", spar = 1, xlim = c(0, max(data[, x]))) {
  if(is.null(z)) {
    z.dat <- rep(1, dim(data)[1])
  } else {
    z.dat <- data[, z]
  }

  symbols(x = data[, x], y = data[, y], circles = (sqrt(z.dat / pi)),
    inches = 1/10, ann = T, bg = "steelblue2", fg = "black",
    xlim = xlim, ylim = c(0, max(data[, y])),
    xlab = x, ylab = y)

  title(main = paste('Cor =', round(cor(data[, x], data[, y]), 2)))

  if(line == "line") {
    linm <- lm(data[, y] ~ data[, x])
    abline(a = linm$coefficients[1], b = linm$coefficients[2], col = "gray50")
  } else if(line == "spline") {
    spl <- smooth.spline(data[, x], data[, y], spar = spar)
    lines(spl, col = "red")
  }

  saveplot <- recordPlot()
  return(saveplot)
}
