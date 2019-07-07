plotData <- function(d, x, y, z = NULL, line = F, spline = F, spar = 1, xlim = c(0, max(d[, x]))) {
  if(is.null(z)) {
    z.dat <- rep(1, dim(d)[1])
  } else {
    z.dat <- d[, z]
  }

  symbols(x = d[, x], y = d[, y], circles = (sqrt(z.dat / pi)),
    inches = 1/10, ann = T, bg = "steelblue2", fg = "black",
    xlim = xlim, ylim = c(0, max(d[, y])),
    xlab = x, ylab = y, las = 1)

  title(main = paste('Cor =', round(cor(d[, x], d[, y]), 2)))

  if(line) {
    linm <- lm(d[, y] ~ d[, x])
    abline(a = linm$coefficients[1], b = linm$coefficients[2], col = "gray50")
  }
  if(spline) {
    spl <- smooth.spline(d[, x], d[, y], spar = spar)
    lines(spl, col = "red")
  }

  saveplot <- recordPlot()
  return(saveplot)
}
