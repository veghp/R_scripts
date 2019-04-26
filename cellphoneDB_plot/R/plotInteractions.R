plotInteractions3 <- function(cpdb.results, type = 1) {
  require("ggplot2")
  cols.use <- c("red", "blue")
  p <- ggplot(data = cpdb.results, mapping = aes(x = Cellpairs,  y = Interactions))
  if(type == 1) {
    p <- p + geom_point(mapping = aes(size = Means, color = Pvalues), shape = 15) +
    scale_color_gradient(low = cols.use[1], high = cols.use[2], limits = c(0, 0.05)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return(p)
  } else if(type == 2) {
    p <- p + geom_point(mapping = aes(size = log2(Means), color = Pvalues)) +
      scale_color_gradient(low = cols.use[1], high = cols.use[2], limits = c(0, 0.05)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank())
    return(p)

  } else { stop("Argument 'type' must be 1 or 2.") }
}
