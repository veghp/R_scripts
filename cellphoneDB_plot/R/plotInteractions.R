plotInteractions <- function(cpdb.results) {
  require("ggplot2")
  cols.use <- c("red", "blue")
  p <- ggplot(data = cpdb.results, mapping = aes(x = Cellpairs,  y = Interactions)) +
    geom_point(mapping = aes(size = Means, color = Pvalues), shape = 15) +
    scale_color_gradient(low = cols.use[1], high = cols.use[2], limits=c(0, 0.05)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}
