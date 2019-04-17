plotInteractions2 <- function(cpdb.results) {
  require("ggplot2")
  cols.use <- c("red", "blue")
  p <- ggplot(data = cpdb.results.subset, mapping = aes(x = Cellpairs, y = Interactions)) +
    geom_point(mapping = aes(size = log2(Means), color = Pvalues)) +
    scale_color_gradient(low = cols.use[1], high = cols.use[2], limits=c(0, 0.05)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank())
  return(p)
}
