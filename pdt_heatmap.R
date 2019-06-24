# Plot heatmap with select genes of results from: https://github.com/haniffalab/Single-cell-RNAseq-data-analysis-bundle/tree/master/pipelines/13_pseudotime
library("ggplot2")
seuratobject <- readRDS(file = "path/to/rds")
path <- "path/to/pdt/output"
plottingmat <- readRDS(file.path(path, "ploting_material.RDS")) # [sic]

# Optionally find TFs in your genelist:
# TF were downloaded from: http://humantfs.ccbr.utoronto.ca/download.php
TF.genes <- scan("TF_names_v_1.01_ccbr_utoronto.txt", what = character())
TF <- plottingmat$beautiful_result_norm$GeneNames[
  plottingmat$beautiful_result_norm$GeneNames %in% TF.genes]

# selected.genes.txt is a text file listing one gene per line.
{
selected.gene.list <- scan("selected.genes.txt", what = character(),
  sep = "\n", blank.lines.skip = T, comment.char = "#", strip.white = T)
selected.gene.list <- unique(selected.gene.list)

subsetplotmat <- plottingmat$beautiful_result_norm[
  plottingmat$beautiful_result_norm$GeneNames %in% selected.gene.list, ]
subsetplotmat$GeneNames <- droplevels(subsetplotmat$GeneNames)
subsetplotmat$GeneNames <- factor(subsetplotmat$GeneNames,
  levels = rev(selected.gene.list))

# The following section is adapted from: https://github.com/haniffalab/Single-cell-RNAseq-data-analysis-bundle/blob/master/pipelines/13_pseudotime/pseudotime.R#L270 commit b86d20dc87d35820daac178a93e46badf99216ab
plot.genes = ggplot(data = subsetplotmat, aes(x = Pseudotime, y = GeneNames))
plot.genes = plot.genes + geom_tile(aes(fill = ExpressionValue),
  width = 1.001, height = 1.001)
plot.genes = plot.genes + scale_fill_gradient2(low = "deepskyblue",
  high = "firebrick3",
  mid = "darkolivegreen3",
  midpoint = 0.5,
  name = "Minmax normalized gene expression")
plot.genes = plot.genes + theme(legend.position = "bottom",
  # legend.text = element_text(size = 25, angle = 90),
  # legend.title = element_text(size = 25),
  # legend.key.width = unit(2, "cm"),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 0),
  axis.text.y = element_text(size = 8))

plot.genes

}

height = 5; width = 2.7
pdf("dpt_heatmap.pdf", height = height, width = width)
plot.genes
dev.off()

svg("dpt_heatmap.svg", height = height, width = width)
plot.genes
dev.off()

postscript("dpt_heatmap.ps", height = height, width = width)
plot.genes
dev.off()

png("dpt_heatmap.png", height = 500, width = 300)
plot.genes
dev.off()
