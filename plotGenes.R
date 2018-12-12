plotGenes <- function(seuratobject, genestoplot, path = ".", dark = F, width = 750, height = 750, reduction.use = "tsne", colours = c("grey", "red")) {
  genesnotpresent <- genestoplot[!(genestoplot %in% rownames(seuratobject@data))]
  if (length(genesnotpresent) != 0) {
    print(paste0("Not present in the dataset: ", genesnotpresent))
  }
  genestoplot <- genestoplot[genestoplot %in% rownames(seuratobject@data)]
  for (gene in genestoplot) {
    geneplotfilepath <- file.path(path, paste0(gene, "-", reduction.use, ".png"))
    if(!file.exists(geneplotfilepath)) {
      png(geneplotfilepath, width=width, height=height)
      FeaturePlot(seuratobject,
        features.plot = gene,
        nCol = 1,
        dark.theme = dark,
        cols.use = colours,
        reduction.use = reduction.use)
      dev.off()
    } else {
      print(paste0(geneplotfilepath, " already exists."))
    }
  }
}
