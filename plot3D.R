#install.packages("rgl")
library(rgl)

maxPC <- 20 # number of PCs to use
seuratobject <- RunTSNE(seuratobject, dims.use = 1:maxPC, dim.embed = 3,
  do.fast = T, reduction.name = "tsne3D")
seuratCoordinates <- seuratobject@dr$tsne3D@cell.embeddings

# Get colours. See https://github.com/satijalab/seurat/issues/257
p <- Seurat::TSNEPlot(seuratobject, do.return = T)
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]]
cellColours <- pdata$colour

# 3D plot
plot3d(seuratCoordinates, col = cellColours, size = 0.1, box = T,
  xlab = "t-SNE x", ylab = "t-SNE y", zlab = "t-SNE z")
rgl.snapshot("tsne_3D.png", fmt = "png", top = T)

# Spin
play3d(spin3d(axis = c(0, 1, 0), rpm = 5), duration = 12)

# Save video
dir.create("video_tsne")
movie3d(spin3d(axis = c(0, 1, 0), rpm = 5), duration = 12, dir = "video_tsne",
  convert = NULL, movie = 'video')
