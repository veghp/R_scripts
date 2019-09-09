# R_scripts

`cellphoneDB_plot` directory: plot CellPhoneDB results, see `cpdbPlot.html` for details.

`compare.R`: Calculate intersect and difference of two strings.

`exportSeurat.R`: Export Seurat object data and metadata as flat txt files.

`export_cellphonedb.R`: Export Seurat object for use with CellPhoneDB.

`markerDescription.R`: Takes a dataframe with a 'gene' column containing HGNC gene IDs (i.e. Seurat::FindAllMarkers() output) and returns it with added Entrez ID, Gene name (description), Cluster occurrences, Summary columns.

`markerDescription_ens97.R`: Version that works with the latest Ensembl database.

`pdt_heatmap.R`: Plot a heatmap from results of https://github.com/haniffalab/scRNA-seq_analysis/tree/master/pipelines/13_pseudotime.

`plot3D.R`: Interactive 3D dimensionality reduction (t-SNE) plot. 

`plotCircles.R`: Plot circles for use in figures.

`plotData.R`: Dotplot for use in https://github.com/veghp/Sequencing.

`plotGenes.R`: Save feature plots of specified genes.

`read10xsummary.R`: Reads 10x Cell Ranger metrics_summary.csv or web_summary.html files into a matrix (see https://github.com/veghp/Sequencing).
