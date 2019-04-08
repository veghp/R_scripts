# Set parameters, then source file.
path.seuratobject <- "" # character string. Path to Seurat object file
annotation <- "" # Seuratobject metadata column name
metadata.filename <- "" # metadata export filename
count.filename <- "" # count export filename
gene.csv <- "" # path to gene.csv downloaded from https://www.cellphonedb.org/downloads/gene.csv
###############################################################################
library(Seurat)
library(plyr)
seuratobject <- readRDS(path.seuratobject)
seuratobject <- SetAllIdent(seuratobject, id = annotation)
if(!all(dim(seuratobject@data) == dim(seuratobject@raw.data))) { stop("Raw and processed data dimensions are not the same.") }

# Export metadata for CellPhoneDB:
##################################
metadata <- FetchData(seuratobject, vars.all = annotation)
metadata <- data.frame(Cell = rownames(metadata), cell_type = metadata[, annotation]) # string "cell_type" is required by CellPhoneDB
write.table(metadata, file = metadata.filename, row.names = FALSE, quote = F, sep = "\t")

# Export count data for CellPhoneDB:
####################################
# Load CellPhoneDB gene list:
cellphonedb.genes <- read.csv(gene.csv, stringsAsFactors = F)
# Get gene expression data:
counts <- as.matrix(seuratobject@raw.data)
genes <- unique(cellphonedb.genes$hgnc_symbol[cellphonedb.genes$hgnc_symbol %in% rownames(counts)])
# Replace data rownames (HGNC symbols) with ENS ids (CellPhoneDB uses ENS ids):
# One HGNC symbol can have multiple ENS ids, but it will be implicitly handled, i.e. first value will be used.
ens.rownames <- mapvalues(rownames(counts),
  from = cellphonedb.genes$hgnc_symbol,
  to = cellphonedb.genes$ensembl)

#grep("ENS", ens.rownames, value = T) # check

# Count normalization
column.sums <- colSums(counts)
cpm <- sweep(counts, 2, column.sums, "/")
cpm <- cpm * 10000
rm(counts)
# Export
data <- cbind("Gene" = ens.rownames, cpm)
write.table(data, count.filename, row.names = FALSE, quote = F, sep = "\t")
