cpdbGetResults <- function(significant_means.txt, converted.pvalues, compared.celltypes, gene.csv, verbose = T) {
  significant_means_txt <- read.delim(significant_means.txt, header = F, sep = "\t", stringsAsFactors = F, na.strings = "")
  gene.ids <- read.csv(file = gene.csv, header = T, stringsAsFactors = F)

  comparisons <- c()
  for(r in 1:nrow(compared.celltypes)) {
    comparisons <- c(comparisons, paste0(compared.celltypes$C1[r], "_", compared.celltypes$C2[r]))
    comparisons <- c(comparisons, paste0(compared.celltypes$C2[r], "_", compared.celltypes$C1[r]))
  }

  selected.columns <- significant_means_txt[, significant_means_txt[1, ] %in% comparisons]
  keep.rows <- apply(selected.columns, 1, function(x) !(all(is.na(x))))
  meansfile <- cbind(significant_means_txt[keep.rows, 1:10], significant_means_txt[keep.rows, significant_means_txt[1, ] %in% comparisons])


  interaction.results <- c()
  cellpairs.results <- c()
  means.results <- c()

  for(comparison in 1:nrow(compared.celltypes)) {
    C1_C2 <- paste0(compared.celltypes$C1[comparison], "_", compared.celltypes$C2[comparison])
    C2_C1 <- paste0(compared.celltypes$C2[comparison], "_", compared.celltypes$C1[comparison])

    if(verbose) {
      print(C1_C2)
      print(C2_C1)
    }

    whichcolumn <- meansfile[1, ] == C1_C2 # results for this cell-pair
    whichcolumn[1] <- T # will be the interactions
    C1_C2_column <- meansfile[, whichcolumn]
    names(C1_C2_column) <- c("interactions", "means")
    C1_C2_column <- C1_C2_column[2:nrow(C1_C2_column), ] # remove header
    C1_C2_column <- C1_C2_column[!(is.na(C1_C2_column$means)), ] # don't need NA rows

    whichcolumn <- meansfile[1, ] == C2_C1 # repeat for reverse cell-pair
    whichcolumn[c(1, 3, 4)] <- T
    C2_C1_column <- meansfile[, whichcolumn]
    names(C2_C1_column) <- c("interactions", "partner_a", "partner_b", "means")
    C2_C1_column <- C2_C1_column[2:nrow(C2_C1_column), ]
    C2_C1_column <- C2_C1_column[!(is.na(C2_C1_column$means)), ]


    for(i in 1:length(C2_C1_column$interactions)) {
      # lookup r and l by uniprot name in gene.csv
      gene1 <- C2_C1_column$partner_a[i]
      gene1 <- sub("simple:", "", gene1)
      gene1 <- sub("complex:", "", gene1)
      gene2 <- C2_C1_column$partner_b[i]
      gene2 <- sub("simple:", "", gene2)
      gene2 <- sub("complex:", "", gene2)

      if (gene1 %in% gene.ids$uniprot) {
        gene1 <- gene.ids$hgnc_symbol[match(gene1, gene.ids$uniprot)]
        # else name remains
      }

      if (gene2 %in% gene.ids$uniprot) {
        gene2 <- gene.ids$hgnc_symbol[match(gene2, gene.ids$uniprot)]
        # else name remains
      }
      C2_C1_column$interactions[i] <- paste0(gene2, "_", gene1)
    }
    C2_C1_column <- C2_C1_column[ -c(2, 3) ]

    # One list for each cell-pair:
    C1_C2_combined <- rbind(C1_C2_column, C2_C1_column)

    v.size <- length(cellpairs.results)
    for(i in 1:length(C1_C2_combined$interactions)) {
      # Deconvolution step skipped. Dataframe will be compiled from these vectors:
      cellpairs.results[v.size + i] <- C1_C2
      interaction.results[v.size + i] <- C1_C2_combined$interactions[i]
      means.results[v.size + i] <- C1_C2_combined$means[i]

     }

  }

  interaction.return <- data.frame(Cellpairs = cellpairs.results, Interactions = interaction.results, Means = as.numeric(means.results))
  # Make a dictionary:
  converted.pvalues$key <- paste0(converted.pvalues[[1]], converted.pvalues[[2]])
  interaction.return$key <- paste0(interaction.return[[1]], interaction.return[[2]])
  converted.pvalues.filtered <- converted.pvalues[converted.pvalues$key %in% interaction.return$key, ]
  converted.pvalues.filtered$Interactions <- droplevels(converted.pvalues.filtered$Interactions)

  interaction.return <- merge(interaction.return, converted.pvalues.filtered)
  interaction.return <- interaction.return[ , !(names(interaction.return) %in% "key")]

  return(interaction.return)
}
