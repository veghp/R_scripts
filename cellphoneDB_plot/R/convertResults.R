convertResults <- function(significant_means.txt, deconvoluted.txt, compared.celltypes, gene.csv, verbose = T) {
  significant_means_txt <- read.delim(significant_means.txt, header = F, sep = "\t", stringsAsFactors = F, na.strings = "")
  deconvoluted <- read.delim(deconvoluted.txt, header = T, stringsAsFactors = F)
  gene.ids <- read.csv(file = gene.csv, header = T, stringsAsFactors = F)

  comparisons <- c()
  for(r in 1:nrow(compared.celltypes)) {
    comparisons <- c(comparisons, paste0(compared.celltypes$C1[r], "_", compared.celltypes$C2[r]))
    comparisons <- c(comparisons, paste0(compared.celltypes$C2[r], "_", compared.celltypes$C1[r]))
  }

  selected.columns <- significant_means_txt[, significant_means_txt[1, ] %in% comparisons]
  keep.rows <- apply(selected.columns, 1, function(x) !(all(is.na(x))))
  meansfile <- cbind(significant_means_txt[keep.rows, 1:10], significant_means_txt[keep.rows, significant_means_txt[1, ] %in% comparisons])


  interaction.results <- c() # final dataframe will be compiled from these vectors
  cellpairs.results <- c()
  means.results <- c()
  rl.name.results <- c()

  for(comparison in 1:nrow(compared.celltypes)) {
    C1_C2 <- paste0(compared.celltypes$C1[comparison], "_", compared.celltypes$C2[comparison])
    C2_C1 <- paste0(compared.celltypes$C2[comparison], "_", compared.celltypes$C1[comparison])

    if(verbose) {
      print(C1_C2)
      print(C2_C1)
    }

    whichcolumn <- meansfile[1, ] == C1_C2 # results for this cell-pair
    whichcolumn[c(1, 3, 4)] <- T # interaction, partner_a and partner_b
    C1_C2_column <- meansfile[, whichcolumn]
    names(C1_C2_column) <- c("interactions", "partner_a", "partner_b", "means")
    C1_C2_column$partner_a_converted <- ""
    C1_C2_column$partner_b_converted <- ""
    C1_C2_column <- C1_C2_column[2:nrow(C1_C2_column), ] # remove header
    C1_C2_column <- C1_C2_column[!(is.na(C1_C2_column$means)), ] # don't need NA rows
    # Add gene ids / complex names:
    for(i in 1:length(C1_C2_column$interactions)) {
      # lookup r and l by uniprot name in gene.csv
      gene1 <- C1_C2_column$partner_a[i]
      gene1 <- sub("simple:", "", gene1)
      gene1 <- sub("complex:", "", gene1)
      gene2 <- C1_C2_column$partner_b[i]
      gene2 <- sub("simple:", "", gene2)
      gene2 <- sub("complex:", "", gene2)

      if (gene1 %in% gene.ids$uniprot) {
        C1_C2_column$partner_a_converted[i] <- gene.ids$hgnc_symbol[match(gene1, gene.ids$uniprot)]
      } else {
        C1_C2_column$partner_a_converted[i] <- gene1 # replace anyway, so string "complex:" is removed
      }

      if (gene2 %in% gene.ids$uniprot) {
        C1_C2_column$partner_b_converted[i] <- gene.ids$hgnc_symbol[match(gene2, gene.ids$uniprot)]
      } else {
        C1_C2_column$partner_b_converted[i] <- gene2 # replace anyway, so string "complex:" is removed
      }
    }

    whichcolumn <- meansfile[1, ] == C2_C1 # repeat for reverse cell-pair
    whichcolumn[c(1, 3, 4)] <- T
    C2_C1_column <- meansfile[, whichcolumn]
    names(C2_C1_column) <- c("interactions", "partner_b", "partner_a", "means") # swap gene names by changing a and b order.
    C2_C1_column$partner_a_converted <- ""
    C2_C1_column$partner_b_converted <- ""
    C2_C1_column <- C2_C1_column[2:nrow(C2_C1_column), ]
    C2_C1_column <- C2_C1_column[!(is.na(C2_C1_column$means)), ]
    # Add gene ids / complex names:
    for(i in 1:length(C2_C1_column$interactions)) {
      # lookup r and l by uniprot name in gene.csv
      gene1 <- C2_C1_column$partner_a[i]
      gene1 <- sub("simple:", "", gene1)
      gene1 <- sub("complex:", "", gene1)
      gene2 <- C2_C1_column$partner_b[i]
      gene2 <- sub("simple:", "", gene2)
      gene2 <- sub("complex:", "", gene2)

      if (gene1 %in% gene.ids$uniprot) {
        C2_C1_column$partner_a_converted[i] <- gene.ids$hgnc_symbol[match(gene1, gene.ids$uniprot)]
        # else name stays that complex
      } else {
        C2_C1_column$partner_a_converted[i] <- gene1 # replace anyway, so string "complex:" is removed
      }

      if (gene2 %in% gene.ids$uniprot) {
        C2_C1_column$partner_b_converted[i] <- gene.ids$hgnc_symbol[match(gene2, gene.ids$uniprot)]
        # else name stays that complex
      } else {
        C2_C1_column$partner_b_converted[i] <- gene2 # replace anyway, so string "complex:" is removed
      }
      C2_C1_column$interactions[i] <- paste0(C2_C1_column$partner_a_converted[i], "::", C2_C1_column$partner_b_converted[i]) # replace with swapped gene names (genes swapped above)
    }


    # One list for each cell-pair:
    C1_C2_combined <- rbind(C1_C2_column, C2_C1_column)

    # Deconvolution:
    for(rl.interaction in 1:length(C1_C2_combined$interactions)) {
      rl.meanvalue <- C1_C2_combined$means[rl.interaction] # to be assigned to the interaction
      rl.name <- C1_C2_combined$interactions[rl.interaction]
      gene_a <- C1_C2_combined$partner_a_converted[rl.interaction]
      gene_b <- C1_C2_combined$partner_b_converted[rl.interaction]

      is.rl1.complex <- gene_a %in% deconvoluted$complex_name
      is.rl2.complex <- gene_b %in% deconvoluted$complex_name

      if(is.rl1.complex) {
        rl1.genes <- unique(deconvoluted[deconvoluted$complex_name == gene_a, ]$gene_name)
      } else {
        rl1.genes <- gene_a
      }

      if(is.rl2.complex) {
        rl2.genes <- unique(deconvoluted[deconvoluted$complex_name == gene_b, ]$gene_name)
      } else {
        rl2.genes <- gene_b
      }

      added.lines <- length(rl1.genes) * length(rl2.genes)

      if(!(is.rl1.complex) & !(is.rl2.complex)) { # none complex
        interaction.results <- c(interaction.results, paste0(rl1.genes, "::", rl2.genes))
      }

      if(is.rl1.complex & !(is.rl2.complex)) { # 1 complex
        interaction.results <- c(interaction.results, paste0(rl1.genes, "::", rl2.genes))
      }

      if(!(is.rl1.complex) & is.rl2.complex) { # 2 complex
        interaction.results <- c(interaction.results, paste0(rl1.genes, "::", rl2.genes))
      }

      if(is.rl1.complex & is.rl2.complex) { # both complex
        for(i in 1:length(rl1.genes)) {
          interaction.results <- c(interaction.results, paste0(rl1.genes[i], "::", rl2.genes))
        }
      }

      rl.name.results <- c(rl.name.results, rep(rl.name, times = added.lines))
      means.results <- c(means.results, rep(rl.meanvalue, times = added.lines))
      cellpairs.results <- c(cellpairs.results, rep(C1_C2, times = added.lines))

     }

  }

  interaction.return <- data.frame(Cellpairs = cellpairs.results, Interactions = interaction.results, Means = as.numeric(means.results), RLnames = rl.name.results)
  return(interaction.return)
}
