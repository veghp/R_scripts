rlTableFromVectors <- function(v.1, v.2) {
  nrows.rl.matrix <- length(v.1) * length(v.2) # N*M pairs
  dummy <- rep("", nrows.rl.matrix)
  compared.celltypes <- data.frame(C1=dummy, C2=dummy, stringsAsFactors=F) # lists comparisons
  counter <- 1
  for(i in 1:length(v.1)) {
    celltype1 <- v.1[i]
    for(j in 1:length(v.2)) {
      celltype2 <- v.2[j]
      compared.celltypes$C1[counter] <- celltype1
      compared.celltypes$C2[counter] <- celltype2
      counter <- counter + 1
    }
  }
  rownames(compared.celltypes) <- paste0(compared.celltypes$C1, "::", compared.celltypes$C2)
  return(compared.celltypes)
}
