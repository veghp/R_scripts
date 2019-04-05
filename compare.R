compare <- function(que, ref, filename = NULL) {

  results <- list()
  results[["intersect"]] <- intersect(que, ref)
  results[["que_only"]] <- setdiff(que, ref)
  results[["ref_only"]] <- setdiff(ref, que)

  if(!is.null(filename)) {
    if(!(is.character(filename) & length(filename) == 1)) {
      stop("Argument 'write' must be a nonempty string.")}
    sink(filename)
    print(results)
    sink()
  }

  return(results)
  
}
