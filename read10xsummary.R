read10xsummary <- function(tenx.samples, path = NULL, what = "csv") {
  # This function reads 10x cellranger metrics_summary.csv or web_summary.html
  # files into a matrix.
  # tenx.samples: vector of 10x sample names (directories)
  # path: vector of csv or html file paths. Specify if current working dir does
  # not contain sample dirs (tenx.samples), or if the files are somewhere else.
  # what: read "csv" or "html" files
  # Functions:
  removeComma <- function(txt) {gsub(",", "", txt, fixed = TRUE)}

  removePercent <- function(txt) {
    if(all(grepl("%", txt))) {
      txt <- gsub("%", "", txt, fixed = TRUE)
      txt <- as.numeric(txt) / 100
      return(txt)
    } else {return(txt)
    }
  }


  if(what == "csv") {
    # Construct filepaths
    if(is.null(path)) {
      path <- tenx.samples
      path <- file.path(tenx.samples, "metrics_summary.csv")
    }
    summary.csv.files <- path

    summary.csv <- do.call(rbind,
      lapply(summary.csv.files,
        function(x) read.csv(x, stringsAsFactors = FALSE)))
    summary.csv <- apply(summary.csv, 2, removeComma)
    summary.csv <- apply(summary.csv, 2, removePercent)
    summary.csv <- apply(summary.csv, 2, as.numeric)
    rownames(summary.csv) <- tenx.samples

    return(summary.csv)


  } else if(what == "html") {
    # Construct filepaths
    if(is.null(path)) {
      path <- tenx.samples
      path <- file.path(tenx.samples, "web_summary.html")
    }
    web_summary.html.files <- path

    to.match <- c(
      "<td>Number of Reads</td>",
      "<td>Sequencing Saturation</td>",
      "<td>Estimated Number of Cells</td>",
      "<td>Fraction Reads in Cells</td>",
      "<td>Mean Reads per Cell</td>",
      "<td>Median Genes per Cell</td>",
      "<td>Total Genes Detected</td>",
      "<td>Median UMI Counts per Cell</td>") # will get these from html
    to.match <- paste(to.match, collapse="|")

    summary.html <- data.frame(row.names = c(
      "Number.of.Reads",
      "Sequencing.Saturation",
      "Estimated.Number.of.Cells",
      "Fraction.Reads.in.Cells",
      "Mean.Reads.per.Cell",
      "Median.Genes.per.Cell",
      "Total.Genes.Detected",
      "Median.UMI.Counts.per.Cell"))

    for(f in 1:length(web_summary.html.files)) {

      command.txt <- paste0(
        "grep -A 1 -E '", to.match, "' ",
        web_summary.html.files[f])

      grep.txt <- system(command.txt, wait = T, intern = T)
      grep.txt <- gsub("      <td>", "", grep.txt)
      grep.txt <- gsub("</td>", "", grep.txt)
      grep.txt <- gsub(",", "", grep.txt) # remove junk txt

      to.df <- grep.txt[seq(2, 23, 3)] # every 4th element is a number we need
      summary.html[, tenx.samples[f]] <- to.df
    }

    summary.html <- apply(summary.html, 1, removePercent)
    summary.html <- apply(summary.html, 2, as.numeric)
    rownames(summary.html) <- tenx.samples

    return(summary.html)


  } else {stop("Argument 'what' must be one of 'csv', 'html'")}
}
