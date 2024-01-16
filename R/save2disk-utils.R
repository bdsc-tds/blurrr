#' @keywords internal
.save2disk_cell <- function(
    xn.sce,
    assigned.xn.cell,
    dirname,
    is.subspot,
    cores = 1,
    verbose = TRUE
) {
  .root <- file.path(dirname, "cell")
  dir.create(.root, showWarnings = FALSE)
  
  .save2disk_cell_helper(
    xn.sce = xn.sce,
    assignment = get_assignment2Spots(assigned.xn.cell),
    assigment.count = get_countOfSpots(assigned.xn.cell),
    dirname = .root,
    is.subspot = FALSE,
    cores = cores,
    verbose = verbose
  )
  
  if (is.subspot) {
    .save2disk_cell_helper(
      xn.sce = xn.sce,
      assignment = get_assignment2Subspots(assigned.xn.cell),
      assigment.count = get_countOfSubspots(assigned.xn.cell),
      dirname = .root,
      is.subspot = TRUE,
      cores = cores,
      verbose = verbose
    )
  }
  
  invisible(NULL)
}

#' @keywords internal
#'
#' @importFrom assertthat assert_that
#' @importFrom pbapply pbsapply
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix Matrix rowSums writeMM
.save2disk_cell_helper <- function(
    xn.sce,
    assignment,
    assigment.count,
    dirname,
    is.subspot,
    cores = 1,
    verbose = TRUE
) {
  assert_that(
    !is.null(xn.sce),
    !is.null(assignment),
    !is.null(assigment.count),
    file.exists(dirname)
  )
  
  .root <- file.path(
    dirname,
    ifelse(
      is.subspot,
      "subspot",
      "spot"
    )
  )
  dir.create(.root, showWarnings = FALSE)
  
  barcode.label <- ifelse(
    is.subspot,
    "subspot_barcode_bayesspace",
    "barcode"
  )
  
  if (verbose) {
    message(paste(
      "Saving",
      ifelse(
        is.subspot,
        "subspot",
        "spot"
      ),
      "level aggregating resluts of Xenium cells to",
      .root
    ))
  }
  
  if (any(file.exists(
    paste(
      c(
        .barcode.fname,
        .feat.fname,
        .count.fname,
        .assignment.fname,
        .assignment.count.fname
      ),
      "gz",
      sep = "."
    )))) {
    stop(paste(
      "Some of the following files already exist:",
      paste(.barcode.fname, "gz", sep = "."),
      paste(.feat.fname, "gz", sep = "."),
      paste(.count.fname, "gz", sep = "."),
      paste(.assignment.fname, "gz", sep = "."),
      paste(.assignment.count.fname, "gz", sep = "."),
      sep = "\n"
    ))
  }
  
  binned <- Matrix(pbsapply(
    setNames(nm = assigment.count[[barcode.label]]),
    function(x) {
      .assignment <- assignment %>%
        filter(
          !!as.name(barcode.label) == x
        )
      
      rowSums(counts(xn.sce[, .assignment[["id"]]]))
    },
    cl = cores
  ))
  
  assert_that(
    identical(colnames(binned), assigment.count[[barcode.label]]),
    identical(rownames(binned), rownames(xn.sce))
  )
  
  .barcode.fname <- file.path(.root, "barcodes.tsv")
  .feat.fname <- file.path(.root, "features.tsv")
  .count.fname <- file.path(.root, "matrix.mtx")
  .assignment.fname <- file.path(.root, "assignment.tsv")
  .assignment.count.fname <- file.path(.root, "assignment_count.tsv")
  
  write.table(
    x = data.frame(colnames(binned)),
    file = .barcode.fname,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  system(paste("gzip", .barcode.fname))
  
  write.table(
    x = data.frame(rowData(xn.sce)),
    file = .feat.fname,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  system(paste("gzip", .feat.fname))
  
  writeMM(
    binned,
    file = .count.fname
  )
  system(paste("gzip", .count.fname))
  
  write.table(
    x = assignment[c("id", barcode.label)],
    file = .assignment.fname,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  system(paste("gzip", .assignment.fname))
  
  write.table(
    x = assigment.count[c(barcode.label, "count")],
    file = .assignment.count.fname,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  system(paste("gzip", .assignment.count.fname))
  
  invisible(NULL)
}
