#' @keywords internal
.save2disk <- function(
    assigned.xn,
    dirname,
    is.cell,
    is.subspot,
    xn.sce = NULL,
    xn.pos = NULL,
    cores = 1,
    verbose = TRUE
) {
  .root <- file.path(dirname, "cell")
  dir.create(.root, showWarnings = FALSE)
  
  if (is.cell) {
    assert_that(
      !is.null(xn.sce)
    )
    
    .save2disk_helper(
      xn = xn.sce,
      assignment = get_assignment2Spots(assigned.xn),
      assigment.count = get_countOfSpots(assigned.xn),
      dirname = .root,
      is.cell = TRUE,
      is.subspot = FALSE,
      cores = cores,
      verbose = verbose
    )
    
    if (is.subspot) {
      .save2disk_helper(
        xn = xn.sce,
        assignment = get_assignment2Subspots(assigned.xn),
        assigment.count = get_countOfSubspots(assigned.xn),
        dirname = .root,
        is.cell = TRUE,
        is.subspot = TRUE,
        cores = cores,
        verbose = verbose
      )
    }
  } else {
    assert_that(
      !is.null(xn.pos)
    )
    
    .save2disk_helper(
      xn = xn.pos,
      assignment = get_assignment2Spots(assigned.xn),
      assigment.count = get_countOfSpots(assigned.xn),
      dirname = .root,
      is.cell = FALSE,
      is.subspot = FALSE,
      cores = cores,
      verbose = verbose
    )
    
    if (is.subspot) {
      .save2disk_helper(
        xn = xn.pos,
        assignment = get_assignment2Subspots(assigned.xn),
        assigment.count = get_countOfSubspots(assigned.xn),
        dirname = .root,
        is.cell = FALSE,
        is.subspot = TRUE,
        cores = cores,
        verbose = verbose
      )
    }
  }
  
  invisible(NULL)
}


#' @keywords internal
.get_output_filename <- function(
    root,
    is.compressed = FALSE
) {
  .barcode.fname <- file.path(root, "barcodes.tsv")
  .feat.fname <- file.path(root, "features.tsv")
  .count.fname <- file.path(root, "matrix.mtx")
  .assignment.fname <- file.path(root, "assignment.tsv")
  .assignment.count.fname <- file.path(root, "assignment_count.tsv")
  
  ret <- c(
    .barcode.fname,
    .feat.fname,
    .count.fname,
    .assignment.fname,
    .assignment.count.fname
  )
  
  if (is.compressed) {
    ret <- paste(
      ret,
      "gz",
      sep = "."
    )
  }
  
  names(ret) <- c(
    "barcodes",
    "features",
    "mtx",
    "assignment",
    "assignment.count"
  )
  
  ret
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom Matrix writeMM
.write2disk <- function(
    barcodes,
    features,
    mtx,
    assignment,
    assignment.count
) {
  lapply(
    list(
      barcodes,
      features,
      mtx,
      assignment,
      assignment.count
    ),
    function(x) {
      assert_that(
        is.list(x),
        length(x) == 2,
        all(
          c(
            "obj",
            "fname"
          ) %in% names(x)
        )
      )
      
      invisible(NULL)
    }
  )
  
  write.table(
    x = barcodes$obj,
    file = barcodes$fname,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  system(paste(
    "gzip",
    barcodes$fname
  ))
  
  write.table(
    x = features$obj,
    file = features$fname,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  system(paste(
    "gzip",
    features$fname
  ))
  
  writeMM(
    mtx$obj,
    file = mtx$fname
  )
  system(paste(
    "gzip",
    mtx$fname
  ))
  
  write.table(
    x = assignment$obj,
    file = assignment$fname,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  system(paste(
    "gzip",
    assignment$fname
  ))
  
  write.table(
    x = assignment.count$obj,
    file = assignment.count$fname,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  system(paste(
    "gzip",
    assignment.count$fname
  ))
  
  invisible(NULL)
}


#' @keywords internal
#'
#' @importFrom assertthat assert_that
#' @importFrom pbapply pbsapply
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix Matrix rowSums
.save2disk_helper <- function(
    xn,
    assignment,
    assigment.count,
    dirname,
    is.cell,
    is.subspot,
    cores = 1,
    verbose = TRUE
) {
  assert_that(
    !is.null(xn),
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
  
  .output.name <- .get_output_filename(.root, is.compressed = FALSE)
  .output.name.compressed <- .get_output_filename(.root, is.compressed = TRUE)
  
  if (any(file.exists(.output.name.compressed))) {
    stop(paste(
      "Some of the following files already exist:",
      paste(
        .output.name.compressed,
        collapse = "\n"
      ),
      sep = "\n"
    ))
  }
  
  if (is.cell) {
    binned <- Matrix(pbsapply(
      setNames(nm = assigment.count[[barcode.label]]),
      function(x) {
        .assignment <- assignment %>%
          filter(
            !!as.name(barcode.label) == x
          )
        
        rowSums(counts(xn[, .assignment[["id"]]]))
      },
      cl = cores
    ))
    
    assert_that(
      identical(colnames(binned), assigment.count[[barcode.label]]),
      identical(rownames(binned), rownames(xn.sce))
    )
    
    .features <- data.frame(rowData(xn.sce))
  } else {
    .unique.spot.names <- assigment.count[[barcode.label]]
    .unique.gene.names <- unique(xn[["feature_name"]])
    
    binned <- Matrix(
      bin_transcripts(
        unique_spot_names = .unique.spot.names,
        unique_gene_names = .unique.gene.names,
        assignment_mole_names = assignment[["id"]],
        assignment_spot_names = assignment[[barcode.label]],
        feature_mole_names = xn[["transcript_id"]],
        feature_gene_names = xn[["feature_name"]],
        thread_num = cores,
        verbose = verbose
      ),
      dimnames = list(
        .unique.gene.names,
        .unique.spot.names
      )
    )
    
    .features <- data.frame(rownames(binned))
  }
  
  .write2disk(
    barcodes = list(
      obj = data.frame(colnames(binned)),
      fname = .output.name["barcodes"]
    ),
    features = list(
      obj = .features,
      fname = .output.name["features"]
    ),
    mtx = list(
      obj = binned,
      fname = .output.name["mtx"]
    ),
    assignment = list(
      obj = assignment[c("id", barcode.label)],
      fname = .output.name["assignment"]
    ),
    assignment.count = list(
      obj = assigment.count[c(barcode.label, "count")],
      fname = .output.name["assignment.count"]
    )
  )
  
  invisible(NULL)
}
