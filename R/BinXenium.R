#' @returns A BinXenium object.
#' 
#' @export
#' 
#' @importFrom assertthat assert_that
BinXenium <- function(
    vs.dir,
    xn.dir,
    bin.mode = c("cell", "trans", "both"),
    bin.level = c("spot", "subspot", "both"),
    rot.mtx = NULL,
    rot.mtx.to = c("fullres", "hires", "lowres"),
    aligned.cells.file = NULL,
    aligned.trans.file = NULL,
    aligned.cells.coords.names = c(x = "x_centroid", y = "y_centroid"),
    aligned.trans.coords.names = c(x = "x_location", y = "y_location"),
    aligned.cells.to = c("fullres", "hires", "lowres"),
    aligned.trans.to = c("fullres", "hires", "lowres")
) {
  assert_that(
    !all(is.null(c(rot.mtx, aligned.cells.file, aligned.trans.file)))
  )
  
  bin.mode <- match.arg(bin.mode)
  bin.level <- match.arg(bin.level)
  rot.mtx.to <- match.arg(rot.mtx.to)
  aligned.cells.to <- match.arg(aligned.cells.to)
  aligned.trans.to <- match.arg(aligned.trans.to)
  
  vs.info <- .load_visium_info(vs.dir)
  
  xn.cell <- NULL
  if (bin.mode %in% c("cell", "both")) {
    xn.cell <- .load_xenium_cell(
      xn.dir,
      rot.mtx,
      rot.mtx.to,
      aligned.cells.file,
      aligned.cells.coords.names,
      aligned.cells.to
    )
  }
  
  xn.trans <- NULL
  if (bin.mode %in% c("trans", "both")) {
    xn.trans <- .load_xenium_transcript(
      xn.dir,
      rot.mtx,
      rot.mtx.to,
      aligned.trans.file,
      aligned.trans.coords.names,
      aligned.trans.to
    )
  }
  
  new(
    "BinXenium",
    vsInfo = vs.info,
    xnCell = xn.cell,
    xnTrans = xn.trans,
    toSpot = ifelse(bin.level %in% c("spot", "both"), TRUE, FALSE),
    toSubspot = ifelse(bin.level %in% c("subspot", "both"), TRUE, FALSE)
  )
}


#' @keywords internal
#' 
#' @importFrom rjson fromJSON
#' @importFrom assertthat assert_that
.load_visium_info <- function(vs.dir) {
  vs.dir <- file.path(vs.dir)
  
  if (grepl(".*outs$", vs.dir)) {
    vs.dir <- file.path(vs.dir, "spatial")
  }
  
  assert_that(grepl(".*spatial$", vs.dir))
  assert_that(file.exists(vs.dir))
  
  # load tissue position
  pos <- .load_visium_tissue_pos(vs.dir)
  
  # load scale factors
  scalef.file <- file.path(vs.dir, "scalefactors_json.json")
  assert_that(file.exists(scalef.file))
  scalef <- fromJSON(file = scalef.file)
  
  new(
    "VisiumInfo",
    pos = pos,
    scalef = scalef
  )
}


#' @keywords internal
#'
#' @importFrom utils read.csv
.load_visium_tissue_pos <- function(vs.dir) {
  if (file.exists(file.path(vs.dir, "tissue_positions_list.csv"))) {
    data <- read.csv(
      file.path(vs.dir, "tissue_positions_list.csv"),
      header = FALSE
    )
    colnames(data) <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
  } else if (file.exists(file.path(vs.dir, "tissue_positions.csv"))) {
    data <- read.csv(file.path(vs.dir, "tissue_positions.csv"))
  } else {
    stop("No file for spot positions found in ", vs.dir)
  }
  
  # Sanity check.
  if (
    abs(cor(data$array_row, data$pxl_row_in_fullres)) < abs(cor(data$array_row, data$pxl_col_in_fullres)) &&
    abs(cor(data$array_col, data$pxl_col_in_fullres)) < abs(cor(data$array_col, data$pxl_row_in_fullres))
  ) {
    message("Warning! The coordinates of arrays and pixels do not match. Swapping the pixel values between the row and column...")
    
    data <- transform(
      data,
      pxl_row_in_fullres = pxl_col_in_fullres,
      pxl_col_in_fullres = pxl_row_in_fullres
    )
  }
  
  rownames(data) <- data$barcode
  data[data$in_tissue > 0, ]
}


#' @keywords internal
#'
#' @importFrom assertthat assert_that
.rotate_raw_coords <- function(
    coords,
    rot.mtx,
    coords.names,
    suffix = "transformed"
) {
  assert_that(all(coords.names %in% colnames(coords)))
  assert_that(unique(dim(rot.mtx)) == 1)
  
  if (is.data.frame(coords)) {
    .coords <- as.matrix(coords[coords.names])
  } else if (is.matrix(coords)) {
    .coords <- coords[, coords.names]
  } else {
    stop(paste("Unknown data structure:", class(coords)))
  }
  
  if (dim(rot.mtx)[1] == 3) {
    .coords$offset <- 1
  }
  
  .transformed.coords <- .coords %*% rot.mtx
  
  if (dim(rot.mtx)[1] == 3) {
    .transformed.coords <- .transformed.coords[, c(1, 2)]
  }
  
  colnames(.transformed.coords) <- paste(coords.names, suffix, sep = "_")
  
  as.data.frame(.transformed.coords)
}


#' @keywords internal
#' 
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment SingleCellExperiment
.load_xenium2sce <- function(xn.dir) {
  assert_that(
    file.exists(file.path(xn.dir, "cell_feature_matrix")),
    file.exists(file.path(xn.dir, "cell_feature_matrix", "barcodes.tsv.gz")),
    file.exists(file.path(xn.dir, "cell_feature_matrix", "features.tsv.gz")),
    file.exists(file.path(xn.dir, "cell_feature_matrix", "matrix.mtx.gz"))
  )
  
  row_names <- read.table(
    file.path(xn.dir, "cell_feature_matrix", "features.tsv.gz"),
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  colnames(row_names) <- c("gene_id", "gene_name", "gene_info")
  
  col_names = read.table(
    file.path(xn.dir, "cell_feature_matrix", "barcodes.tsv.gz"),
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  counts <- readMM(file.path(xn.dir, "cell_feature_matrix", "matrix.mtx.gz"))
  dimnames(counts) <- list(row_names[[1]], col_names[[1]])
  
  SingleCellExperiment(
    assays = list(
      counts = counts
    ),
    rowData = row_names,
    colData = col_names
  )
}


#' @keywords internal
#' 
#' @importFrom utils read.csv
#' @importFrom assertthat assert_that
#' @importFrom data.table fread
#' @importFrom magrittr `%>%`
#' @importFrom dplyr left_join join_by rename
.load_xenium_cell <- function(
    xn.dir,
    rot.mtx,
    rot.mtx.to,
    aligned.cells.file,
    aligned.cells.coords.names,
    aligned.cells.to
) {
  sce <- .load_xenium2sce(xn.dir)
  
  if (!is.null(rot.mtx)) {
    if (!is.null(aligned.cells.file)) {
      warning("Both a rotation matrix and a file of aligned cells are provided. Using the rotation matrix.")
    }
    
    assert_that(
      file.exists(file.path(xn.dir, "cells.csv.gz"))
    )
    
    raw.coords <- fread(
      file = file.path(xn.dir, "cells.csv.gz"),
      stringsAsFactors = FALSE,
      data.table = FALSE
    )
    
    .transformed.coords <- .rotate_raw_coords(
      coords = raw.coords,
      rot.mtx = rot.mtx,
      coords.names = c("x_centroid", "y_centroid"),
      suffix = rot.mtx.to
    )
    
    .transformed.coords$cell_id <- raw.coords$cell_id
    
    coords <- raw.coords %>%
      left_join(
        .transformed.coords,
        by = join_by(cell_id)
      )
  } else {
    assert_that(!is.null(aligned.cells.file))
    
    new.coords.names <- paste(
      c("x_centroid", "y_centroid"),
      aligned.cells.to,
      sep = "_"
    )
    
    aligned.cells.coords.names.x <- aligned.cells.coords.names[["x"]]
    aligned.cells.coords.names.y <- aligned.cells.coords.names[["y"]]
    
    coords <- read.csv(aligned.cells.file) %>%
      rename(
        "{new.coords.names[1]}" := {{aligned.cells.coords.names.x}},
        "{new.coords.names[2]}" := {{aligned.cells.coords.names.y}}
      )
  }
  
  new(
    "XeniumCell",
    pos = coords,
    sce = sce
  )
}
