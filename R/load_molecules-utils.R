#' @name load_molecules-utils
#' 
#' @title Load molecules as a Molecule object
#' 
NULL


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr rename
.load_molecules <- function(
    dirname,
    image.rad,
    scalef = list(
      hires = 1,
      lowres = 1
    ),
    mole = c("xenium_cell", "xenium_transcript", "visium_hd"),
    rm.mole.feats.pat = c("^NegControl.*", "^BLANK.*"),
    rot.mtx = NULL,
    rot.mtx.to = c("fullres", "hires", "lowres"),
    aligned.file = NULL,
    aligned.coords.names = c(x = NULL, y = NULL),
    aligned.to = c("fullres", "hires", "lowres"),
    dbg.load.mole.nrows = Inf
) {
  # sanity check
  assert_that(
    file.exists(dirname),
    !is.null(scalef),
    all(c("hires", "lowres") %in% names(scalef)),
    any(c(is.null(rot.mtx), is.null(aligned.file))),
    is.numeric(dbg.load.mole.nrows)
  )
  
  mole <- match.arg(mole)
  
  coords.names <- get_coords_names(
    type = mole,
    use.names = "r"
  )[[1]]
  
  barcode.id <- .get_identifier_name(mole)
  
  # load spatial coordinates and clean column names
  if (!is.null(aligned.file)) {
    # if the molecules have been aligned manually to a template
    
    assert_that(
      file.exists(aligned.file),
      length(aligned.coords.names) == 2,
      all(c("x", "y") %in% names(aligned.coords.names))
    )
    
    aligned.to <- match.arg(aligned.to)
    
    new.coords.names <- paste(
      coords.names,
      aligned.to,
      sep = "_"
    )
    
    aligned.coords.names.x <- aligned.coords.names[["x"]]
    aligned.coords.names.y <- aligned.coords.names[["y"]]
    
    aligned.coords <- .load_spatial_coords(
      type = mole,
      filename = aligned.file,
      nrows = dbg.load.mole.nrows
    ) %>%
      rename(
        "{new.coords.names[1]}" := {{aligned.coords.names.x}},
        "{new.coords.names[2]}" := {{aligned.coords.names.y}}
      )
  } else {
    # load raw spatial coordinates
    
    coords <- .load_spatial_coords(
      type = mole,
      filename = .get_spatial_coords_filename(
        type = mole,
        dirname = dirname
      ),
      nrows = dbg.load.mole.nrows
    )
    
    if (!is.null(rot.mtx)) {
      assert_that(
        is.matrix(rot.mtx),
        nrow(rot.mtx) == ncol(rot.mtx),
        nrow(rot.mtx) == 2 || nrow(rot.mtx) == 3
      )
      
      rot.mtx.to <- match.arg(rot.mtx.to)
      
      new.coords.names <- paste(
        coords.names,
        rot.mtx.to,
        sep = "_"
      )
      
      # rotate the raw spatial coordinates
      aligned.coords <- .rotate_raw_coords(
        coords = coords,
        rot.mtx = rot.mtx,
        coords.names = coords.names,
        suffix = rot.mtx.to
      ) %>%
        mutate(
          id = coords[[barcode.id]]
        )
      
      coords <- coords %>%
        left_join(
          aligned.coords,
          by = setNames(nm = barcode.id, "id")
        )
    } else {
      new.coords.names <- paste(
        coords.names,
        "fullres",
        sep = "_"
      )
      
      coords.names.x <- coords.names[["x"]]
      coords.names.y <- coords.names[["y"]]
      
      coords <- coords %>%
        rename(
          "{new.coords.names[1]}" := {{coords.names.x}},
          "{new.coords.names[2]}" := {{coords.names.y}}
        )
    }
  }
  
  # scale coordinates
  coords <- .scale_coords(
    data = coords,
    coord.names = new.coords.names,
    avail = ifelse(
      !is.null(aligned.file),
      aligned.to,
      ifelse(
        !is.null(rot.mtx),
        rot.mtx.to,
        "fullres"
      )
    ),
    scalef = scalef
  )
  
  # align the coordinate system of image to that of array
  coords <- .align_img2array(
    data = coords,
    radians = image.rad,
    prev.cols = get_coords_names(
      type = mole,
      prefix = NULL,
      use.names = "f"
    ),
    transformed.cols = get_coords_names(
      type = mole,
      prefix = "array_aligned",
      use.names = "f"
    )
  )
  
  # create a SingleCellExperiment object
  sce <- .create_sce(
    dirname = dirname,
    type = mole,
    col.data = coords,
    rm.mole.feats.pat = rm.mole.feats.pat
  )
  assert_that(!is.null(sce))
  
  new(
    "Molecule",
    type = mole,
    sce = sce,
    rotMtx = rot.mtx,
    alignedFile = aligned.file
  )
}


#' @keywords internal
#' 
#' @importFrom arrow read_parquet
#' @importFrom data.table fread
#' @importFrom assertthat assert_that
.load_spatial_coords <- function(type, filename, nrows) {
  if (grepl("^xenium", type)) {
    ret <- fread(
      file = filename,
      nrows = nrows,
      stringsAsFactors = FALSE,
      data.table = FALSE
    )
  } else if (type == "visium_hd") {
    ret <- read_parquet(filename)
    
    if (nrows > 0 && nrows < nrow(ret)) {
      ret <- ret[1:nrows, ]
    }
  }
  
  assert_that(!is.null(ret))
  
  barcode.id <- .get_identifier_name(type)
  
  ret[[barcode.id]] <- as.character(ret[[barcode.id]])
  
  ret
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
.get_spatial_coords_filename <- function(type, dirname) {
  if (type == "xenium_cell") {
    fname <- file.path(dirname, "cells.csv.gz")
  } else if (type == "xenium_transcript") {
    fname <- file.path(dirname, "transcripts.csv.gz")
  } else if (mole == "visium_hd") {
    fname <- file.path(dirname, "tissue_positions.parquet")
  }
  
  assert_that(file.exists(fname))
  
  fname
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
  assert_that(
    all(coords.names %in% colnames(coords)),
    nrow(rot.mtx) == ncol(rot.mtx)
  )
  
  if (is.data.frame(coords)) {
    .coords <- as.matrix(coords[coords.names])
  } else if (is.matrix(coords)) {
    .coords <- coords[, coords.names]
  } else {
    stop(paste("Unsupported data structure:", class(coords)))
  }
  
  if (nrow(rot.mtx) == 3) {
    .coords <- cbind(.coords, 1)
  }
  
  .transformed.coords <- .coords %*% rot.mtx
  
  if (nrow(rot.mtx) == 3) {
    .transformed.coords <- .transformed.coords[, c(1, 2)]
  }
  
  colnames(.transformed.coords) <- paste(coords.names, suffix, sep = "_")
  
  as.data.frame(.transformed.coords)
}


#' Extract row and column indices of the count matrix from h5 file.
#'
#' @param idx Row index of corresponding element in the non-zero count matrix.
#' @param new.start Index of the start of each column corresponding to
#'     \code{idx} and the non-zero count matrix.
#' @param zero.based Whether the \code{} and \code{} are zero-based or not.
#'     (By default is TRUE)
#'
#' @return List of row (i) and column (j) indices of the non-zero elements
#'     in the count matrix.
#'
#' @keywords internal
#'
#' @importFrom tibble as_tibble
#' @importFrom tidyr uncount
.extract_indices <- function(idx, new.start, zero.based = TRUE) {
  if (length(idx) < 1) {
    return(NULL)
  }
  
  idx.cnts <- do.call(
    rbind,
    lapply(
      seq_len(length(new.start))[-1],
      function(x) c(x - ifelse(zero.based, 2, 1), new.start[[x]] - new.start[[x - 1]])
    )
  )
  colnames(idx.cnts) <- c("id", "n")
  
  return(
    list(
      i = idx,
      j = as.integer(uncount(as_tibble(idx.cnts), n)[[1]]),
      new.start = new.start
    )
  )
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom rhdf5 h5read
#' @importFrom magrittr `%>%`
#' @importFrom tibble column_to_rownames
#' @importFrom Matrix sparseMatrix
.read10Xh5 <- function(filename) {
  assert_that(file.exists(filename))
  
  .non.zero.indices <- .extract_indices(
    h5read(filename, "matrix/indices"),
    h5read(filename, "matrix/indptr")
  )
  
  .row.data <- data.frame(
    gene_id = h5read(filename, "matrix/features/id"),
    gene_name = h5read(filename, "matrix/features/name")
  ) %>%
    column_to_rownames("gene_id")
  
  .counts <- sparseMatrix(
    i = .non.zero.indices$i,
    j = .non.zero.indices$j,
    x = h5read(filename, "matrix/data"),
    dims = h5read(filename, "matrix/shape"),
    dimnames = list(
      rownames(.row.data),
      h5read(filename, "matrix/barcodes")
    ),
    index1 = FALSE
  )
  
  list(
    mtx = .counts,
    feat = .row.data
  )
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom SingleCellExperiment SingleCellExperiment
.create_sce <- function(dirname, type, col.data, rm.mole.feats.pat) {
  if (type == "xenium_cell") {
    sce <- .create_sce_h5(
      filename = file.path(
        dirname,
        "cell_feature_matrix.h5"
      ),
      col.data,
      rm.mole.feats.pat
    )
  } else if (type == "xenium_transcript") {
    sce <- .create_sce_xenium_transcripts(
      filename = file.path(
        dirname,
        "cell_feature_matrix",
        "features.tsv.gz"
      ),
      col.data,
      rm.mole.feats.pat
    )
  } else if (type == "visium_hd") {
    sce <- .create_sce_h5(
      filename = file.path(
        dirname,
        "filtered_feature_bc_matrix.h5"
      ),
      col.data,
      rm.mole.feats.pat
    )
  }
  
  assert_that(!is.null(sce))
  
  sce
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate filter select
.create_sce_h5 <- function(filename, col.data, rm.mole.feats.pat) {
  assert_that(
    file.exists(filename)
  )
  
  data <- .read10Xh5(filename)
  
  # remove some features
  if (!is.null(rm.mole.feats.pat) && length(rm.mole.feats.pat) > 0) {
    data$feat <- data$feat %>%
      mutate(
        keep = !grepl(
          paste(rm.mole.feats.pat, collapse = "|"),
          gene_name
        )
      ) %>%
      filter(
        keep
      ) %>%
      select(
        -keep
      )
    
    data$mtx <- data$mtx[rownames(data$feat), ]
  }
  
  SingleCellExperiment(
    assays = list(
      counts = data$mtx
    ),
    rowData = data$feat,
    colData = col.data
  )
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr rename select mutate
#' @importFrom tibble column_to_rownames
#' @importFrom Matrix sparseMatrix
.create_sce_xenium_transcripts <- function(
    filename,
    col.data,
    rm.mole.feats.pat
) {
  assert_that(
    file.exists(filename)
  )
  
  feat <- read.csv(
    file = filename,
    header = FALSE,
    sep = "\t"
  ) %>%
    rename(
      gene_id = V1,
      gene_name = V2
    ) %>%
    select(
      -V3
    ) %>%
    column_to_rownames("gene_id")
  
  # remove some features
  if (!is.null(rm.mole.feats.pat) && length(rm.mole.feats.pat) > 0) {
    feat <- feat %>%
      mutate(
        keep = !grepl(
          paste(rm.mole.feats.pat, collapse = "|"),
          gene_name
        )
      ) %>%
      filter(
        keep
      ) %>%
      select(
        -keep
      )
  }
  
  col.data <- col.data %>%
    mutate(
      gene_name_idx = match(col.data[["feature_name"]], feat$gene_name)
    ) %>%
    filter(
      !is.na(gene_name_idx)
    )
  
  
  mtx <- sparseMatrix(
    i = col.data[["gene_name_idx"]],
    j = seq_len(nrow(col.data)),
    x = 1,
    dims = c(
      length(feat$gene_name),
      nrow(col.data)
    ),
    dimnames = list(
      feat$gene_id,
      col.data[["transcript_id"]]
    ),
    index1 = TRUE
  )
  
  SingleCellExperiment(
    assays = list(
      counts = mtx
    ),
    rowData = feat,
    colData = col.data %>%
      select(
        -gene_name_idx
      )
  )
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom purrr compact
.parse_extra_args_mole <- function(args) {
  assert_that(is.list(args))
  
  .args <- compact(args[names(formals(.load_molecules))])
}
