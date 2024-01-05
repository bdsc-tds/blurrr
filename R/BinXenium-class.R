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
    aligned.trans.to = aligned.cells.to
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
    xn.cell <- .load_xenium_molecular(
      xn.dir = xn.dir,
      rot.mtx = rot.mtx,
      rot.mtx.to = rot.mtx.to,
      aligned.file = aligned.cells.file,
      aligned.coords.names = aligned.cells.coords.names,
      aligned.to = aligned.cells.to,
      is.cell = TRUE,
      scalef = list(
        hires = get_scalef(vs.info, "h"),
        lowres = get_scalef(vs.info, "l")
      )
    )
  }
  
  xn.trans <- NULL
  if (bin.mode %in% c("trans", "both")) {
    xn.trans <- .load_xenium_molecular(
      xn.dir = xn.dir,
      rot.mtx = rot.mtx,
      rot.mtx.to = rot.mtx.to,
      aligned.file = aligned.trans.file,
      aligned.coords.names = aligned.trans.coords.names,
      aligned.to = aligned.trans.to,
      is.cell = FALSE,
      scalef = list(
        hires = get_scalef(vs.info, "h"),
        lowres = get_scalef(vs.info, "l")
      )
    )
  }
  
  new(
    "BinXenium",
    vsInfo = vs.info,
    xnCell = xn.cell,
    xnTrans = xn.trans,
    is2Spot = ifelse(bin.level %in% c("spot", "both"), TRUE, FALSE),
    is2Subspot = ifelse(bin.level %in% c("subspot", "both"), TRUE, FALSE)
  )
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr rename
#' @importFrom purrr compact
.scale_coords <- function(
    data,
    coord.names,
    avail = c("fullres", "hires", "lowres"),
    scalef = list(
      hires = 1,
      lowres = 1
    )
) {
  assert_that(
    is.character(coord.names),
    length(coord.names) == 2,
    all(coord.names %in% colnames(data))
  )
  
  assert_that(
    is.list(scalef),
    length(scalef) == 2,
    "hires" %in% names(scalef),
    "lowres" %in% names(scalef)
  )

  avail <- match.arg(avail)
  
  # properly label the names of the coordinates
  new.coord.names <- coord.names
  is.unlabelled <- !grepl(avail, coord.names)
  if (any(is.unlabelled)) {
    new.coord.names[is.unlabelled] <- paste(
      coord.names[is.unlabelled],
      "avail",
      sep = "_"
    )
    
    coord.names.1 <- coord.names[1]
    coord.names.2 <- coord.names[2]
    
    data <- data %>%
      rename(
        "{new.coord.names[1]}" := {{coord.names.1}},
        "{new.coord.names[2]}" := {{coord.names.2}}
      )
  }
  
  all.coord.names <- sapply(
    setNames(nm = c(
      "fullres",
      "hires",
      "lowres"
    )),
    function(x) {
      if (x == avail) {
        return(new.coord.names)
      } else {
        return(gsub(
          avail,
          x,
          new.coord.names
        ))
      }
    },
    simplify = FALSE
  )
  
  if (avail != "fullres") {
    data[[all.coord.names$fullres[1]]] <- data[[new.coord.names[1]]] / scalef[[avail]]
    data[[all.coord.names$fullres[2]]] <- data[[new.coord.names[2]]] / scalef[[avail]]
  }
  
  do.call(
    cbind,
    c(
      data,
      compact(lapply(
        c("hires", "lowres"),
        function(x) {
          if (avail == "fullres" || avail != x) {
            ret <- data.frame(
              V1 = data[[all.coord.names$fullres[1]]] * scalef[[x]],
              V2 = data[[all.coord.names$fullres[2]]] * scalef[[x]]
            )
            colnames(ret) <- all.coord.names[[x]]
            
            return(ret)
          }
          
          return(NULL)
        }
      ))
    )
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
  
  # scale coordinates
  pos <- .scale_coords(
    data = pos,
    coord.names = c("pxl_row_in_fullres", "pxl_col_in_fullres"),
    avail = "fullres",
    scalef = list(
      hires = scalef$tissue_hires_scalef,
      lowres = scalef$tissue_lowres_scalef
    )
  )
  
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
  assert_that(length(unique(dim(rot.mtx))) == 1)
  
  if (is.data.frame(coords)) {
    .coords <- as.matrix(coords[coords.names])
  } else if (is.matrix(coords)) {
    .coords <- coords[, coords.names]
  } else {
    stop(paste("Unknown data structure:", class(coords)))
  }
  
  if (dim(rot.mtx)[1] == 3) {
    .coords <- cbind(.coords, 1)
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
  colnames(col_names) <- c("cell_id")
  
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
#' @importFrom dplyr left_join rename filter mutate
#' @importFrom SummarizedExperiment rowData
.load_xenium_molecular <- function(
    xn.dir,
    rot.mtx,
    rot.mtx.to,
    aligned.file,
    aligned.coords.names,
    aligned.to,
    is.cell,
    scalef = list(
      hires = 1,
      lowres = 1
    )
) {
  sce <- .load_xenium2sce(xn.dir)
  
  if (is.cell) {
    coords.names <- c("x_centroid", "y_centroid")
  } else {
    coords.names <- c("x_location", "y_location")
  }
  
  new.coords.names <- paste(
    coords.names,
    ifelse(
      !is.null(rot.mtx),
      rot.mtx.to,
      aligned.to
    ),
    sep = "_"
  )
  
  if (!is.null(rot.mtx)) {
    if (!is.null(aligned.file)) {
      warning("Both a rotation matrix and a file of aligned cells are provided. Using the rotation matrix.")
      
      aligned.file <- NULL
    }
    
    molecular_file <- file.path(xn.dir, ifelse(
      is.cell,
      "cells.csv.gz",
      "transcripts.csv.gz"
    ))
    
    assert_that(
      file.exists(molecular_file)
    )
    
    raw.coords <- fread(
      file = molecular_file,
      stringsAsFactors = FALSE,
      data.table = FALSE
    )
    
    id <- ifelse(is.cell, "cell_id", "transcript_id")
    
    if (!is.cell) {
      raw.coords <- raw.coords %>%
        filter(
          feature_name %in% rowData(sce)$gene_name
        ) %>%
        mutate(
          transcript_id = as.character(transcript_id)
        )
    }
    
    .transformed.coords <- .rotate_raw_coords(
      coords = raw.coords,
      rot.mtx = rot.mtx,
      coords.names = coords.names,
      suffix = rot.mtx.to
    ) %>%
      mutate(
        id = raw.coords[[id]]
      )
    
    coords <- raw.coords %>%
      left_join(
        .transformed.coords,
        by = setNames(nm = id, "id")
      )
  } else {
    aligned.file <- normalizePath(aligned.file)
    
    assert_that(file.exists(aligned.file))
    assert_that(!is.null(names(aligned.coords.names)))
    
    aligned.coords.names.x <- aligned.coords.names[["x"]]
    aligned.coords.names.y <- aligned.coords.names[["y"]]
    
    coords <- read.csv(aligned.file) %>%
      rename(
        "{new.coords.names[1]}" := {{aligned.coords.names.x}},
        "{new.coords.names[2]}" := {{aligned.coords.names.y}}
      )
  }
  
  # scale coordinates
  coords <- .scale_coords(
    data = coords,
    coord.names = new.coords.names,
    avail = ifelse(
      !is.null(rot.mtx),
      rot.mtx.to,
      aligned.to
    ),
    scalef = scalef
  )
  
  if (is.cell) {
    return(
      new(
        "XeniumCell",
        pos = coords,
        rotMtx = rot.mtx,
        alignedFile = aligned.file,
        sce = sce
      )
    )
  } else {
    return(
      new(
        "XeniumMolecular",
        pos = coords,
        rotMtx = rot.mtx,
        alignedFile = aligned.file
      )
    )
  }
}
