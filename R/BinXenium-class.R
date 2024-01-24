#' @returns A BinXenium object.
#' 
#' @export
#' 
#' @importFrom assertthat assert_that
BinXenium <- function(
    vs.dir,
    xn.dir,
    bin.mole = c("cell", "trans", "both"),
    bin.level = c("spot", "subspot"),
    rm.feat.in.trans.pat = c("^NegControl.*", "^BLANK.*"),
    rot.mtx = NULL,
    rot.mtx.to = c("fullres", "hires", "lowres"),
    aligned.cells.file = NULL,
    aligned.trans.file = NULL,
    aligned.cells.coords.names = c(x = "x_centroid", y = "y_centroid"),
    aligned.trans.coords.names = c(x = "x_location", y = "y_location"),
    aligned.cells.to = c("fullres", "hires", "lowres"),
    aligned.trans.to = aligned.cells.to,
    dbg.load.trans = Inf
) {
  assert_that(
    !all(is.null(c(rot.mtx, aligned.cells.file, aligned.trans.file)))
  )
  
  bin.mole <- match.arg(bin.mole)
  bin.level <- match.arg(bin.level)
  rot.mtx.to <- match.arg(rot.mtx.to)
  aligned.cells.to <- match.arg(aligned.cells.to)
  aligned.trans.to <- match.arg(aligned.trans.to)
  
  vs.info <- .load_visium_info(vs.dir)
  
  xn.cell <- NULL
  if (bin.mole %in% c("cell", "both")) {
    xn.cell <- .load_xenium_molecular(
      xn.dir = xn.dir,
      rm.feat.in.trans.pat = NULL,
      rot.mtx = rot.mtx,
      rot.mtx.to = rot.mtx.to,
      aligned.file = aligned.cells.file,
      aligned.coords.names = aligned.cells.coords.names,
      aligned.to = aligned.cells.to,
      is.cell = TRUE,
      image.rad = get_imgRad(vs.info),
      scalef = list(
        hires = get_scalef(vs.info, "h"),
        lowres = get_scalef(vs.info, "l")
      )
    )
  }
  
  xn.trans <- NULL
  if (bin.mole %in% c("trans", "both")) {
    xn.trans <- .load_xenium_molecular(
      xn.dir = xn.dir,
      rm.feat.in.trans.pat = rm.feat.in.trans.pat,
      rot.mtx = rot.mtx,
      rot.mtx.to = rot.mtx.to,
      aligned.file = aligned.trans.file,
      aligned.coords.names = aligned.trans.coords.names,
      aligned.to = aligned.trans.to,
      is.cell = FALSE,
      image.rad = get_imgRad(vs.info),
      scalef = list(
        hires = get_scalef(vs.info, "h"),
        lowres = get_scalef(vs.info, "l")
      ),
      dbg.load.row.num = dbg.load.trans
    )
  }
  
  new(
    "BinXenium",
    vsInfo = vs.info,
    xnCell = xn.cell,
    xnTrans = xn.trans,
    is2Subspot = (bin.level == "subspot"),
    assignedXnCell = NULL,
    assignedXnTrans = NULL
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
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter
.load_visium_info <- function(vs.dir) {
  vs.dir <- file.path(vs.dir)
  
  if (grepl(".*outs$", vs.dir)) {
    vs.dir <- file.path(vs.dir, "spatial")
  }
  
  assert_that(grepl(".*spatial$", vs.dir))
  assert_that(file.exists(vs.dir))
  
  # load tissue position
  pos <- .load_visium_tissue_pos(vs.dir)
  
  # compute the directions of the coordinate system of array w.r.t. that of image
  array.direc <- .compute_array_direction(pos)
  
  # compute the angle in radians between the coordinate systems of array and image
  image.rad <- .compute_image_rad(pos)
  
  # load scale factors
  scalef.file <- file.path(vs.dir, "scalefactors_json.json")
  assert_that(file.exists(scalef.file))
  scalef <- fromJSON(file = scalef.file)
  
  # only keep spots that are in tissue
  pos <- pos %>%
    filter(
      in_tissue == 1
    )
  
  # scale coordinates
  pos <- .scale_coords(
    data = pos,
    coord.names = get_coords_names(
      is.xn = FALSE,
      prefix = NULL,
      use.names = "f"
    )[[1]],
    avail = "fullres",
    scalef = list(
      hires = scalef$tissue_hires_scalef,
      lowres = scalef$tissue_lowres_scalef
    )
  )
  
  # align the coordinate system of image to that of array
  pos <- .align_img2array(
    data = pos,
    radians = image.rad,
    prev.cols = get_coords_names(
      is.xn = FALSE,
      use.names = "f"
    ),
    transformed.cols = get_coords_names(
      is.xn = FALSE,
      prefix = "array_aligned",
      use.names = "f"
    )
  )
  
  # create subspots in the style of BayesSpace under the coordinate system of image
  subspot.pos <- .create_subspot_pos(
    radians = image.rad,
    radius = scalef$spot_diameter_fullres / 2,
    array.direc = array.direc,
    spot.coords = pos[c(
      "barcode",
      get_coords_names(
        is.xn = FALSE,
        use.names = "f"
      )[[1]]
    )]
  )
  
  new(
    "VisiumInfo",
    pos = pos,
    arrayDirec = array.direc,
    imgRad = image.rad,
    scalef = scalef,
    subspotPos = subspot.pos
  )
}

#' @keywords internal
#'
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate
.create_subspot_pos <- function(
    radians,
    radius,
    array.direc,
    spot.coords
) {
  .spot.coords.names <- get_coords_names(
    is.xn = FALSE,
    use.names = "f"
  )[[1]]
  
  assert_that(
    is.list(array.direc),
    all(c("col", "row") %in% names(array.direc)),
    all(c("barcode", .spot.coords.names) %in% colnames(spot.coords))
  )
  
  .spot.coords <- as.matrix(spot.coords[.spot.coords.names])
  colnames(.spot.coords) <- .spot.coords.names
  rownames(.spot.coords) <- spot.coords[["barcode"]]
  
  subspot.pos <- .create_subspots()
  
  .new.coords.names <- get_coords_names(
    is.xn = FALSE,
    prefix = "subspot",
    use.names = "f"
  )[[1]]
  .new.coords.names.col <- .new.coords.names[["x"]]
  .new.coords.names.row <- .new.coords.names[["y"]]
  
  # 1. make the coordinate system of array left-handed
  subspot.pos <- subspot.pos %>%
    mutate(
      "{.new.coords.names.col}" := subspot_rela_pxl_col * array.direc$col,
      "{.new.coords.names.row}" := subspot_rela_pxl_row * array.direc$row
    )
  
  # 2. rotate the adjusted coordinate system in the other direction
  .rotated.coords <- as.matrix(
    subspot.pos[.new.coords.names]
  ) %*% t(
    get_rot_mtx(radians, is.rh = FALSE, is.clock = FALSE)
  )
  colnames(.rotated.coords) <- .new.coords.names
  subspot.pos[.new.coords.names] <- .rotated.coords
  
  # 3. scale with the radius of the spot
  subspot.pos[[.new.coords.names.col]] <- subspot.pos[[.new.coords.names.col]] * radius
  subspot.pos[[.new.coords.names.row]] <- subspot.pos[[.new.coords.names.row]] * radius
  
  # 4. create coordinates of subspots for all spots
  .all.subspot.pos <- lapply(
    rownames(.spot.coords),
    function(.row) {
      .subspot.pos <- subspot.pos
      .subspot.pos[[.new.coords.names.col]] <- .subspot.pos[[.new.coords.names.col]] + .spot.coords[.row, .spot.coords.names[["x"]]]
      .subspot.pos[[.new.coords.names.row]] <- .subspot.pos[[.new.coords.names.row]] + .spot.coords[.row, .spot.coords.names[["y"]]]
      .subspot.pos[["subspot_barcode"]] <- paste(.row, .subspot.pos[["subspot_id"]], sep = ":")
      .subspot.pos[["subspot_barcode_bayesspace"]] <- paste(.row, .subspot.pos[["subspot_id_bayesspace"]], sep = ":")
      .subspot.pos[["barcode"]] <- .row
      
      .subspot.pos
    }
  )
  
  do.call(rbind, .all.subspot.pos)
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
  data
}


#' Compute the directions of the two axes of coordinate system of array with
#' respect to those of image. The coordinate system of image is assumed left-
#' handed.
#' 
#' @keywords internal
#'
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr group_by ungroup mutate
.compute_array_direction <- function(data) {
  assert_that(
    all(
      c("array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres") %in% colnames(data)
    )
  )
  
  .data <- data %>%
    group_by(array_row) %>%
    mutate(
      cor_col = cor(array_col, pxl_col_in_fullres)
    ) %>%
    ungroup() %>%
    group_by(array_col) %>%
    mutate(
      cor_row = cor(array_row, pxl_row_in_fullres)
    ) %>%
    ungroup()
  
  assert_that(
    min(.data$cor_row) * max(.data$cor_row) > 0,
    min(.data$cor_col) * max(.data$cor_col) > 0
  )
  
  list(
    row = ifelse(min(.data$cor_row) > 0, 1, -1),
    col = ifelse(min(.data$cor_col) > 0, 1, -1)
  )
}


#' Compute the degree in radians to align the coordinate systems of image to 
#' that of array. The coordinate system of image is assumed left-handed.
#' 
#' For a left-handed coordinate system, rotating by a positive angle means to
#' rotate clock-wisely, and vice versa. For a right-handed coordinate system,
#' rotating by a positive angle means to rotate counterclock-wisely, and vice
#' versa.
#' 
#' @keywords internal
#'
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter arrange
.compute_image_rad <- function(data) {
  assert_that(
    all(
      c("array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres") %in% colnames(data)
    )
  )
  
  radians <- vapply(
    unique(data$array_row),
    function(x) {
      .data <- data %>%
        filter(
          array_row == x
        ) %>%
        arrange(
          array_col
        )
      
      atan(
        (.data$pxl_row_in_fullres[dim(.data)[1]] - .data$pxl_row_in_fullres[1]) / (.data$pxl_col_in_fullres[dim(.data)[1]] - .data$pxl_col_in_fullres[1])
      )
    },
    FUN.VALUE = numeric(1)
  )
  
  assert_that(min(radians) * max(radians) > 0)
  
  mean(radians)
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
#' @importFrom dplyr left_join rename mutate
#' @importFrom SummarizedExperiment rowData
.load_xenium_molecular <- function(
    xn.dir,
    rm.feat.in.trans.pat,
    rot.mtx,
    rot.mtx.to,
    aligned.file,
    aligned.coords.names,
    aligned.to,
    is.cell,
    image.rad,
    scalef = list(
      hires = 1,
      lowres = 1
    ),
    dbg.load.row.num = Inf
) {
  if (is.cell) {
    sce <- .load_xenium2sce(xn.dir)
  }
  
  id <- ifelse(is.cell, "cell_id", "transcript_id")
  
  coords.names <- get_coords_names(
    is.xn = TRUE,
    is.cell = is.cell,
    use.names = "r"
  )[[1]]
  
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
      nrows = dbg.load.row.num,
      stringsAsFactors = FALSE,
      data.table = FALSE
    )
    
    if (!is.cell) {
      if (!is.null(rm.feat.in.trans.pat) && length(rm.feat.in.trans.pat) > 0) {
        raw.coords <- raw.coords[!grepl(
          paste(rm.feat.in.trans.pat, collapse = "|"),
          raw.coords$feature_name
        ), ]
      }
      
      raw.coords <- raw.coords %>%
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
    
    if (!is.cell) {
      if (!is.null(rm.feat.in.trans.pat) && length(rm.feat.in.trans.pat) > 0) {
        coords <- coords[!grepl(
          paste(rm.feat.in.trans.pat, sep = "|"),
          coords$feature_name
        ), ]
      }
      
      coords <- coords %>%
        mutate(
          transcript_id = as.character(transcript_id)
        )
    }
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
  
  # align the coordinate system of image to that of array
  coords <- .align_img2array(
    data = coords,
    radians = image.rad,
    prev.cols = get_coords_names(
      is.xn = TRUE,
      is.cell = is.cell,
      prefix = NULL,
      use.names = "f"
    ),
    transformed.cols = get_coords_names(
      is.xn = TRUE,
      is.cell = is.cell,
      prefix = "array_aligned",
      use.names = "f"
    )
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
        "XeniumMolecule",
        pos = coords,
        rotMtx = rot.mtx,
        alignedFile = aligned.file
      )
    )
  }
}
