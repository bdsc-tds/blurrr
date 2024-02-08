#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom rjson fromJSON
.load_blur_template <- function(
  dirname,
  blur = c("visium", "visium_hd"),
  blur2subspot = FALSE
) {
  blur <- match.arg(blur)
  
  if (blur == "visium") {
    assert_that(grepl(".*outs$", dirname))
  } else {
    assert_that(grepl(".*square_\\d+um$", dirname))
  }
  
  dirname <- file.path(dirname, "spatial")
  assert_that(file.exists(dirname))
  
  # load tissue position
  pos <- .load_tissue_pos(
    dirname = dirname,
    blur = blur
  )
  
  # compute the directions of the coordinate system of array w.r.t. that of image
  array.direc <- .compute_array_direction(pos)
  
  # compute the angle in radians between the coordinate systems of array and image
  image.rad <- .compute_image_rad(pos)
  
  # load scale factors
  scalef.file <- file.path(dirname, "scalefactors_json.json")
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
      type = blur,
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
      type = blur,
      use.names = "f"
    ),
    transformed.cols = get_coords_names(
      type = blur,
      prefix = "array_aligned",
      use.names = "f"
    )
  )
  
  # create subspots in the style of BayesSpace under the coordinate system of image
  subspot.pos <- NULL
  if (blur2subspot) {
    subspot.pos <- .create_subspot_pos(
      radians = image.rad,
      diameter = scalef$spot_diameter_fullres,
      array.direc = array.direc,
      spot.coords = pos[c(
        "barcode",
        get_coords_names(
          type = blur,
          use.names = "f"
        )[[1]]
      )],
      blur = blur
    )
  }
  
  # constructor
  new(
    "Template",
    type = blur,
    pos = pos,
    arrayDirec = array.direc,
    imgRad = image.rad,
    scalef = scalef,
    subspotPos = subspot.pos
  )
}


#' @keywords internal
#' 
#' TODO
.create_blur_template <- function() {}


#' @keywords internal
.load_tissue_pos <- function(
    dirname,
    blur = c("visium", "visium_hd")
) {
  blur <- match.arg(blur)
  
  if (blur == "visium") {
    data <- .load_tissue_pos_visium(dirname)
  } else {
    data <- .load_tissue_pos_visium_hd(dirname)
  }
  
  # sanity check
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


#' @keywords internal
#'
#' @importFrom utils read.csv
.load_tissue_pos_visium <- function(dirname) {
  if (file.exists(file.path(dirname, "tissue_positions_list.csv"))) {
    data <- read.csv(
      file.path(dirname, "tissue_positions_list.csv"),
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
  } else if (file.exists(file.path(dirname, "tissue_positions.csv"))) {
    data <- read.csv(file.path(dirname, "tissue_positions.csv"))
  } else {
    stop("No file for spot positions found in ", dirname)
  }
  
  data
}


#' @keywords internal
#'
#' @importFrom assertthat assert_that
#' @importFrom arrow read_parquet
.load_tissue_pos_visium_hd <- function(dirname) {
  filename <- file.path(dirname, "tissue_positions.parquet")
  
  assert_that(
    file.exists(filename)
  )
  
  read_parquet(filename)
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


#' Compute the degree in radians to align the coordinate system of image to 
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
#' TODO
.create_subspot_pos <- function(
    radians,
    diameter,
    array.direc,
    spot.coords,
    blur = c("visium", "visium_hd")
) {
  blur <- match.arg(blur)
  
  if (blur == "visium") {
    data <- .create_subspot_pos_visium(
      radians = radians,
      radius = diameter / 2,
      array.direc = array.direc,
      spot.coords = spot.coords
    )
  } else {
    data <- .create_subspot_pos_visium_hd()
  }
  
  data
}


#' @keywords internal
#'
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate
.create_subspot_pos_visium <- function(
    radians,
    radius,
    array.direc,
    spot.coords
) {
  .spot.coords.names <- get_coords_names(
    type = "visium",
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
  
  subspot.pos <- .create_subspots_visium()
  
  .new.coords.names <- get_coords_names(
    type = "visium",
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
#' TODO
.create_subspot_pos_visium_hd <- function() {}
