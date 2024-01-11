#' Create coordinates of points of a hexagon with respect to the coordinate
#' system of array, assuming unit radius of the spot.
#'
#' Always starting from the positive direction of the y axis to the positive
#' direction of the x axis, i.e., clockwise w.r.t. a right-handed coordinate
#' system, or counterclockwise w.r.t. a left-handed coordinate system.
#' 
#' @keywords internal
#' 
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate lead
.create_subspots <- function() {
  data.frame(
    subspot_id = seq_len(6),
    subspot_id_bayesspace = c(1, 5, 3, 4, 6, 2),
    subspot_rela_pxl_col = c(
      0,
      sqrt(3) / 2,
      sqrt(3) / 2,
      0,
      -sqrt(3) / 2,
      -sqrt(3) / 2
    ),
    subspot_rela_pxl_row = c(
      1,
      0.5,
      -0.5,
      -1,
      -0.5,
      0.5
    )
  ) %>%
    mutate(
      third_pt_id = lead(subspot_id, default = min(subspot_id))
    )
}


#' @export
#'
#' @importFrom assertthat assert_that
get_coords_names <- function(
    is.xn,
    is.cell = TRUE,
    is.subspot = FALSE,
    prefix = NULL,
    use.names = c("r", "f", "h", "l")
) {
  assert_that(is.null(prefix) || is.character(prefix))
  
  .names <- c("raw", "fullres", "hires", "lowres")
  use.names <- match.arg(use.names, .names, several.ok = TRUE)
  
  if (is.xn) {
    if (is.cell) {
      coords.names <- c(x = "x_centroid", y = "y_centroid")
    } else {
      coords.names <- c(x = "x_location", y = "y_location")
    }
  } else {
    coords.names <- c(x = "pxl_col_in", y = "pxl_row_in")
    
    if (is.subspot) {
      prefix = "subspot"
      use.names = "fullres"
    }
  }
  
  sapply(
    setNames(nm = use.names),
    function(x) {
      if (x == "raw") return(coords.names)
      
      ret <- paste(coords.names, x, sep = "_")
      if (!is.null(prefix) && length(prefix) > 0) {
        ret <- paste(prefix, ret, sep = "_")
      }
      names(ret) <- c("x", "y")
      
      return(ret)
    },
    simplify = FALSE
  )
}


#' Get a rotation matrix for rotating a coordinate system by a degree.
#' 
#' Usage:
#' 1. new_coordinates (2xN) = rotation_matrix (2x2) %*% old_coordinates (2xN)
#' 2. new_coordinates (Nx2) = old_coordinates (Nx2) %*% t(rotation_matrix) (2x2)
#' 
#' @param radians  The radians to rotate. Positive for clock-wise rotation, and
#'   negative for counterclock-wise rotation.
#' @param is.rh    Whether the coordinate system is right-handed (by default
#'   TRUE).
#' @param is.clock Whether to rotate clock-wisely (by default TRUE).
#' 
#' @returns A 2-by-2 rotation matrix.
#'
#' @export
get_rot_mtx <- function(radians, is.rh = TRUE, is.clock = TRUE) {
  if (!is.clock) radians <- -radians
  if (!is.rh) radians <- -radians
  
  matrix(
    c(cos(radians), sin(radians), -sin(radians), cos(radians)),
    nrow = 2
  )
}
