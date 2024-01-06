#' @export
#'
#' @importFrom assertthat assert_that
get_coords_names <- function(
    is.xn,
    is.cell = TRUE,
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
