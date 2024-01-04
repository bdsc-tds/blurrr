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
#' @export
get_rot_mtx <- function(radians, is.rh = TRUE, is.clock = TRUE) {
  if (!is.clock) radians <- -radians
  if (!is.rh) radians <- -radians
  
  matrix(
    c(cos(radians), sin(radians), -sin(radians), cos(radians)),
    nrow = 2
  )
}
