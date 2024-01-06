#' 
#' 
#' @keywords internal
#' 
#' @importFrom assertthat assert_that
.align_img2array <- function(
    data,
    radians,
    prev.cols,
    transformed.cols
) {
  assert_that(
    is.list(prev.cols),
    is.list(transformed.cols),
    all(names(prev.cols) %in% names(transformed.cols))
  )
  
  transformed.coords <- lapply(
    names(prev.cols),
    function(.type) {
      assert_that(
        length(prev.cols[[.type]]) == 2,
        length(transformed.cols[[.type]]) == 2,
        all(prev.cols[[.type]] %in% colnames(data)),
        all(c("x", "y") %in% names(prev.cols[[.type]])),
        all(c("x", "y") %in% names(transformed.cols[[.type]]))
      )
      
      .transformed.coords <- as.matrix(
        data[prev.cols[[.type]][c("x", "y")]]
      ) %*% t(
        get_rot_mtx(radians, is.rh = FALSE)
      )
      colnames(.transformed.coords) <- transformed.cols[[.type]][c("x", "y")]
      
      as.data.frame(.transformed.coords)
    }
  )
  
  do.call(
    cbind,
    c(data, transformed.coords)
  )
}
