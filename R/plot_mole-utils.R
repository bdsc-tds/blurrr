#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 ggplot aes coord_fixed scale_x_continuous scale_y_reverse labs theme_bw
#' @importFrom ggforce geom_circle
.plot_bg <- function(vs.info, res = c("hires", "lowres", "fullres")) {
  assert_that(identical(class(vs.info)[[1]], "VisiumInfo"))
  
  res <- match.arg(res)
  
  .names <- get_coords_names(
    is.xn = FALSE,
    use.names = res
  )[[1]]
  
  .radius <- get_scalef(vs.info, "diameter") / 2
  if (res != "fullres") .radius <- .radius * get_scalef(vs.info, res)
  
  ggplot() +
    geom_circle(
      aes(
        x0 = !!as.name(.names[["x"]]),
        y0 = !!as.name(.names[["y"]]),
        r = .radius
      ),
      colour = "salmon",
      data = get_vsPos(vs.info)
    ) +
    coord_fixed() +
    scale_x_continuous(
      position = "top"
    ) +
    scale_y_reverse() +
    labs(
      x = "x (column)",
      y = "y (row)"
    ) +
    theme_bw()
}


#' @keywords internal
#' 
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 geom_point aes
.plot_mole <- function(mole.data, is.cell, res = c("hires", "lowres", "fullres")) {
  res <- match.arg(res)
  
  .names <- get_coords_names(
    is.xn = TRUE,
    is.cell = is.cell,
    use.names = res
  )[[1]]
  
  geom_point(
    aes(
      x = !!as.name(.names[["x"]]),
      y = !!as.name(.names[["y"]])
    ),
    size = 1e-3,
    colour = "steelblue",
    alpha = 0.2,
    data = mole.data
  )
}
