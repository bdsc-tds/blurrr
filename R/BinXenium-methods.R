# getters ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_vsInfo",
  "BinXenium",
  function(x) x@vsInfo
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xnCell",
  "BinXenium",
  function(x) x@xnCell
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xnTrans",
  "BinXenium",
  function(x) x@xnTrans
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_is2Subspot",
  "BinXenium",
  function(x) x@is2Subspot
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignedXnCell",
  "BinXenium",
  function(x) x@assignedXnCell
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignedXnTrans",
  "BinXenium",
  function(x) x@assignedXnTrans
)


# others ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "assign2visium",
  c("BinXenium", "missing"),
  function(x) assign2visium(x, force = FALSE)
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "assign2visium",
  c("BinXenium", "logical"),
  function(x, force) {
    if (!is.null(get_xnCell(x))) {
      x@assignedXnCell <- .assign2visium(
        assigned.xn = get_assignedXnCell(x),
        vs.pos = get_vsPos(x),
        xn.pos = get_xnPos(x, is.cell = TRUE),
        is2subspot = get_is2Subspot(x),
        subspot.pos = get_subspotPos(x),
        is.cell = TRUE,
        spot_radius = get_scalef(x, "d") / 2,
        force = force
      )
    }
    
    if (!is.null(get_xnTrans(x))) {
      x@assignedXnTrans <- .assign2visium(
        assigned.xn = get_assignedXnTrans(x),
        vs.pos = get_vsPos(x),
        xn.pos = get_xnPos(x, is.cell = FALSE),
        is2subspot = get_is2Subspot(x),
        subspot.pos = get_subspotPos(x),
        is.cell = FALSE,
        spot_radius = get_scalef(x, "d") / 2,
        force = force
      )
    }
    
    x
  }
)
