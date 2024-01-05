#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "vsInfo",
  "BinXenium",
  function(x) x@vsInfo
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "xnCell",
  "BinXenium",
  function(x) x@xnCell
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "xnTrans",
  "BinXenium",
  function(x) x@xnTrans
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "is2Spot",
  "BinXenium",
  function(x) x@is2Spot
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "is2Subspot",
  "BinXenium",
  function(x) x@is2Subspot
)
