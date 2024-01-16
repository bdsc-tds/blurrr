# getters ---------------------------------------------

## For `XeniumCell` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xn_sce",
  c("XeniumCell"),
  function(x) x@sce
)


## For `BinXenium` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xn_sce",
  c("BinXenium"),
  function(x) get_xn_sce(get_xnCell(x))
)
