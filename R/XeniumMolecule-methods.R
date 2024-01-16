# getters ---------------------------------------------

## For `XeniumMolecule` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xnPos",
  c("XeniumMolecule", "missing"),
  function(x, is.cell) x@pos
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xn_id",
  c("XeniumMolecule", "missing"),
  function(x, is.cell) get_xnPos(x)$transcript_id
)


## For `XeniumCell` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xn_id",
  c("XeniumCell", "missing"),
  function(x, is.cell) get_xnPos(x)$cell_id
)


## For `BinXenium` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xnPos",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_xnPos(get_xnCell(x)))
    } else {
      return(get_xnPos(get_xnTrans(x)))
    }
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xn_id",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_xn_id(get_xnCell(x)))
    } else {
      return(get_xn_id(get_xnTrans(x)))
    }
  }
)
