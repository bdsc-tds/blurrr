# `BinXenium` methods ---------------------------------------------

## getters ---------------------------------------------

#' @export
setGeneric(
  "get_vsInfo",
  function(x) standardGeneric("get_vsInfo")
)


#' @export
setGeneric(
  "get_xnCell",
  function(x) standardGeneric("get_xnCell")
)


#' @export
setGeneric(
  "get_xnTrans",
  function(x) standardGeneric("get_xnTrans")
)


#' @export
setGeneric(
  "get_is2Subspot",
  function(x) standardGeneric("get_is2Subspot")
)


#' @export
setGeneric(
  "get_assignedXnCell",
  function(x) standardGeneric("get_assignedXnCell")
)


#' @export
setGeneric(
  "get_assignedXnTrans",
  function(x) standardGeneric("get_assignedXnTrans")
)


## others ---------------------------------------------

#' @export
setGeneric(
  "assign2visium",
  function(x, cores, force = FALSE, verbose = TRUE) standardGeneric("assign2visium")
)


# `VisiumInfo` methods ---------------------------------------------

#' @export
setGeneric(
  "get_vsPos",
  function(x) standardGeneric("get_vsPos")
)

#' @export
setGeneric(
  "get_arrayDirec",
  function(x, key) standardGeneric("get_arrayDirec")
)


#' @export
setGeneric(
  "get_imgRad",
  function(x) standardGeneric("get_imgRad")
)


#' @export
setGeneric(
  "get_scalef",
  function(x, key) standardGeneric("get_scalef")
)


#' @export
setGeneric(
  "get_subspotPos",
  function(x) standardGeneric("get_subspotPos")
)


# `XeniumMolecular` methods ---------------------------------------------

#' @export
setGeneric(
  "get_xnPos",
  function(x, is.cell) standardGeneric("get_xnPos")
)


#' @export
setGeneric(
  "get_xn_id",
  function(x, is.cell) standardGeneric("get_xn_id")
)


# `AssignedXeniumMolecular` methods ---------------------------------------------

#' @export
setGeneric(
  "get_assignment2Spots",
  function(x, is.cell) standardGeneric("get_assignment2Spots")
)


#' @export
setGeneric(
  "get_ambiAssignment2Spots",
  function(x, is.cell) standardGeneric("get_ambiAssignment2Spots")
)


#' @export
setGeneric(
  "get_countOfSpots",
  function(x, is.cell) standardGeneric("get_countOfSpots")
)

#' @export
setGeneric(
  "get_assignment2Subspots",
  function(x, is.cell) standardGeneric("get_assignment2Subspots")
)

#' @export
setGeneric(
  "get_ambiAssignment2Subspots",
  function(x, is.cell) standardGeneric("get_ambiAssignment2Subspots")
)

#' @export
setGeneric(
  "get_countOfSubspots",
  function(x, is.cell) standardGeneric("get_countOfSubspots")
)

#' @export
setGeneric(
  "get_propAssigned",
  function(x, is.cell) standardGeneric("get_propAssigned")
)
