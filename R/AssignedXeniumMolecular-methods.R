# getters ---------------------------------------------

## For `AssignedXeniumMolecular` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignment2Spots",
  c("AssignedXeniumMolecular", "missing"),
  function(x, is.cell) x@assignment2Spots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_ambiAssignment2Spots",
  c("AssignedXeniumMolecular", "missing"),
  function(x, is.cell) x@ambiAssignment2Spots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_countOfSpots",
  c("AssignedXeniumMolecular", "missing"),
  function(x, is.cell) x@assignment2Spots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignment2Subspots",
  c("AssignedXeniumMolecular", "missing"),
  function(x, is.cell) x@assignment2Spots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_ambiAssignment2Subspots",
  c("AssignedXeniumMolecular", "missing"),
  function(x, is.cell) x@ambiAssignment2Spots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_countOfSubspots",
  c("AssignedXeniumMolecular", "missing"),
  function(x, is.cell) x@assignment2Spots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_propAssigned",
  c("AssignedXeniumMolecular", "missing"),
  function(x, is.cell) x@assignment2Spots
)


## For `BinXenium` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignment2Spots",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_assignment2Spots(get_xnCell(x)))
    } else {
      return(get_assignment2Spots(get_xnTrans(x)))
    }
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_ambiAssignment2Spots",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_ambiAssignment2Spots(get_xnCell(x)))
    } else {
      return(get_ambiAssignment2Spots(get_xnTrans(x)))
    }
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_countOfSpots",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_countOfSpots(get_xnCell(x)))
    } else {
      return(get_countOfSpots(get_xnTrans(x)))
    }
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignment2Subspots",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_assignment2Subspots(get_xnCell(x)))
    } else {
      return(get_assignment2Subspots(get_xnTrans(x)))
    }
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_ambiAssignment2Subspots",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_ambiAssignment2Subspots(get_xnCell(x)))
    } else {
      return(get_ambiAssignment2Subspots(get_xnTrans(x)))
    }
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_countOfSubspots",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_countOfSubspots(get_xnCell(x)))
    } else {
      return(get_countOfSubspots(get_xnTrans(x)))
    }
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_propAssigned",
  c("BinXenium", "logical"),
  function(x, is.cell) {
    if (is.cell) {
      return(get_propAssigned(get_xnCell(x)))
    } else {
      return(get_propAssigned(get_xnTrans(x)))
    }
  }
)
