# getters ---------------------------------------------

## For `AssignedXeniumMolecule` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignment2Spots",
  c("AssignedXeniumMolecule", "missing"),
  function(x, is.cell) x@assignment2Spots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_ambiAssignment2Spots",
  c("AssignedXeniumMolecule", "missing"),
  function(x, is.cell) x@ambiAssignment2Spots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_countOfSpots",
  c("AssignedXeniumMolecule", "missing"),
  function(x, is.cell) x@countOfSpots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignment2Subspots",
  c("AssignedXeniumMolecule", "missing"),
  function(x, is.cell) x@assignment2Subspots
)

#' @include generic-def.R class-def.R
#' 
#' @export
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate filter select
#' @importFrom stringr str_split_i
setMethod(
  "get_assignment2Subspot",
  c("AssignedXeniumMolecule", "missing", "numeric"),
  function(x, is.cell, subspot.idx) {
    assert_that(subspot.idx %in% seq_len(6))
    
    x@assignment2Subspots %>%
      mutate(
        subspot_idx = str_split_i(
          subspot_barcode_bayesspace,
          ":",
          -1
        )
      ) %>%
      filter(
        subspot_idx == !!subspot.idx
      ) %>%
      select(
        -subspot_idx
      )
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_ambiAssignment2Subspots",
  c("AssignedXeniumMolecule", "missing"),
  function(x, is.cell) x@ambiAssignment2Subspots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_countOfSubspots",
  c("AssignedXeniumMolecule", "missing"),
  function(x, is.cell) x@countOfSubspots
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_propAssigned",
  c("AssignedXeniumMolecule", "missing"),
  function(x, is.cell) x@propAssigned
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
      return(get_assignment2Spots(get_assignedXnCell(x)))
    } else {
      return(get_assignment2Spots(get_assignedXnTrans(x)))
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
      return(get_ambiAssignment2Spots(get_assignedXnCell(x)))
    } else {
      return(get_ambiAssignment2Spots(get_assignedXnTrans(x)))
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
      return(get_countOfSpots(get_assignedXnCell(x)))
    } else {
      return(get_countOfSpots(get_assignedXnTrans(x)))
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
      return(get_assignment2Subspots(get_assignedXnCell(x)))
    } else {
      return(get_assignment2Subspots(get_assignedXnTrans(x)))
    }
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignment2Subspot",
  c("BinXenium", "logical", "numeric"),
  function(x, is.cell, subspot.idx) {
    if (is.cell) {
      return(get_assignment2Subspot(
        get_assignedXnCell(x),
        subspot.idx = subspot.idx
      ))
    } else {
      return(get_assignment2Subspot(
        get_assignedXnTrans(x),
        subspot.idx = subspot.idx
      ))
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
      return(get_ambiAssignment2Subspots(get_assignedXnCell(x)))
    } else {
      return(get_ambiAssignment2Subspots(get_assignedXnTrans(x)))
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
      return(get_countOfSubspots(get_assignedXnCell(x)))
    } else {
      return(get_countOfSubspots(get_assignedXnTrans(x)))
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
      return(get_propAssigned(get_assignedXnCell(x)))
    } else {
      return(get_propAssigned(get_assignedXnTrans(x)))
    }
  }
)
