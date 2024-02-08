# getters ---------------------------------------------

## For `Molecule` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_moleculeType",
  c("Molecule"),
  function(x) x@type
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_sce",
  c("Molecule"),
  function(x) x@sce
)


## For `Blurrr` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_moleculeType",
  c("Blurrr"),
  function(x) get_moleculeType(get_molecule(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_sce",
  c("Blurrr"),
  function(x) get_sce(get_molecule(x))
)
