#' @name Blurrr-methods
#' 
#' @title Methods for class `Blurrr`
#' 
#' @aliases
#' get_template
#' get_molecule
#' get_is2Subspot
#' get_mappedMolecule
#' get_mappedMolecule2subspot
#' assign2visium
#' subset2assigned
#' plot_mole
#' save2disk
#' 
#' 
NULL


# getters ---------------------------------------------

#' @rdname Blurrr-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_template",
  "Blurrr",
  function(x) x@template
)


#' @rdname Blurrr-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_molecule",
  "Blurrr",
  function(x) x@molecule
)


#' @rdname Blurrr-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_is2Subspot",
  "Blurrr",
  function(x) x@is2Subspot
)


#' @rdname Blurrr-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_mappedMolecule",
  "Blurrr",
  function(x) x@mappedMolecule
)


#' @rdname Blurrr-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_mappedMolecule2subspot",
  "Blurrr",
  function(x) x@mappedMolecule2subspot
)
