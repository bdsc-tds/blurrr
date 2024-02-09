# `Blurrr` methods ---------------------------------------------

## getters ---------------------------------------------

#' @export
setGeneric(
  "get_template",
  function(x) standardGeneric("get_template")
)


#' @export
setGeneric(
  "get_molecule",
  function(x) standardGeneric("get_molecule")
)


#' @export
setGeneric(
  "get_is2Subspot",
  function(x) standardGeneric("get_is2Subspot")
)


#' @export
setGeneric(
  "get_mappedMolecule",
  function(x) standardGeneric("get_mappedMolecule")
)


#' @export
setGeneric(
  "get_spot_map",
  function(x) standardGeneric("get_spot_map")
)


#' @export
setGeneric(
  "get_spot_mapInDoubt",
  function(x) standardGeneric("get_spot_mapInDoubt")
)


#' @export
setGeneric(
  "get_spot_cnt",
  function(x) standardGeneric("get_spot_cnt")
)


#' @export
setGeneric(
  "get_mappedMolecule2subspot",
  function(x) standardGeneric("get_mappedMolecule2subspot")
)


#' @export
setGeneric(
  "get_subspot_map",
  function(x, idx) standardGeneric("get_subspot_map")
)


#' @export
setGeneric(
  "get_subspot_mapInDoubt",
  function(x) standardGeneric("get_subspot_mapInDoubt")
)


#' @export
setGeneric(
  "get_subspot_cnt",
  function(x) standardGeneric("get_subspot_cnt")
)


#' @export
setGeneric(
  "get_propAssigned",
  function(x) standardGeneric("get_propAssigned")
)


## others ---------------------------------------------

#' @export
setGeneric(
  "assign2template",
  function(x, cores, force = FALSE, verbose = TRUE) standardGeneric("assign2template")
)


#' @export
setGeneric(
  "subset2assigned",
  function(
    x,
    mole = c("cell", "trans"),
    level = c("spot", "subspot"),
    subspot.idx = seq_len(6)
  ) standardGeneric("subset2assigned")
)


#' @export
setGeneric(
  "plot_mole",
  function(
    x,
    mole = c("cell", "trans"),
    mode = c("raw", "assigned", "both"),
    level = c("spot", "subspot"),
    subspot.idx = seq_len(6),
    res = c("fullres", "hires", "lowres"),
    ...
  ) standardGeneric("plot_mole")
)


#' @export
setGeneric(
  "save2disk",
  function(x, dirname, cores = 1, cache = FALSE, verbose = TRUE) standardGeneric("save2disk")
)


# `Template` methods ---------------------------------------------

#' @export
setGeneric(
  "get_templateType",
  function(x) standardGeneric("get_templateType")
)

#' @export
setGeneric(
  "get_pos",
  function(x) standardGeneric("get_pos")
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


# `Molecule` methods ---------------------------------------------

#' @export
setGeneric(
  "get_moleculeType",
  function(x) standardGeneric("get_moleculeType")
)


#' @export
setGeneric(
  "get_sce",
  function(x) standardGeneric("get_sce")
)


# `MappedMolecule` methods ---------------------------------------------

#' @export
setGeneric(
  "get_map",
  function(x, idx) standardGeneric("get_map")
)


#' @export
setGeneric(
  "get_mapInDoubt",
  function(x) standardGeneric("get_mapInDoubt")
)


#' @export
setGeneric(
  "get_cnt",
  function(x) standardGeneric("get_cnt")
)


#' @export
setGeneric(
  "get_propAssigned",
  function(x) standardGeneric("get_propAssigned")
)
