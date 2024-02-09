# getters ---------------------------------------------

## For `MappedMolecule` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate filter select
#' @importFrom stringr str_split_i
setMethod(
  "get_map",
  c("MappedMolecule", "numeric"),
  function(x, idx) {
    if (is.na(idx)) {
      ret <- x@map
    } else {
      assert_that(idx %in% seq_len(6))
      
      ret <- x@map %>%
        mutate(
          subspot_idx = str_split_i(
            subspot_barcode_bayesspace,
            ":",
            -1
          )
        ) %>%
        filter(
          subspot_idx == !!idx
        ) %>%
        select(
          -subspot_idx
        )
    }
    
    return(ret)
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_map",
  c("MappedMolecule", "missing"),
  function(x, idx = NA_integer_) get_map(x, idx = NA_integer_)
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_mapInDoubt",
  c("MappedMolecule"),
  function(x) x@mapInDoubt
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_cnt",
  c("MappedMolecule"),
  function(x) x@cnt
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_propAssigned",
  c("MappedMolecule"),
  function(x) x@propAssigned
)


## For `Blurrr` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_spot_map",
  c("Blurrr"),
  function(x) get_map(get_mappedMolecule(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_spot_mapInDoubt",
  c("Blurrr"),
  function(x) get_mapInDoubt(get_mappedMolecule(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_spot_cnt",
  c("Blurrr"),
  function(x) get_cnt(get_mappedMolecule(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_subspot_map",
  c("Blurrr", "numeric"),
  function(x, idx) get_map(get_mappedMolecule2subspot(x), idx)
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_subspot_mapInDoubt",
  c("Blurrr"),
  function(x) get_mapInDoubt(get_mappedMolecule2subspot(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_subspot_cnt",
  c("Blurrr"),
  function(x) get_cnt(get_mappedMolecule2subspot(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_propAssigned",
  c("Blurrr"),
  function(x) get_propAssigned(get_mappedMolecule(x))
)

