# getters ---------------------------------------------

## For `Template` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_templateType",
  c("Template"),
  function(x) x@type
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_pos",
  c("Template"),
  function(x) x@pos
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_arrayDirec",
  c("Template", "character"),
  function(x, key = c("row", "col", "both")) {
    key = match.arg(key)
    
    switch(
      key,
      "row" = {
        return(x@arrayDirec$row)
      },
      "col" = {
        return(x@arrayDirec$col)
      },
      {
        return(x@arrayDirec)
      }
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_imgRad",
  "Template",
  function(x) x@imgRad
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_scalef",
  c("Template", "character"),
  function(x, key = c("hires", "lowres", "diameter")) {
    key = match.arg(key)
    
    switch(
      key,
      "hires" = {
        return(x@scalef$tissue_hires_scalef)
      },
      "lowres" = {
        return(x@scalef$tissue_lowres_scalef)
      },
      "diameter" = {
        return(x@scalef$spot_diameter_fullres)
      },
      {
        return(NULL)
      }
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_subspotPos",
  "Template",
  function(x) x@subspotPos
)


## For `Blurrr` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_templateType",
  c("Blurrr"),
  function(x) get_templateType(get_template(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_pos",
  c("Blurrr"),
  function(x) get_pos(get_template(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_arrayDirec",
  c("Blurrr", "character"),
  function(x, key) get_arrayDirec(get_template(x), key)
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_imgRad",
  "Blurrr",
  function(x) get_imgRad(get_template(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_scalef",
  c("Blurrr", "character"),
  function(x, key) get_scalef(get_template(x), key)
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_subspotPos",
  "Blurrr",
  function(x) get_subspotPos(get_template(x))
)
